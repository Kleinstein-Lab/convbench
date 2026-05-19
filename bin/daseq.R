#!/usr/bin/env Rscript
message(paste0('Starting run: ', Sys.time()))

suppressPackageStartupMessages({
  library(argparse)
  library(DAseq)
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(RColorBrewer)
  library(scales)
  library(cowplot)
  library(patchwork)
  library(ggrastr)
  library(pracma)
  library(uwot)
  require(R.utils)
  # library(alakazam)
})

################################################################################
# options that will be hard-coded for now

set.seed(37)
options(future.globals.maxSize = 16 * 1024^3)

################################################################################
# user input parameters

# Create a parser object
# need to provide inputs for:
#############################
# python location
# data location
# metadata location
# output location
# k hyperparameters
# V(D)J info present or not
# single cell V(D)J info provided or not
# simulated data vs. real data
#############################

parser <- ArgumentParser(description = "Data location and DAseq algorithm hyperparameters.")

# NOTE: as currently written, assumes that data and metadata have matching cell_id or sequence_id columns
# and that each cell ID is unique -- i.e. no repeated cell IDs will be found in different
# subjects. If this is a possibility, unique cell IDs that map 1:1 between data and metadata
# should be created before running DA tools, which may require some pre-processing.
# A subject_id and sample_id column are also both assumed. If no sample id provided, will 
# automatically copy the subject_id into a sample_id column.
parser$add_argument('-d', '--data_loc', type = 'character', default = 'data',
                    help = 'File path for the embedding or RNA-Seq data location.')

parser$add_argument('-md', '--metadata_loc', type = 'character', default = 'metadata',
                    help = 'File path for the metadata location. Metadata and data files should have 1:1 matching sequence identifiers.')

parser$add_argument('-o', '--output_dir', type = 'character', default = 'DAseq_output',
                    help = 'Specify an output directory location.')

parser$add_argument('-da', '--da_variable', type = 'character', default = 'status',
                    help = 'Stratification variable that should be used to determine for differential abundance. There should be two levels in this factor/categorical variable.')

parser$add_argument('-m', '--k_min', type = 'integer', default = 50,
                    help = 'Minimum number of neighbors to use in KNN algorithm.')

parser$add_argument('-t', '--k_step', type = 'integer', default = 50,
                    help = 'Step value for KNN neighbor values.')

parser$add_argument('-x', '--k_max', type = 'integer', default = 300,
                    help = 'Maximum number of neighbors to use in KNN algorithm.')

parser$add_argument('-v', '--vdj_info', type = 'logical', default = TRUE,
                    help = 'Is v call and j call information included in the metadata? Can apply to expression or embedding data.')

# TODO: can change to be more granular/option to plot at gene, family etc. level
# right now defaults to v_call and j_call columns and removes allele info
parser$add_argument('-sc', '--single_cell', type = 'logical', default = FALSE,
                    help = 'Input true if V(D)J info is present and contains paired heavy and light chain info.')

parser$add_argument('-si', '--simulated', type = 'logical', default = FALSE,
                    help = 'Specify whether input data is simulated or real.')

parser$add_argument('-r', '--remove_dups', type = 'logical', default = FALSE,
                    help = 'Will remove duplicate embeddings within an individual if TRUE.')

parser$add_argument('-w', '--overwrite', type = 'logical', default = FALSE,
                    help = 'Specify whether to recalculate UMAP and da_cells.')

################################################################################

# Parse the arguments
args <- parser$parse_args()

# specify which dataset we are analyzing
DATA_LOC <- args$data_loc
MD_LOC <- args$metadata_loc

OUTPUT_DIR <- args$output_dir

message(paste0('Data will be saved to ', OUTPUT_DIR, '.'))

# get the variable for DA calculations
DA_VAR <- args$da_variable

KVEC <- seq(args$k_min, args$k_max, args$k_step)

message(paste0('K nearest neighbor values: ', paste0(KVEC, collapse = ' ')))

VDJ <- args$vdj_info
SINGLE_CELL <- args$single_cell
SIMULATED <- args$simulated

OVERWRITE <- args$overwrite
REMOVE_DUPS <- args$remove_dups

################################################################################

if (VDJ){
  message('V(D)J calls included in metadata.')
} else{
  message('V(D)J calls not provided.')
}

if (SINGLE_CELL){
  message('Paired heavy and light chain info provided.')
} else{
  message('Bulk V(D)J info only available.')
}

if (SIMULATED){
  message('Simulated data present.')
} else{
  message('Simulated data not present.')
}

if (REMOVE_DUPS){
  message('Duplicate embeddings within a subject will be collapsed.')
}

################################################################################
# create locations for figures and results to be saved within output dir

if(!dir.exists(file.path(OUTPUT_DIR))){
  dir.create(file.path(OUTPUT_DIR))
}

if(!dir.exists(file.path(OUTPUT_DIR, 'figures'))){
  dir.create(file.path(OUTPUT_DIR, 'figures'))
}

if(!dir.exists(file.path(OUTPUT_DIR, 'tables'))){
  dir.create(file.path(OUTPUT_DIR, 'tables'))
}

################################################################################
#################
### LOAD DATA ###
#################

# load embeddings or expr data
message(paste0('Loading data: ', DATA_LOC))

tryCatch(
  
  {
    data <- data.table::fread(DATA_LOC, sep = '\t', header = T)
  }, error = function(e){
    
    stop(e)
    
  }
)

# metadata
message(paste0('Loading metadata: ', MD_LOC))

tryCatch(
  
  {
    md <- readr::read_tsv(MD_LOC)
  }, error = function(e){
    
    stop(e)
    
  }
  
)

# standardize column names
colnames(md) <- tolower(colnames(md))

# create artificial sample_id copies from subject ID if not present
if (!'sample_id' %in% colnames(md)){
  md$sample_id <- md$subject_id
}

# change to a generic id column
if ('sequence_id' %in% colnames(md) & 'sequence_id' %in% colnames(data)){
  
  ID_COL_NAME <- 'sequence_id'
  names(md)[names(md) == 'sequence_id'] <- 'id_col'
  names(data)[names(data) == 'sequence_id'] <- 'id_col'
  
} else if ('cell_id' %in% colnames(md) & 'cell_id' %in% colnames(data)){
  
  ID_COL_NAME <- 'cell_id'
  names(md)[names(md) == 'cell_id'] <- 'id_col'
  names(data)[names(data) == 'cell_id'] <- 'id_col'
  
} else{
  
  stop('Matching cell_id or sequence_id columns not found in data and metadata files.')
  
}

if (REMOVE_DUPS){
  
  old_seq_num <- nrow(data)
  
  # add subject info
  data <- data %>%
    dplyr::left_join(md[c('id_col', 'subject_id')], by = 'id_col')
  
  # get distinct sequences within individuals
  row.names(data) <- data$id_col
  data <- data %>% dplyr::select(-id_col)
  data <- distinct(data)
  data <- data %>% dplyr::select(-subject_id)
  data$id_col <- row.names(data) # add back in for now
  
  new_seq_num <- nrow(data)
  
  seqs_removed <- old_seq_num - new_seq_num
  message(paste0('Duplicates removed. ', seqs_removed, ' sequences removed. New total: ', new_seq_num))
  
}

# not all of the seqs in the data will necessarily result in successful embeddings
# so we can filter the metadata for only the relevant cell info
md <- md %>%
  dplyr::filter(id_col %in% data$id_col) %>%
  as.data.frame()

row.names(md) <- md$id_col
row.names(data) <- data$id_col

# make sure data and md aligned properly
md <- md[row.names(data),]

################################################################################
#################
### PREP DATA ###
#################

# get label info
X.label.info <- md %>%
  dplyr::select(sample_id, !!sym(DA_VAR)) %>%
  distinct()

colnames(X.label.info) <- c('label', 'condition')

# ensure the condition is a factor
X.label.info$condition <- factor(X.label.info$condition)

# throw an error if there aren't two levels
if (length(levels(X.label.info$condition)) != 2){
  stop('Provided DA variable is not a factor with 2 levels.')
}

# get labels for both conditions
label_gps <- lapply(levels(X.label.info$condition), function(current_cond){
  label_choice_cells <- X.label.info %>% dplyr::filter(condition == current_cond) %>%
    dplyr::pull(label)
  return(label_choice_cells)
})

# now, prep the embedding data
# made a guide of cell labels to sample labels
X.cells <- md %>%
  dplyr::select(subject_id, sample_id, id_col) %>%
  distinct() %>%
  data.frame()

print(table(X.cells$subject_id))

# put cell ids in the rownames
row.names(X.cells) <- X.cells$id_col
row.names(data) <- data$id_col

# should not have this column for making UMAP
data <- data %>%
  dplyr::select(-c(id_col))

# check that IDs are unique
if (nrow(data) != length(unique(row.names(data)))){
  warning('Beware: Not all cell IDs in provided data are unique.')
}

# narrow to the common cells and line up the cell IDs by matching X.cells
# rownames to data rownames
X.cells <- X.cells[row.names(data),]

# retrieve sample labels - will need these to run DAseq later!
X.label.embeddings <- X.cells$sample_id

# make a UMAP (switched from tSNE)
Sys.time()
message('Data loaded and prepped for DAseq. Generating PCA & UMAP...')
# tSNE_embeddings <- Rtsne(data,
#                          check_duplicates = FALSE) # there will be duplicated embed. so just turn off

# try with python 3.7? Or use virtual environment instead of conda if that doesn't work
# NOT WORKING - next step try virtual env instead, and if that doesn't work try uwot version
# reticulate::use_condaenv("daseq_v5", required = TRUE) 
# reticulate::use_virtualenv('~/project/convergence/tools/DAseq/.venv', required=TRUE)
# umap_embeddings <- umap(data, method = 'umap-learn')

# run PCA
if (!file.exists(file.path(OUTPUT_DIR, 'tables', 'UMAP_embeddings.rds')) | OVERWRITE == T){
  
  if (nrow(data) >= 200){
    pca <- prcomp(data, center = T, scale. = T)
    umap_embeddings <- uwot::umap(pca$x[, 1:200]) # use 200 PCs
  } else{
    umap_embeddings <- uwot::umap(data)
  }
  
  saveRDS(umap_embeddings, file.path(OUTPUT_DIR, 'tables', 'UMAP_embeddings.rds'))
  
} else{
  cat(paste0('Loading embeddings from: ', file.path(OUTPUT_DIR, 'tables', 'UMAP_embeddings.rds')))
  umap_embeddings <- readRDS(file.path(OUTPUT_DIR, 'tables', 'UMAP_embeddings.rds'))
}

################################################################################
################
### DATA VIZ ###
################

Sys.time()
message('UMAP embeddings generated. Making visualizations...')

make_umap_viz <- function(var, var_name, custom_pal = NULL){
  
  default_h <- 12
  default_w <- 14
  
  # let's get some visualizations first
  p <- ggplot(umap_embeddings_coords, aes(x = UMAP1, y = UMAP2, color = !!sym(var))) +
    geom_point(alpha = 0.7, size = 0.8) +
    theme_minimal(base_size = 15) +
    labs(
      x = "UMAP 1",
      y = "UMAP 2",
      color = var_name
    ) +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    theme_cowplot()
  
  num_vars <- length(unique(umap_embeddings_coords[[var]]))
  
  if (num_vars <= 9){
    
    if(!is.null(custom_pal)){
      p + scale_color_manual(values = custom_pal)
      
    } else{
      p + scale_color_brewer(palette = "Set1") 
    }
    
    ggsave(file.path(OUTPUT_DIR, 'figures', paste0('UMAP_', var, '.png')), 
           device = 'png', width = default_w, height = default_h, units = 'in')
  } else{
    
    if(!is.null(custom_pal)){
      p + scale_color_manual(values = custom_pal)
    }
    
    # scale width according to how big the legend is going to be
    long_w <- default_w + (0.5 * num_vars%/%18)
    ggsave(file.path(OUTPUT_DIR, 'figures', paste0('UMAP_', var, '.png')), 
           device = 'png', width = long_w, height = default_h, units = 'in')
    
  }
}

# prep UMAP data for vizualization
umap_embeddings_coords <- data.frame(umap_embeddings)
colnames(umap_embeddings_coords) <- c('UMAP1', 'UMAP2')
row.names(umap_embeddings_coords) <- row.names(data)
umap_embeddings_coords$id_col <- row.names(umap_embeddings_coords)
umap_embeddings_coords$sample_id <- X.cells$sample_id
umap_embeddings_coords$subject_id <- X.cells$subject_id

# add label information
umap_embeddings_coords <- umap_embeddings_coords %>%
  dplyr::left_join(X.label.info, by=join_by(sample_id == label))

# add other metadata for plotting
if (VDJ){
  
  if(!('v_gene' %in% colnames(md) & 'j_gene' %in% colnames(md))) warning('v_gene and j_gene columns not provided. UMAP plots for V and J gene will not be generated.')
  
  else {
    
    if (!SINGLE_CELL){
      
      umap_embeddings_coords <- umap_embeddings_coords %>%
        dplyr::left_join(md[c('v_gene', 'j_gene', 'id_col')], by = 'id_col')
      
      # get rid of alleles for v and j call
      # umap_embeddings_coords$v_gene <- getGene(umap_embeddings_coords$v_call)
      # umap_embeddings_coords$j_gene <- getGene(umap_embeddings_coords$j_call)
      
      make_umap_viz('v_gene', 'V gene')
      make_umap_viz('j_gene', 'J gene')
      
    } else if (SINGLE_CELL){
      
      ##### APPLIES TO SINGLE CELL ONLY #####
      # get heavy and light chain V/J assignments
      heavy_info <- md %>%
        dplyr::filter(id_col %in% row.names(data)) %>%
        dplyr::filter(locus == 'IGH') %>%
        dplyr::select(id_col, v_gene, j_gene) %>%
        distinct() %>%
        data.frame(check.names = F)
      
      row.names(heavy_info) <- heavy_info$id_col
      
      light_info <- md %>%
        dplyr::filter(id_col %in% row.names(data)) %>%
        dplyr::filter(locus == 'IGK' | locus == 'IGL') %>%
        dplyr::group_by(id_col) %>%
        dplyr::arrange(desc(consensus_count)) %>%
        dplyr::slice_head(n = 1) %>%
        dplyr::ungroup() %>%
        dplyr::select(id_col, v_gene, j_gene) %>%
        distinct() %>%
        data.frame(check.names = F)
      
      row.names(light_info) <- light_info$id_col
      
      umap_embeddings_coords <- umap_embeddings_coords %>%
        dplyr::left_join(heavy_info, by = 'id_col') %>%
        dplyr::rename(v_gene_heavy = v_gene,
                      j_gene_heavy = j_gene)
      
      umap_embeddings_coords <- umap_embeddings_coords %>%
        dplyr::left_join(light_info, by = 'id_col') %>%
        dplyr::rename(v_gene_light = v_gene,
                      j_gene_light = j_gene)
      
      # umap_coords$v_gene_heavy <- getGene(umap_coords$v_call_heavy)
      # umap_coords$j_gene_heavy <- getGene(umap_coords$j_call_heavy)
      # umap_coords$v_gene_light <- getGene(umap_coords$v_call_light)
      # umap_coords$j_gene_light <- getGene(umap_coords$j_call_light)
      
      ##### APPLIES TO SINGLE CELL ONLY #####
      make_umap_viz('v_gene_heavy', 'V Gene - \nHeavy Chain')
      make_umap_viz('j_gene_heavy', 'J Gene - \nHeavy Chain')
      make_umap_viz('v_gene_light', 'V Gene - \nLight Chain')
      make_umap_viz('j_gene_light', 'J Gene - \nLight Chain')
      
    }
  }
}

# include info if simulated
if (SIMULATED){
  umap_embeddings_coords <- umap_embeddings_coords %>%
    dplyr::left_join(md[c('id_col', 'simulated')], by = 'id_col')
  
  make_umap_viz('simulated', 'Simulated', custom_pal = c('TRUE' = "red", 'FALSE' = "gray"))
}

make_umap_viz('condition', DA_VAR)
make_umap_viz('sample_id', 'Sample ID')
make_umap_viz('subject_id', 'Subject ID')

# matrix version for downstream steps
X.embed <- as.matrix(umap_embeddings_coords %>% dplyr::select(c('UMAP1', 'UMAP2')))

#################################################################################
#############
### DASEQ ###
#############

# NOTE: change the k values if they are too small
# if (nrow(data) <= 100){
#   new_max_k <- min(25, round(nrow(data)/2))
#   KVEC <- seq(5, new_max_k, 2)
#   warning(paste0('Fewer than 100 cells - forcing K values to ', paste(KVEC, collapse = ', ')))
# } else if(nrow(data) > 100 & nrow(data) <= 1000){
#   new_max_k <- min(200, round(nrow(data)/4))
#   KVEC <- seq(5, new_max_k, 20)
#   warning(paste0('Fewer than 1,000 cells - forcing K values to ', paste(KVEC, collapse = ', ')))
# }

Sys.time()
message('Finding DA cells...')

# measure how long the DASeq process itself takes
start_time <- Sys.time()

# save memory for now
if (!file.exists(file.path(OUTPUT_DIR, 'tables', 'da_cells.rds')) | OVERWRITE == T){
  da_cells <- DAseq::getDAcells(
    X = data,
    cell.labels = X.label.embeddings,
    labels.1 = label_gps[[1]],
    labels.2 = label_gps[[2]],
    k.vector = KVEC, # can tweak this
    plot.embedding = X.embed
  )
  
  # saveRDS(da_cells, file.path(OUTPUT_DIR, 'tables', 'da_cells.rds'))
  
  da_cells$pred.plot
  ggsave(file.path(OUTPUT_DIR, 'figures', 'UMAP_pred_plot.png'), 
         device = 'png', width = 8, height = 7, units = 'in')
  
  da_cells$rand.plot
  
  da_cells$da.cells.plot
  ggsave(file.path(OUTPUT_DIR, 'figures', 'UMAP_da_cells.png'), 
         device = 'png', width = 8, height = 7, units = 'in')
  
  
} else{
  cat(paste0('Loading da_cells from: ', file.path(OUTPUT_DIR, 'tables', 'da_cells.rds')))
  da_cells <- readRDS(file.path(OUTPUT_DIR, 'tables', 'da_cells.rds'))
}

Sys.time()
message('Clustering DA regions...')

# clustering
# account for scenario in which no da regions are found
tryCatch({
  
  da_regions <- DAseq::getDAregion(
    X = data,
    da.cells = da_cells,
    cell.labels = X.label.embeddings,
    labels.1 = label_gps[[1]],
    labels.2 = label_gps[[2]],
    resolution = 0.01,
    plot.embedding = X.embed,
  )
  
}, error = function(e){
  
  warning('No DA regions found. Run ended without completing final analysis.')
  message(paste0('Finishing run: ', Sys.time()))
  print(sessionInfo())
  
  stop(e)
  
})

da_regions$da.region.plot
ggsave(file.path(OUTPUT_DIR, 'figures', 'UMAP_da_regions.png'), 
       device = 'png', width = 10, height = 8, units = 'in')

# get ending time after getting DA Regions and making basic figures
end_time <- Sys.time()
time_taken <- end_time - start_time

# preserve the original time taken if NOT re-creating DA cells and have a run table
if (OVERWRITE == F & file.exists(file.path(OUTPUT_DIR, 'tables', 'run_stats.tsv'))){
  run_stat_existing <- read.csv(file.path(OUTPUT_DIR, 'tables', 'run_stats.tsv'), sep = '\t', check.names = F)
  time_taken <- run_stat_existing$`time (min)`
}

################################################################################
X.cells$da.region.label <- da_regions$da.region.label

# add embedding and simulated info if applicable
if (VDJ){
  
  if (!SINGLE_CELL & 'v_gene' %in% colnames(umap_embeddings_coords)){ # assume if v gene not available, j gene not either
    
    X.cells <- X.cells %>%
      dplyr::left_join(umap_embeddings_coords[c('v_gene', 'j_gene', 'id_col')], by = 'id_col')
    
  } else if (SINGLE_CELL & 'v_gene_heavy' %in% colnames(umap_embeddings_coords)){ # assume if v gene heavy not available, all others are not
    
    X.cells <- X.cells %>%
      dplyr::left_join(umap_embeddings_coords[c('v_gene_heavy', 'j_gene_heavy', 'v_gene_light', 'j_gene_light', 'id_col')], by = 'id_col')
    
  }
}

if (SIMULATED){
  X.cells <- X.cells %>%
    dplyr::left_join(md[c('simulated', 'id_col')], by = 'id_col')
}

# add da var info
X.cells <- X.cells %>%
  dplyr::left_join(md[c(DA_VAR, 'id_col')], by = 'id_col')

# cluster-level stats
write.table(da_regions[["DA.stat"]], 
            file.path(OUTPUT_DIR, 'tables', 'region_stats.tsv'), 
            sep='\t', quote = F, row.names = F)

################################################################################
##################
### EVALUATION ###
##################

Sys.time()
message('Final data analysis...')

# prep the DA region stats for later use
DA.stat <- data.frame(da_regions[["DA.stat"]])
DA.stat$da.region.label <- row.names(DA.stat)

DA.stat %>%
  ggplot(aes(x = pval.ttest)) + 
  geom_histogram(color = 'white', binwidth = 0.01) + 
  theme_bw() +
  labs(title = 'DA-Seq T-test P-Value Distribution') +
  coord_cartesian(xlim = c(0, 1))

ggsave(file.path(OUTPUT_DIR, 'figures', 'pvalue_hist.png'),
       device = 'png', width = 8, height = 6, units = 'in')

# figure out which clusters are significant - nominal and adjust p val
DA.stat$wilcox.adj.BH <- p.adjust(DA.stat$pval.wilcoxon, method = 'BH')
DA.stat$ttest.adj.BH <- p.adjust(DA.stat$pval.ttest, method = 'BH')

# will need info about statistical tests later as well
tests <- c('pval.wilcoxon', 'pval.ttest', 'wilcox.adj.BH', 'ttest.adj.BH')

# add pred info to X.cells
X.cells$pred <- da_cells$da.pred

# add the stats
X.cells$da.region.label <- as.character(X.cells$da.region.label)
X.cells <- X.cells %>%
  dplyr::left_join(DA.stat, by = 'da.region.label')

# cell-level info
# restore ID column for writing
names(X.cells)[names(X.cells) == 'id_col'] <- ID_COL_NAME

write.table(X.cells, 
            file.path(OUTPUT_DIR, 'tables', 'da_seqs.tsv'), 
            sep='\t', quote = F, row.names = F)

names(X.cells)[names(X.cells) == ID_COL_NAME] <- 'id_col'

#######
# AUC #
#######

if (SIMULATED & 'simulated' %in% colnames(X.cells)){ 
  
  # make the AUC plot
  
  # get thresholds 
  # min_bg <- min(unlist(da_cells$rand.pred))
  # max_bg <- max(unlist(da_cells$rand.pred))
  # 
  # min_pred <- min(unlist(da_cells$da.pred))
  # max_pred <- max(unlist(da_cells$da.pred))
  
  da_score_auc_thresholds <- sort(unique(abs(da_cells$da.pred)))
  # da_score_auc_thresholds <- quantile(abs(da_cells$da.pred), seq(0, 1, 0.01), names=F)
  # auc_thresholds[1] <- auc_thresholds[1] - 1e-8
  tot_da_thresh <- length(da_score_auc_thresholds)
  da_score_auc_thresholds[tot_da_thresh] <- da_score_auc_thresholds[tot_da_thresh] + 1e-3
  
  da_score_auc_data <- lapply(da_score_auc_thresholds, function(thresh){
    
    DA_cells <- X.cells %>%
      dplyr::filter(abs(pred) >= thresh)
    
    non_DA_cells <- X.cells %>%
      dplyr::filter(abs(pred) < thresh)
    
    true_pos <- sum(DA_cells$simulated == TRUE)
    false_neg <- sum(non_DA_cells$simulated == TRUE)
    true_neg <- sum(non_DA_cells$simulated == FALSE)
    false_pos <- sum(DA_cells$simulated == FALSE)
    
    return(data.frame('TPR' = true_pos / (true_pos + false_neg),
                      'FPR' = 1 - (true_neg / (true_neg + false_pos))))
    
  })
  
  da_score_auc_df <- do.call(rbind, da_score_auc_data)
  da_score_auc_df$da_score_threshold <- da_score_auc_thresholds
  
  write.table(da_score_auc_df, 
              file.path(OUTPUT_DIR, 'tables', 'da_score_auc_curve_vals.tsv'), 
              sep = '\t', row.names = F, quote = F)
  
  # get auroc
  da_auroc <- pracma::trapz(rev(da_score_auc_df$FPR), rev(da_score_auc_df$TPR))
  
  da_score_auc_df %>%
    ggplot(aes(x = FPR, y = TPR)) +
    geom_point() +
    geom_line() +
    labs(x = 'FPR',
         y = 'TPR',
         title = paste0('DA Threshold ', round(min(da_score_auc_thresholds)), ' to ', round(max(da_score_auc_thresholds), 3)),
         subtitle = paste0('AUC: ', round(da_auroc, 3))) + 
    theme_minimal()
  
  ggsave(file.path(OUTPUT_DIR, 'figures', 'da_score_AUC_curve.png'),
         device = 'png',
         width = 7,
         height = 6)
  
  ##############################
  # ALTERNATIVE AUC CALCUATION #
  ##############################
  
  total_cells <- nrow(X.cells)
  unassigned_cells <- sum(is.na(X.cells$ttest.adj.BH))
  assigned_cells <- total_cells - unassigned_cells
  X.cells[is.na(X.cells$ttest.adj.BH), 'ttest.adj.BH'] <- 1
  
  p_auc_thresholds <- sort(unique(X.cells$ttest.adj.BH))
  # p_auc_thresholds <- quantile(X.cells$ttest.adj.BH, seq(0, 1, 0.01), names=F)
  # auc_thresholds[1] <- auc_thresholds[1] - 1e-8
  
  # add to the largest to make sure the entire curve is captured
  tot_p_thresh <- length(p_auc_thresholds)
  p_auc_thresholds[tot_p_thresh] <- p_auc_thresholds[tot_p_thresh] + 1e-3
  
  p_auc_data <- lapply(p_auc_thresholds, function(thresh){
    
    DA_cells <- X.cells %>%
      dplyr::filter(ttest.adj.BH < thresh)
    
    non_DA_cells <- X.cells %>%
      dplyr::filter(ttest.adj.BH >= thresh)
    
    true_pos <- sum(DA_cells$simulated == TRUE)
    false_neg <- sum(non_DA_cells$simulated == TRUE)
    true_neg <- sum(non_DA_cells$simulated == FALSE)
    false_pos <- sum(DA_cells$simulated == FALSE)
    
    return(data.frame('TPR' = true_pos / (true_pos + false_neg),
                      'FPR' = 1 - (true_neg / (true_neg + false_pos))))
    
  })
  
  p_auc_df <- do.call(rbind, p_auc_data)
  p_auc_df$da_score_threshold <- p_auc_thresholds
  
  write.table(p_auc_df, 
              file.path(OUTPUT_DIR, 'tables', 'p_auc_curve_vals.tsv'), 
              sep = '\t', row.names = F, quote = F)
  
  # get auroc
  p_auroc <- pracma::trapz(p_auc_df$FPR, p_auc_df$TPR)
  
  p_auc_df %>%
    ggplot(aes(x = FPR, y = TPR)) +
    geom_point() +
    geom_line() +
    labs(x = 'FPR',
         y = 'TPR',
         title = paste0('T-test FDR threshold ', round(min(p_auc_thresholds)), ' to ', round(max(p_auc_thresholds), 3)),
         subtitle = paste0('AUC: ', round(p_auroc, 3), '; ', 
                           prettyNum(assigned_cells, big.mark = ",", scientific = FALSE), '/', 
                           prettyNum(total_cells, big.mark = ",", scientific = FALSE), ' cells in DA clusters')) + 
    theme_minimal()
  
  ggsave(file.path(OUTPUT_DIR, 'figures', 'p_AUC_curve.png'),
         device = 'png',
         width = 7,
         height = 6)
  
  ###########
  # JACCARD #
  ###########
  jaccard_df <- X.cells
  
  jaccard_df <- jaccard_df %>%
    dplyr::mutate(p_under_0.005 = ttest.adj.BH <= 0.005,
                  p_under_0.05 = ttest.adj.BH <= 0.05,
                  p_under_0.1 = ttest.adj.BH <= 0.1)
  
  # calc jaccard index
  jaccard_005 <- sum(jaccard_df$simulated & jaccard_df$p_under_0.005, na.rm = T) / sum(jaccard_df$simulated | jaccard_df$p_under_0.005, na.rm = T)
  jaccard_05 <- sum(jaccard_df$simulated & jaccard_df$p_under_0.05, na.rm = T) / sum(jaccard_df$simulated | jaccard_df$p_under_0.05, na.rm = T)
  jaccard_1 <- sum(jaccard_df$simulated & jaccard_df$p_under_0.1, na.rm = T) / sum(jaccard_df$simulated | jaccard_df$p_under_0.1, na.rm = T)
  
  # get Jaccard across a range
  # jaccard_thresholds <- quantile(jaccard_df$ttest.adj.BH, seq(0, 1, 0.01), names=F, na.rm = T)
  jaccard_thresholds <- sort(unique(jaccard_df$ttest.adj.BH))
  jaccard_thresholds <- jaccard_thresholds[!is.na(jaccard_thresholds)]
  
  # auc_thresholds[1] <- auc_thresholds[1] - 1e-8
  # jaccard_thresholds[101] <- jaccard_thresholds[101] + 1e-8
  
  jaccards <- sapply(jaccard_thresholds, function(thresh){
    j <- sum(jaccard_df$simulated & jaccard_df$ttest.adj.BH <= thresh, na.rm = T) / sum(jaccard_df$simulated | jaccard_df$ttest.adj.BH <= thresh, na.rm = T)
  })
  
  # get max Jaccard and its corresponding p-value
  Jaccard_max <- max(jaccards, na.rm = T)
  Jaccard_max_p <- jaccard_thresholds[which.max(jaccards)]
  
  jaccard_plot_df <- data.frame('Adjusted P-Value Threshold' = jaccard_thresholds,
                                'Jaccard Similarity Index' = jaccards,
                                check.names = F)
  
  write.table(jaccard_plot_df, 
              file.path(OUTPUT_DIR, 'tables', 'jaccard_plot_vals.tsv'), 
              sep = '\t', row.names = F, quote = F)
  
  jaccard_plot <- jaccard_plot_df %>%
    ggplot(aes(x = !!sym('Adjusted P-Value Threshold'), y = !!sym('Jaccard Similarity Index'))) +
    geom_point() +
    geom_line() +
    theme_bw() +
    labs(title = 'Jaccard Similarity Across Adjusted P Thresholds',
         subtitle = paste0('Max Jaccard: ', round(Jaccard_max, 3), 
                           ' at adjusted P-value ', round(Jaccard_max_p, 3)))
  
  ggsave(filename = file.path(OUTPUT_DIR, 'figures', 'jaccard_plot.png'),
         plot = jaccard_plot,
         device = 'png',
         width = 7,
         height = 5)
  
} else{
  p_auroc <- NA
}

################################################################################
###################
## DATA SUMMARY ###
###################

# where do the sequences come from i.e. are they dominated by one patient?
subj_cts <- X.cells %>%
  dplyr::select('subject_id', 'da.region.label') %>%
  dplyr::group_by(da.region.label, subject_id) %>%
  dplyr::summarise(seqs_per_subj = n()) %>%
  dplyr::filter(da.region.label != 0)

cluster_cts <- X.cells %>%
  dplyr::group_by(da.region.label) %>%
  dplyr::summarise(seqs_per_cluster = n())

subj_cts <- dplyr::left_join(subj_cts, cluster_cts, by = 'da.region.label')

subj_cts <- subj_cts %>%
  dplyr::mutate(pct_subj = seqs_per_subj / seqs_per_cluster)

subj_cts$da.region.label <- as.character(subj_cts$da.region.label)

if (SIMULATED & 'simulated' %in% colnames(X.cells)){
  # also add what percentage of the cluster is simulated sequences
  sim_cts <- X.cells %>%
    dplyr::select('da.region.label', 'simulated') %>%
    dplyr::group_by(da.region.label) %>%
    dplyr::summarise(sim_seqs = sum(simulated == TRUE)) %>%
    dplyr::filter(da.region.label != 0)
  
  sim_cts$da.region.label <- as.character(sim_cts$da.region.label)
  
  subj_cts <- subj_cts %>%
    dplyr::left_join(sim_cts, by = 'da.region.label') %>%
    dplyr::mutate(pct_sim = sim_seqs / seqs_per_cluster)
}

# add the statistical info
subj_cts <- subj_cts %>%
  dplyr::left_join(DA.stat, by = 'da.region.label')

for (test_choice in tests){
  
  subj_cts[paste0(test_choice, '_significant')] <- subj_cts[test_choice] < 0.05
  
}

subj_cts <- subj_cts %>%
  dplyr::arrange(da.region.label, desc(pct_subj))

write.table(subj_cts, 
            file.path(OUTPUT_DIR, 'tables', 'cluster_cts.tsv'), 
            sep='\t', quote = F, row.names = F)

# make a summary of stats
stat_table <- data.frame('tool' = c('DAseq'),
                         'total_seqs' = nrow(X.cells),
                         'total_subj' = length(unique(X.cells$subject_id)),
                         'time (min)' = as.numeric(time_taken, units = "mins"),
                         'subjects' = paste(names(table(X.cells$subject_id)), collapse = ', '),
                         'depths' = paste(table(X.cells$subject_id), collapse = ', '),
                         check.names = F)

if (SIMULATED == T & 'simulated' %in% colnames(X.cells)){
  
  stat_table$pct_simulated <- c(mean(X.cells$simulated, na.rm = T) * 100)
  stat_table$Jaccard_0.005 = jaccard_005
  stat_table$Jaccard_0.05 = jaccard_05
  stat_table$Jaccard_0.1 = jaccard_1
  stat_table$Jaccard_max = Jaccard_max
  stat_table$Jaccard_max_p = Jaccard_max_p
  stat_table$AUC <- c(p_auroc)
  
  stat_table <- stat_table[c('tool', 'total_seqs', 'total_subj', 'pct_simulated',
                             'AUC', 'Jaccard_0.005', 'Jaccard_0.05',
                             'Jaccard_0.1', 'Jaccard_max', 'Jaccard_max_p',
                             'time (min)', 'subjects', 'depths')]
}

write.table(stat_table, 
            file.path(OUTPUT_DIR, 'tables', 'run_stats.tsv'), 
            sep = '\t', row.names = F, quote = F)

################################################################################
# DA region plot but only show significant

# add coordinates for DA regions for UMAP plotting
umap_embeddings_coords <- umap_embeddings_coords %>%
  dplyr::left_join(X.cells[c('id_col', 'da.region.label')], by = 'id_col')

umap_embeddings_coords$da.region.label <- as.character(umap_embeddings_coords$da.region.label)

prefix <- 'da.region.'

for (stat_choice in tests){
  
  # create new columns only highlighting significant clusters
  sig_clusters <- DA.stat %>%
    dplyr::filter(!!sym(stat_choice) < 0.05) %>%
    dplyr::pull(da.region.label)
  
  umap_embeddings_coords <- umap_embeddings_coords %>%
    dplyr::mutate('{prefix}{stat_choice}' := case_when(da.region.label %in% sig_clusters ~ da.region.label,
                                                       !da.region.label %in% sig_clusters ~ '0'))
  
}

# now generate a list of plots
plot_list <- lapply(paste0(prefix, tests), function(col_choice){
  
  # create a color palette - keep 0 gray always
  num_clusters <- max(as.numeric(umap_embeddings_coords$da.region.label), na.rm = T)
  all_clusters <- 0:num_clusters
  my_pal <- c("gray", hue_pal()(num_clusters))
  names(my_pal) <- as.character(all_clusters)
  
  p <- ggplot(umap_embeddings_coords, aes(x = UMAP1, 
                                          y = UMAP2, 
                                          color = as.factor(!!sym(col_choice)))) +
    rasterize(geom_point(alpha = 0.4, size = 0.5)) +           # semi-transparent points, moderate size
    theme_minimal(base_size = 18) +                # clean minimal theme with larger text
    labs(
      x = "UMAP 1",
      y = "UMAP 2",
      color = 'Region',
      title = str_replace(col_choice, prefix, '')
    ) +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    scale_color_manual(values = my_pal) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    theme_cowplot()
  
})

p_combined <- (plot_list[[1]] | plot_list[[2]]) / (plot_list[[3]] | plot_list[[4]]) +
  plot_layout(axis_titles = 'collect')

ggsave(file.path(OUTPUT_DIR, 'figures', 'da_regions_significant.png'), 
       device = 'png', width = 24, height = 24, units = 'in')

message(paste0('Finishing run: ', Sys.time()))
sessionInfo()
