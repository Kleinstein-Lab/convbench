#!/usr/bin/env Rscript
message(paste0('Starting run: ', Sys.time()))

suppressPackageStartupMessages({
  library(argparse)
  library(miloR)
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(dplyr)
  library(patchwork)
  library(stringr)
  library(Matrix)
  library(matrixStats)
  library(BiocNeighbors)
  library(cowplot)
  library(alakazam)
  library(RColorBrewer)
})

set.seed(37)


########################
### PREP ENVIRONMENT ###
########################

# prepare to take input parameters
parser <- ArgumentParser(description = "Data location and Milo algorithm hyperparameters.")

parser$add_argument('-d', '--data_loc', type = 'character', default = 'data',
                    help = 'File path for the embedding or RNA-Seq data location.')

parser$add_argument('-md', '--metadata_loc', type = 'character', default = 'metadata',
                    help = 'File path for the metadata location. Metadata and data files should have 1:1 matching sequence identifiers.')

parser$add_argument('-g', '--group_id', type = 'character', default = 'all',
                    help = 'Group ID (i.e. the ASC) used to prefix neighborhood IDs and output files.')

parser$add_argument('-o', '--output_dir', type = 'character', default = 'Milo_output',
                    help = 'Specify an output directory location.')

parser$add_argument('-da', '--da_variable', type = 'character', default = 'status',
                    help = 'Stratification variable used for differential abundance. Only used for visualization here; testing is done in evaluate.R.')

parser$add_argument('-k', '--k_val', type = 'integer', default = 50,
                    help = 'Number of neighbors to use in KNN algorithm.')

parser$add_argument('-pr', '--prop', type = 'double', default = 0.1,
                    help = 'Proportion of vertices to randomly sample.')

parser$add_argument('-a', '--auc_variable', type = 'character', default = FALSE,
                    help = 'Specify which column contains the ground truth for evaluation (i.e. "simulated" or "binder"). Column type should be logical. If no AUC variable, set to FALSE.')

parser$add_argument('-v', '--vdj_info', type = 'logical', default = FALSE,
                    help = 'Is v call and j call information included in the metadata? Can apply to expression or embedding data.')

# TODO: can change to be more granular/option to plot at gene, family etc. level
# right now defaults to v_call and j_call columns and removes allele info
parser$add_argument('-sc', '--single_cell', type = 'logical', default = FALSE,
                    help = 'Input true if V(D)J info is present and contains paired heavy and light chain info.')

parser$add_argument('-r', '--remove_dups', type = 'logical', default = FALSE,
                    help = 'Will remove duplicate embeddings within an individual if TRUE.')

# parser$add_argument('-g', '--use_glmm', type = 'logical', default = FALSE,
#                     help = 'Specify whether to account for subject in design formula.')

################################################################################

# Parse the arguments
args <- parser$parse_args()

# specify which dataset we are analyzing
DATA_LOC <- args$data_loc
MD_LOC <- args$metadata_loc
GROUP_ID <- args$group_id

OUTPUT_DIR <- args$output_dir

message(paste0('Data will be saved to ', OUTPUT_DIR, '.'))

# get the variable for DA calculations
DA_VAR <- args$da_variable

K_VAL <- args$k_val

PROP <- args$prop
message(paste0('Using index proportion: ', PROP))

message(paste0('K nearest neighbor value: ', K_VAL))

VDJ <- args$vdj_info
SINGLE_CELL <- args$single_cell
AUC_VAR <- args$auc_variable
REMOVE_DUPS <- args$remove_dups

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

if (AUC_VAR != FALSE){
  message(paste0('AUC variable ', AUC_VAR, ' will be saved for evaluation.'))
} else{
  message('No AUC variable provided.')
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
# FOR NOW, USE SUBJECT IDS AS SAMPLE IDS
# may need to change later on if these are different from each other!
if (!'sample_id' %in% colnames(md)){
  md$sample_id <- md$subject_id
} else{
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

################################################################################

# now make sure that metadata and embeddings are in same order
data <- data.frame(data, check.names = F)

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

  new_seq_num <- nrow(data)
  
  seqs_removed <- old_seq_num - new_seq_num
  message(paste0('Duplicates removed. ', seqs_removed, ' sequences removed. New total: ', new_seq_num))
  
} else{
  row.names(data) <- data$id_col
  data <- data %>% dplyr::select(-id_col)
}

# not all of the seqs in the data will necessarily result in successful embeddings
# so we can filter the metadata for only the relevant cell info
md <- md %>%
  dplyr::filter(id_col %in% row.names(data))

# the input for Milo is a SingleCellExperiment object, so we will create one
# make sure the rows are the same in the data and metadata
md <- data.frame(md, check.names = F)

# add gene, allele info
if (!'v_gene' %in% colnames(md)){
  md$v_gene <- alakazam::getGene(md$v_call, strip_d = F, omit_nl = F)
}

if (!'v_allele' %in% colnames(md)){
  md$v_allele <- alakazam::getAllele(md$v_call, strip_d = F, omit_nl = F)
}

if (!'j_gene' %in% colnames(md)){
  md$j_gene <- alakazam::getGene(md$j_call, strip_d = F, omit_nl = F)
}

if (!'j_allele' %in% colnames(md)){
  md$j_allele <- alakazam::getAllele(md$j_call, strip_d = F, omit_nl = F)
}

reduced_md_cols <- c(DA_VAR, 'subject_id', 'sample_id', 'id_col')

if (AUC_VAR != FALSE){
  reduced_md_cols <- c(reduced_md_cols, AUC_VAR)
}

md_reduced <- md %>%
  dplyr::select(all_of(reduced_md_cols)) %>%
  distinct()

row.names(md_reduced) <- md_reduced$id_col

# get a single list of cells found in the embeddings and the 
# filtered metadata so we can standardize the order
cells <- intersect(row.names(data), row.names(md_reduced))

md_reduced <- md_reduced[cells,]
data <- data[cells,]

# sanity check
message('Metadata and data rows aligned?')
message(all(row.names(md_reduced) == row.names(data)))

# need embeddings to be in matrix form for SCE
data_input <- as.matrix(data)

# NOTE: change the k value if it is too small
# if (length(cells) <= 100 & K_VAL > 5){
#   warning('Fewer than 100 cells - forcing K value to 5')
#   K_VAL <- 5
# } else if(length(cells) <= 500 & K_VAL > 10){
#   warning('Fewer than 500 cells - forcing K value to 10')
#   K_VAL <- 10
# }

# create the SCE
sce <- SingleCellExperiment(list(counts = t(data_input)))
colData(sce) <- DataFrame(md_reduced)
colnames(sce) <- colData(sce)$id_col

# Milo wants to used a reduction for the graph construction, so I will just feed
# the embedding information in
message('Single cell experiment object properly formatted?')
message(all(row.names(data_input) == colnames(sce))) # sanity check

# make umap for viz - runUMAP looks for log counts but we will just use 
# the embedding value PCs
# run PCA if enough data
if (nrow(data) >= 200){
  message('Using first 200 PCs to generate UMAP...')
  pca <- prcomp(data_input, center = T, scale. = T)
  reducedDim(sce, 'PCA') <- pca$x[, 1:200] # use 200 PCs
  sce <- runUMAP(sce, dimred = 'PCA', n_neighbors = K_VAL)
} else{
  message('Generating UMAP from all data...')
  sce <- runUMAP(sce, exprs_values = "counts", n_neighbors = K_VAL)
}

reducedDim(sce, 'embedding') <- data_input

################################################################################
# measure how long the Milo process itself takes
start_time <- Sys.time()

# make it into a Milo object
milo <- Milo(sce)
reducedDim(milo, "UMAP") <- reducedDim(sce, "UMAP")

message(Sys.time())
message('Building KNN graph')

# next, build KNN graph
milo <- buildGraph(milo,
                   k = K_VAL,
                   d = length(colnames(data)),
                   reduced.dim = 'embedding')

message(Sys.time())
message('Defining representative neighborhoods')

# now make neighborhoods of indices
milo <- makeNhoods(milo, 
                   prop = PROP, 
                   k = K_VAL, 
                   d = length(colnames(data)), 
                   refined = TRUE,
                   reduced_dims = 'embedding')

# distribution should peak between at a point that makes sense for the 
# neighborhood sizes we are anticipating
nhood_dist <- plotNhoodSizeHist(milo)

ggsave(file.path(OUTPUT_DIR, 'figures', 'neighborhood_size_dist.png'),
       plot = nhood_dist, device = 'png', width = 10, height = 8, units = 'in')

message('Counting cells')

milo <- countCells(milo, 
                   meta.data = data.frame(colData(milo)), 
                   samples="sample_id")

Sys.time()
message('Calculating distances between nearest neighbors')
milo <- calcNhoodDistance(milo,
                          d = length(colnames(data)),
                          reduced.dim = 'embedding')

end_time <- Sys.time()
time_taken <- end_time - start_time

##################################
### SAVE INPUTS FOR EVALUATION ###
##################################
# the GLM, Fisher and Wilcoxon tests are done in evaluate.R

nhood_ids <- as.character(unlist(milo@nhoodIndex))
cluster_ids <- paste0(GROUP_ID, '_', nhood_ids)

# neighborhood x subject count matrix
nhood_counts <- as.matrix(nhoodCounts(milo))
if (!all(as.numeric(row.names(nhood_counts)) == seq_along(nhood_ids))){
  stop('Neighborhood count matrix rows are not in neighborhood order.')
}

write.table(data.frame(cluster_id = cluster_ids, nhood_counts, check.names = F),
            file.path(OUTPUT_DIR, 'tables', paste0(GROUP_ID, '_cluster_counts.tsv')),
            sep = '\t', row.names = F, quote = F)

# sequence-to-neighborhood memberships, in neighborhood order. Sequences can be in
# multiple neighborhoods, and sequences in no neighborhood are kept with an NA cluster ID
memberships <- Matrix::summary(as(nhoods(milo), 'CsparseMatrix'))
memberships <- memberships[memberships$x != 0, ]
memberships <- memberships[order(memberships$j, memberships$i), ]

memberships <- data.frame(id_col = row.names(nhoods(milo))[memberships$i],
                          cluster_id = cluster_ids[memberships$j])

memberships <- md_reduced %>%
  dplyr::left_join(memberships, by = 'id_col', relationship = 'one-to-many') %>%
  dplyr::mutate(sequence_id = id_col, group_id = GROUP_ID)

membership_cols <- c('sequence_id', 'subject_id', 'group_id', 'cluster_id')
if (AUC_VAR != FALSE){
  membership_cols <- c(membership_cols, AUC_VAR)
}

write.table(memberships[membership_cols],
            file.path(OUTPUT_DIR, 'tables', paste0(GROUP_ID, '_memberships.tsv')),
            sep = '\t', row.names = F, quote = F)

write.table(data.frame(group_id = GROUP_ID, time_min = as.numeric(time_taken, units = 'mins')),
            file.path(OUTPUT_DIR, 'tables', paste0(GROUP_ID, '_tool_info.tsv')),
            sep = '\t', row.names = F, quote = F)

# save the Milo object for the GLM and spatial FDR. To save space, drop what evaluation
# does not need: the counts assay (a copy of the embedding, which is kept as a reduced
# dimension), the PCA, and the neighborhood distances (not used by neighbour-distance weighting)
milo_eval <- milo
assay(milo_eval, 'counts') <- Matrix(0, nrow = nrow(milo_eval), ncol = ncol(milo_eval), sparse = TRUE,
                                     dimnames = dimnames(milo_eval))
reducedDim(milo_eval, 'PCA') <- NULL
milo_eval@nhoodDistances <- list()
saveRDS(milo_eval, file.path(OUTPUT_DIR, 'tables', paste0(GROUP_ID, '_milo.rds')))
rm(milo_eval)

message(paste0(length(cluster_ids), ' neighborhoods saved for evaluation.'))

###################
### DIAGNOSTICS ###
###################

dispersion_df <- data.frame(nhood_id = unlist(milo@nhoodIndex),
                            Mean_Counts = rowMeans(milo@nhoodCounts),
                            Var_Counts = rowVars(milo@nhoodCounts),
                            Total_Counts = rowSums(milo@nhoodCounts),
                            Sharing_Num = rowSums(milo@nhoodCounts > 0))

dispersion_df %>%
  ggplot(aes(x = Mean_Counts, y = Var_Counts)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  theme_bw()

ggsave(file.path(OUTPUT_DIR, 'figures', 'dispersion_fig_plain.png'),
       device = 'png', width = 6, height = 6, units = 'in')


dispersion_df %>%
  ggplot(aes(x = Mean_Counts, y = Var_Counts)) +
  geom_point(aes(color = factor(Sharing_Num),
                 size = Total_Counts), alpha = 0.65) +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  theme_bw()

ggsave(file.path(OUTPUT_DIR, 'figures', 'dispersion_fig_fancy1.png'),
       device = 'png', width = 8, height = 6, units = 'in')

dispersion_df %>%
  ggplot(aes(x = Mean_Counts, y = Var_Counts)) +
  geom_point(aes(size = Total_Counts), alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  theme_bw()

ggsave(file.path(OUTPUT_DIR, 'figures', 'dispersion_fig_fancy2.png'),
       device = 'png', width = 8, height = 6, units = 'in')

################################################################################
###############
### SUMMARY ###
###############

# percent subject in various neighborhoods
nhood_sizes <- data.frame(colSums(milo@nhoods))
colnames(nhood_sizes) <- c('cells_per_nhood')
nhood_sizes$nhood_id <- row.names(nhood_sizes)

# tabulate cells from each person in each neighborhood
subj_nhood_cts <- lapply(unique(md$subject_id), function(subj){
  
  subj_cells <- md_reduced %>%
    dplyr::filter(subject_id == subj) %>%
    row.names()
  
  nhood_cts <- lapply(nhood_sizes$nhood_id, function(nhood){
    
    cells_per_subj <- sum(milo@nhoods[subj_cells,nhood])
    
    return(data.frame('nhood_id' = nhood,
                      'subject_id' = subj,
                      'cells_per_subj' = cells_per_subj))
    
  })
  
  return(do.call(rbind, nhood_cts))
  
})

subj_nhood_cts <- do.call(rbind, subj_nhood_cts)

subj_nhood_cts <- subj_nhood_cts %>%
  dplyr::left_join(nhood_sizes,
                   by = 'nhood_id')

subj_nhood_cts <- subj_nhood_cts %>%
  dplyr::mutate(pct_subj = cells_per_subj / cells_per_nhood)

# do the same for samples if not redundant
if (!all(sort(unique(md$subject_id)) == sort(unique(md$sample_id)))){
  
  samp_nhood_cts <- lapply(unique(md$sample_id), function(samp){
    
    samp_cells <- md_reduced %>%
      dplyr::filter(sample_id == samp) %>%
      row.names()
    
    nhood_cts <- lapply(nhood_sizes$nhood_id, function(nhood){
      
      cells_per_samp <- sum(milo@nhoods[samp_cells,nhood])
      
      return(data.frame('nhood_id' = nhood,
                        'sample_id' = samp,
                        'cells_per_samp' = cells_per_samp))
      
    })
    
    return(do.call(rbind, nhood_cts))
    
  })
  
  samp_nhood_cts <- do.call(rbind, samp_nhood_cts)
  
  samp_nhood_cts <- samp_nhood_cts %>%
    dplyr::left_join(nhood_sizes,
                     by = 'nhood_id')
  
  samp_nhood_cts <- samp_nhood_cts %>%
    dplyr::mutate(pct_samp = cells_per_samp / cells_per_nhood)
  
  samp_nhood_cts <- samp_nhood_cts %>%
    dplyr::left_join(distinct(md_reduced[c('subject_id', 'sample_id')]), by = 'sample_id')
  
  # TODO: figure out what to do with sample level information
  # join on subj id and nhood id?
  subj_nhood_cts <- samp_nhood_cts %>% 
    dplyr::left_join(subj_nhood_cts, by = c('nhood_id', 'subject_id', 'cells_per_nhood'))
  
}

# also get percent sim
if (AUC_VAR != FALSE){
  
  hit_cells <- md %>%
    dplyr::filter(!!sym(AUC_VAR) == TRUE) %>%
    pull(id_col)
  
  hit_pct <- lapply(nhood_sizes$nhood_id, function(nhood){
    
    total_hits <- sum(milo@nhoods[hit_cells, nhood])
    
    return(data.frame('nhood_id' = nhood,
                      'hit_seqs' = total_hits))
    
  })
  
  hit_pct <- do.call(rbind, hit_pct)
  
  subj_nhood_cts <- subj_nhood_cts %>%
    dplyr::left_join(hit_pct,
                     by = 'nhood_id')
  
  # get percent of nhood that is hits seqs
  subj_nhood_cts <- subj_nhood_cts %>%
    dplyr::mutate(pct_hits = hit_seqs / cells_per_nhood)
  
}

# add additional info: pairwise dist stats
dist_mat_stats_list <- lapply(names(milo@nhoodDistances), function(nh){
  mat <- milo@nhoodDistances[[nh]]
  
  # get upper tri vals
  upper_tri <- mat[lower.tri(mat)]
  
  data.frame(
    nhood_id = nh,
    n_pairs = length(upper_tri),
    mean_dist = mean(upper_tri),
    median_dist = median(upper_tri),
    sd_dist = sd(upper_tri),
    min_dist = min(upper_tri),
    max_dist = max(upper_tri)
  )
  
})

dist_mat_stats <- do.call(rbind, dist_mat_stats_list)

# add to subj_nhood_cts, do re-runs of bg
subj_nhood_cts <- subj_nhood_cts %>%
  dplyr::left_join(dist_mat_stats, by = 'nhood_id', relationship = 'many-to-one')

write.table(subj_nhood_cts, 
            file.path(OUTPUT_DIR, 'tables', 'nhood_stats.tsv'), 
            sep = '\t', row.names = F, quote = F)

################
### DATA VIZ ###
################
# FOR UMAP VIZ
umap_coords <- data.frame(reducedDim(milo, 'UMAP'))
umap_coords$id_col <- row.names(umap_coords)

umap_coords <- umap_coords %>%
  dplyr::left_join(md_reduced, by = 'id_col')

make_UMAP_viz <- function(var, var_name, custom_pal = NULL){
  
  default_h <- 12
  default_w <- 14
  
  umap_coords[[var]] <- as.factor(umap_coords[[var]])
  
  # let's get some visualizations first
  p <- ggplot(umap_coords, aes(x = UMAP1, y = UMAP2, color = !!sym(var))) +
    geom_point(alpha = 0.6, size = 0.6) +
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
  
  num_vars <- length(unique(umap_coords[[var]]))
  
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

if (nrow(umap_coords) >= 400000){
  message('Skipping UMAP visualization because dataset too large.')
} else{
  message('Making UMAP visualizations...')
}

if (VDJ && nrow(umap_coords) < 400000){
  
  if(!('v_gene' %in% colnames(md) & 'j_gene' %in% colnames(md))) warning('v_gene and j_gene columns not provided. UMAP plots for V and J gene will not be generated.')
  
  else {
    
    if (!SINGLE_CELL){
      
      umap_coords <- umap_coords %>%
        dplyr::left_join(md[c('v_gene', 'j_gene', 'id_col')], by = 'id_col')
      # 
      # # get rid of alleles for v and j call
      # umap_coords$v_gene <- str_replace(umap_coords$v_call, '\\*.*', '')
      # umap_coords$j_gene <- str_replace(umap_coords$j_call, '\\*.*', '')
      
      message('Visualize V and J genes.')
      make_UMAP_viz('v_gene', 'V gene')
      make_UMAP_viz('j_gene', 'J gene')
      
    } else if (SINGLE_CELL){
      
      ##### APPLIES TO SINGLE CELL ONLY #####
      # get heavy and light chain V/J assignments
      heavy_info <- umap_coords %>%
        dplyr::left_join(md[c('v_gene', 'j_gene', 'locus', 'id_col')], by = 'id_col') %>%
        dplyr::filter(id_col %in% row.names(data)) %>%
        dplyr::filter(locus == 'IGH') %>%
        dplyr::select(id_col, v_gene, j_gene) %>%
        distinct() %>%
        data.frame(check.names = F)
      
      row.names(heavy_info) <- heavy_info$id_col
      
      con_ct_col <- NA
      
      if ('consensus_count' %in% colnames(md)){
        con_ct_col <- 'consensus_count'
      } else if ('conscount' %in% colnames(md)){
        con_ct_col <- 'conscount'
      } else{
        message('Consensus count column not found. light chain plot will not be generated.')
      }
      
      if (!is.na(con_ct_col)){
        
        light_info <- umap_coords %>%
          dplyr::left_join(md[c('v_gene', 'j_gene', 'id_col', 'locus', con_ct_col)], by = 'id_col') %>%
          dplyr::filter(id_col %in% row.names(data)) %>%
          dplyr::filter(locus == 'IGK' | locus == 'IGL') %>%
          dplyr::group_by(id_col) %>%
          dplyr::arrange(desc(!!sym(con_ct_col))) %>%
          dplyr::slice_head(n = 1) %>%
          dplyr::ungroup() %>%
          dplyr::select(id_col, v_gene, j_gene) %>%
          distinct() %>%
          data.frame(check.names = F)
        
        if (nrow(light_info) > 0){
          
          row.names(light_info) <- light_info$id_col
          
          umap_coords <- umap_coords %>%
            dplyr::left_join(light_info, by = 'id_col') %>%
            dplyr::rename(v_gene_light = v_gene,
                          j_gene_light = j_gene)
          
          make_UMAP_viz('v_gene_light', 'V Gene - \nLight Chain')
          make_UMAP_viz('j_gene_light', 'J Gene - \nLight Chain')
          
        }
        
      }
      
      umap_coords <- umap_coords %>%
        dplyr::left_join(heavy_info, by = 'id_col') %>%
        dplyr::rename(v_gene_heavy = v_gene,
                      j_gene_heavy = j_gene)
      
      # umap_coords$v_gene_heavy <- getGene(umap_coords$v_call_heavy)
      # umap_coords$j_gene_heavy <- getGene(umap_coords$j_call_heavy)
      # umap_coords$v_gene_light <- getGene(umap_coords$v_call_light)
      # umap_coords$j_gene_light <- getGene(umap_coords$j_call_light)
      
      ##### APPLIES TO SINGLE CELL ONLY #####
      make_UMAP_viz('v_gene_heavy', 'V Gene - \nHeavy Chain')
      make_UMAP_viz('j_gene_heavy', 'J Gene - \nHeavy Chain')

      
    }
  }
}

# include info if simulated
if (AUC_VAR != FALSE && nrow(umap_coords) < 400000){
  
  message('Visualize hit sequences.')
  
  make_UMAP_viz(AUC_VAR, 'Hits', custom_pal = c('TRUE' = "red", 'FALSE' = "gray"))
  
}

if (nrow(umap_coords) < 400000){
  message('Visualize subject and sample information.')

  make_UMAP_viz(DA_VAR, DA_VAR)
  make_UMAP_viz('sample_id', 'Sample ID')
  make_UMAP_viz('subject_id', 'Subject ID')
}

message(paste0('Ending run: ', Sys.time()))

sessionInfo()
