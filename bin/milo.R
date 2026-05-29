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
  library(pracma)
  library(cowplot)
  library(alakazam)
  library(RColorBrewer)
})

set.seed(37)

# prepare to take input parameters
parser <- ArgumentParser(description = "Data location and Milo algorithm hyperparameters.")

parser$add_argument('-d', '--data_loc', type = 'character', default = 'data',
                    help = 'File path for the embedding or RNA-Seq data location.')

parser$add_argument('-md', '--metadata_loc', type = 'character', default = 'metadata',
                    help = 'File path for the metadata location. Metadata and data files should have 1:1 matching sequence identifiers.')

parser$add_argument('-o', '--output_dir', type = 'character', default = 'DAseq_output',
                    help = 'Specify an output directory location.')

parser$add_argument('-da', '--da_variable', type = 'character', default = 'status',
                    help = 'Stratification variable that should be used to determine for differential abundance. There should be two levels in this factor/categorical variable.')

parser$add_argument('-k', '--k_val', type = 'integer', default = 50,
                    help = 'Number of neighbors to use in KNN algorithm.')

parser$add_argument('-pr', '--prop', type = 'double', default = 0.1,
                    help = 'Proportion of vertices to randomly sample.')

parser$add_argument('-a', '--auc_variable', type = 'character', default = FALSE,
                    help = 'Specify which column should be used for generating AUC curve (i.e. "simulated" or "binder"). Column type should be logical. If no AUC variable, set to FALSE.')

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

parser$add_argument('-w', '--overwrite', type = 'logical', default = FALSE,
                    help = 'Specify whether to re-create and write new Milo object.')

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

K_VAL <- args$k_val

PROP <- args$prop
message(paste0('Using index proportion: ', PROP))

message(paste0('K nearest neighbor value: ', K_VAL))

VDJ <- args$vdj_info
SINGLE_CELL <- args$single_cell
AUC_VAR <- args$auc_variable
OVERWRITE <- args$overwrite
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

if (!'v_gene' %in% colnames(md)){
  # assume if v_gene not included, J probably isn't either
  md$v_gene <- alakazam::getGene(md$v_call, strip_d = F, omit_nl = F)
  md$v_allele <- alakazam::getAllele(md$v_call, strip_d = F, omit_nl = F)
  md$j_gene <- alakazam::getGene(md$j_call, strip_d = F, omit_nl = F)
  
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

if (!file.exists(file.path(OUTPUT_DIR, 'tables', 'milo.RDS')) | OVERWRITE == T){
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
    pca <- prcomp(data_input, center = T, scale. = T)
    reducedDim(sce, 'PCA') <- pca$x[, 1:200] # use 200 PCs
    sce <- runUMAP(sce, dimred = 'PCA', n_neighbors = K_VAL)
  } else{
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
  
  design_df <- data.frame(colData(milo))[,c('sample_id', 'subject_id', DA_VAR)]
  design_df <- distinct(design_df)
  rownames(design_df) <- design_df$sample_id
  ## Reorder rownames to match columns of nhoodCounts(milo)
  design_df <- design_df[colnames(nhoodCounts(milo)), , drop=FALSE]
  
  design_df$sample_id <- as.factor(design_df$sample_id)
  design_df$subject_id <- as.factor(design_df$subject_id)
  design_df[,DA_VAR] <- as.factor(design_df[,DA_VAR])
  
  print('Design:')
  print(design_df)
  
  print(table(data.frame(colData(milo))$subject_id))
  
  Sys.time()
  message('Calculating distances between nearest neighbors')
  milo <- calcNhoodDistance(milo,
                            d = length(colnames(data)),
                            reduced.dim = 'embedding')
  
  # NOTE: not using GLMM currently, but could be implemented if needed
  # if (USE_GLMM){
  #   
  #   formula_string <- paste0(' ~ ', DA_VAR, ' + (1|subject_id)')
  #   
  #   design_formula <- as.formula(formula_string)
  #   
  #   message(paste0('Using formula: ', formula_string))
  #   
  #   da_results <- testNhoods(milo,
  #                            design = design_formula, 
  #                            design.df = design_df,
  #                            reduced.dim = 'embedding',
  #                            fdr.weighting = 'neighbour-distance',
  #                            glmm.solver = 'Fisher',
  #                            norm.method = 'TMM',
  #                            REML = TRUE)
  #   
  # } else{
  
  formula_string <- paste0('~ ', DA_VAR)
  
  design_formula <- as.formula(formula_string)
  
  message(paste0('Using formula: ', formula_string))
  
  da_results <- testNhoods(milo, 
                           norm.method = 'logMS',
                           design = design_formula, 
                           design.df = design_df,
                           reduced.dim = 'embedding',
                           fdr.weighting = 'neighbour-distance')
  
  # }
  
  print('Top DA results:')
  da_results %>%
    arrange(SpatialFDR) %>%
    head() %>%
    print()
  
  ########### ADDED DIAGNOSTIC VISUALS ########### 
  da_results %>%
    ggplot(aes(x = PValue)) + 
    geom_histogram(color = 'white', binwidth = 0.01) + 
    theme_bw() +
    labs(title = 'Milo P-Value Distribution',
         subtitle = 'NON-Permuted labels') +
    coord_cartesian(xlim = c(0, 1))
  
  ggsave(file.path(OUTPUT_DIR, 'figures', 'pvalue_hist.png'),
         device = 'png', width = 8, height = 6, units = 'in')
  
  da_results %>%
    ggplot(aes(x = SpatialFDR)) + 
    geom_histogram(color = 'white', binwidth = 0.01) + 
    theme_bw() +
    labs(title = 'Milo Spatial FDR Distribution',
         subtitle = 'NON-Permuted labels') +
    coord_cartesian(xlim = c(0, 1))
  
  ggsave(file.path(OUTPUT_DIR, 'figures', 'spatialFDR_hist.png'),
         device = 'png', width = 8, height = 6, units = 'in')
  
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
  
  #####################################################################
  # viz
  milo <- buildNhoodGraph(milo)
  
  plotUMAP(milo) + 
    plotNhoodGraphDA(milo, da_results, alpha=0.05) +
    plot_layout(guides="collect")
  
  ggsave(file.path(OUTPUT_DIR, 'figures', 'final_UMAP.png'), 
         device = 'png',  width = 12, height = 6, units = 'in')
  
  # add neighborhood ids to results
  da_results$nhood_id <- as.character(unlist(milo@nhoodIndex))
  
  ################################################################################
  end_time <- Sys.time()
  time_taken <- end_time - start_time

  # save the milo obj for later use
  # saveRDS(milo, file.path(OUTPUT_DIR, 'tables', 'milo.RDS'))
  saveRDS(milo@nhoods, file.path(OUTPUT_DIR, 'tables', 'milo_nhoods.RDS'))
  
} else{
  
  cat(paste0('Loading Milo object from: ', file.path(OUTPUT_DIR, 'tables', 'da_cells.rds')))
  milo <- readRDS(file.path(OUTPUT_DIR, 'tables', 'milo.RDS'))
  da_results <- read.csv(file.path(OUTPUT_DIR, 'tables', 'da_results.tsv'), sep = '\t')
  da_results$nhood_id <- as.character(da_results$nhood_id)
  
  # preserve the original time taken...if somehow no run stat table, just set to NA
  if (file.exists(file.path(OUTPUT_DIR, 'tables', 'run_stats.tsv'))){
    run_stat_existing <- read.csv(file.path(OUTPUT_DIR, 'tables', 'run_stats.tsv'), sep = '\t', check.names = F)
    time_taken <- run_stat_existing[['time (min)']]
  } else{
    time_taken <- NA
  } 
  
}

nhoods_match <- all(colnames(milo@nhoods) == da_results$nhood_id)
if (!nhoods_match){
  message('WARNING: neighborhood matrix does not match results table. AUC values will not be accurate.')
}


# match each cell with the lowest p-value of all the neighborhoods it occupies
if (!file.exists(file.path(OUTPUT_DIR, 'tables', 'seq_results.tsv')) | OVERWRITE == T){
  min_p_nhoods <- lapply(row.names(milo@nhoods), function(current_seq){
    
    test <- milo@nhoods[current_seq,]*da_results$SpatialFDR
    test <- test[test != 0]
    
    if (length(test) == 0){
      result <- NA
      min_p_clust <- NA
      min_p_logFC <- NA
    } else{
      result <- min(test)
      min_p_clust <- names(which.min(test))
      min_p_logFC <- da_results %>% 
        dplyr::filter(nhood_id == min_p_clust) %>% 
        dplyr::pull(logFC)
    }
    
    return(data.frame('id_col' = current_seq,
                      'min_nhood_id' = min_p_clust,
                      'min_nhood_FDR' = result,
                      'min_nhood_logFC' = min_p_logFC
    ))
    
  })
  
  min_p_nhoods_df <- do.call(rbind, min_p_nhoods)
  
  write.table(min_p_nhoods_df, 
              file.path(OUTPUT_DIR, 'tables', 'seq_results.tsv'), 
              sep = '\t', row.names = F, quote = F)
} else{
  min_p_nhoods_df <- read.csv(file.path(OUTPUT_DIR, 'tables', 'seq_results.tsv'), sep = '\t')
}

# get a continuous DA measure - copy of benchmark - sum of logFC of all neighborhoods
# da.cell.mat <- milo@nhoods %*% da_results$logFC
# da.cell <- da.cell.mat[,1]

# return with cell_id and UMAP stats

# but it looks like we need to do the alphas version?? Need to double check code:
# https://github.com/CompCy-lab/benchmarkDA/blob/09c4b20a6b36a3633d327b551374915edc27108d/scripts/benchmark_utils.R#L472
# and language here: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03143-0#Sec11

# I am not really sure what the purpose of the continuous version is tbh
# but it seems like for the other methods not daseq (and meld I think) we are in fact using an FDR threshold method

##################
### EVALUATION ###
##################

#######
# AUC #
#######

# updated to deal with cells that have no nhood

if (AUC_VAR != FALSE){
  invalid_cells <- sum(is.na(min_p_nhoods_df$min_nhood_FDR))
  total_cells <- nrow(min_p_nhoods_df)
  valid_cells <- total_cells - invalid_cells
  
  # change the sequences with no nhood to a min p of 1
  min_p_nhoods_df[is.na(min_p_nhoods_df$min_nhood_FDR), 'min_nhood_FDR'] <- 1
  
  # add AUC var info
  min_p_nhoods_df <- min_p_nhoods_df %>%
    dplyr::inner_join(md_reduced, by = 'id_col')
  
  auc_thresholds <- sort(unique(min_p_nhoods_df$min_nhood_FDR))
  # auc_thresholds <- quantile(min_p_nhoods_df$min_nhood_FDR, seq(0, 1, 0.01), names=F)
  # auc_thresholds <- quantile(min_p_nhoods_df$min_nhood_FDR, seq(0, 1, 0.01), names=F)
  
  # auc_thresholds[1] <- auc_thresholds[1] - 1e-8
  
  # add to the largest to make sure the entire curve is captured
  tot_thresh <- length(auc_thresholds)
  auc_thresholds[tot_thresh] <- auc_thresholds[tot_thresh] + 1e-3
  
  auc_data <- lapply(auc_thresholds, function(thresh){
    
    # get cells with min nhood p below threshold
    da.cell.list <- min_p_nhoods_df$min_nhood_FDR < thresh
    
    true_pos <- sum(da.cell.list == T & min_p_nhoods_df[[AUC_VAR]] == T)
    false_neg <- sum(da.cell.list == F & min_p_nhoods_df[[AUC_VAR]] == T)
    true_neg <- sum(da.cell.list == F & min_p_nhoods_df[[AUC_VAR]] == F)
    false_pos <- sum(da.cell.list == T & min_p_nhoods_df[[AUC_VAR]] == F)
    
    TPR <- true_pos / (true_pos + false_neg)
    FPR <- 1 - (true_neg / (true_neg + false_pos))
    
    return(data.frame('TPR' = TPR,
                      'FPR' = FPR))
    
  })
  
  auc_df <- do.call(rbind, auc_data)
  auc_df$spatialFDR_threshold <- auc_thresholds
  
  write.table(auc_df, 
              file.path(OUTPUT_DIR, 'tables', 'auc_curve_vals.tsv'), 
              sep = '\t', row.names = F, quote = F)
  
  # get auroc
  auroc <- pracma::trapz(auc_df$FPR, auc_df$TPR)
  
  auc_df %>%
    ggplot(aes(x = FPR, y = TPR)) +
    geom_point() +
    geom_line() +
    labs(title = paste0('Alpha Threshold ', round(min(auc_thresholds)), ' to ', round(max(auc_thresholds), 3)),
         subtitle = paste0('AUC: ', round(auroc, 3), '; ', 
                           prettyNum(sum(valid_cells), big.mark = ",", scientific = FALSE), '/', 
                           prettyNum(total_cells, big.mark = ",", scientific = FALSE), ' cells in DA neighborhoods')) + 
    theme_minimal()
  
  ggsave(file.path(OUTPUT_DIR, 'figures', 'AUC_curve.png'),
         device = 'png',
         width = 7,
         height = 6)
  
  ###########
  # JACCARD #
  ###########
  jaccard_df <- min_p_nhoods_df %>%
    dplyr::mutate(p_under_0.005 = min_nhood_FDR <= 0.005,
                  p_under_0.05 = min_nhood_FDR <= 0.05,
                  p_under_0.1 = min_nhood_FDR <= 0.1) 
  
  # calc jaccard index
  jaccard_005 <- sum(jaccard_df[[AUC_VAR]] & jaccard_df$p_under_0.005, na.rm = T) / sum(jaccard_df[[AUC_VAR]] | jaccard_df$p_under_0.005, na.rm = T)
  jaccard_05 <- sum(jaccard_df[[AUC_VAR]] & jaccard_df$p_under_0.05, na.rm = T) / sum(jaccard_df[[AUC_VAR]] | jaccard_df$p_under_0.05, na.rm = T)
  jaccard_1 <- sum(jaccard_df[[AUC_VAR]] & jaccard_df$p_under_0.1, na.rm = T) / sum(jaccard_df[[AUC_VAR]] | jaccard_df$p_under_0.1, na.rm = T)
  
  jaccard_thresholds <- sort(unique(jaccard_df$min_nhood_FDR))
  jaccard_thresholds <- jaccard_thresholds[!is.na(jaccard_thresholds)]
  
  # get Jaccard across a range
  jaccards <- sapply(jaccard_thresholds, function(thresh){
    j <- sum(jaccard_df[[AUC_VAR]] & jaccard_df$min_nhood_FDR <= thresh, na.rm = T) / sum(jaccard_df[[AUC_VAR]] | jaccard_df$min_nhood_FDR <= thresh, na.rm = T)
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
  auroc <- NA
}

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
    dplyr::filter(`AUC_VAR` == TRUE) %>%
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

write.table(da_results, 
            file.path(OUTPUT_DIR, 'tables', 'da_results.tsv'), 
            sep = '\t', row.names = F, quote = F)

subj_nhood_cts <- subj_nhood_cts %>%
  dplyr::left_join(da_results[c('nhood_id', 'logFC', 'SpatialFDR')],
                   by = 'nhood_id')

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

# make a summary of stats
stat_table <- data.frame('tool' = c('Milo'),
                         'total_seqs' = c(ncol(milo)),
                         'total_subj' = ncol(milo@nhoodCounts),
                         'time (min)' = as.numeric(time_taken, units = "mins"),
                         'subjects' = paste(names(table(milo@colData$subject_id)), collapse = ', '),
                         'depths' = paste(table(milo@colData$subject_id), collapse = ', '),
                         check.names = F)

if (AUC_VAR != FALSE){
  
  stat_table$pct_hits <- c(mean(milo@colData[[AUC_VAR]], na.rm = T) * 100)
  stat_table$Jaccard_0.005 = jaccard_005
  stat_table$Jaccard_0.05 = jaccard_05
  stat_table$Jaccard_0.1 = jaccard_1
  stat_table$Jaccard_max = Jaccard_max
  stat_table$Jaccard_max_p = Jaccard_max_p
  stat_table$AUC <- c(auroc)
  
  stat_table <- stat_table[c('tool', 'total_seqs', 'total_subj', 'pct_hits',
                             'AUC', 'Jaccard_0.005', 'Jaccard_0.05',
                             'Jaccard_0.1', 'Jaccard_max', 'Jaccard_max_p',
                             'time (min)', 'subjects', 'depths')]
}

write.table(stat_table, 
            file.path(OUTPUT_DIR, 'tables', 'run_stats.tsv'), 
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

if (VDJ){
  
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
if (AUC_VAR != FALSE){
  
  message('Visualize hit sequences.')
  
  make_UMAP_viz(AUC_VAR, 'Hits', custom_pal = c('TRUE' = "red", 'FALSE' = "gray"))
  
}

message('Visualize subject and sample information.')

make_UMAP_viz(DA_VAR, DA_VAR)
make_UMAP_viz('sample_id', 'Sample ID')
make_UMAP_viz('subject_id', 'Subject ID')

# make viz for all cells in neighborhoods significant at alpha 0.05
sig_nhoods <- da_results %>%
  dplyr::filter(SpatialFDR < 0.05) %>%
  dplyr::pull(nhood_id)

sig_nhood_cells <- milo@nhoods[,as.character(sig_nhoods)]
sig_nhood_idx <- rowSums(sig_nhood_cells) > 0

sig_nhood_cell_ids <- row.names(sig_nhood_cells[sig_nhood_idx,])

umap_coords <- umap_coords %>%
  dplyr::mutate(da_cell = if_else(id_col %in% sig_nhood_cell_ids, TRUE, FALSE))

message('Visualize significant neighborhoods.')

make_UMAP_viz('da_cell', 'in significant neighborhood \n (spatial FDR < 0.05)', custom_pal = c('TRUE' = "red", 'FALSE' = "gray"))

message(paste0('Ending run: ', Sys.time()))

sessionInfo()
