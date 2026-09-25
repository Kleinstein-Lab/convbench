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

########################
### HELPER FUNCTIONS ###
########################

make_purity_plot <- function(purity_data, cluster_id_col, pct_hit_col, total_seq_col, auc_variable){

  ggplot(purity_data, aes(x = !!sym(cluster_id_col), y = !!sym(pct_hit_col))) +
    geom_col(fill = 'dodgerblue3') +
    geom_text(
      aes(label = !!sym(total_seq_col)),
      vjust = -0.5
    ) +
    scale_y_continuous(
      labels = function(x) paste0(x * 100, "%"),
      limits = c(0, 1.05)
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    ) +
    labs(x = 'Cluster ID',
         y = paste0('Percent ', auc_variable))
    
  ggsave(file.path(OUTPUT_DIR, 'figures', 'cluster_purity.png'), 
         device="png", width=10, height=5, units="in")

}

run_nhood_fisher <- function(nhood_counts, subj_info, da_variable, disease_group){
  # perform a one-sided Fisher Exact Test on each neighborhood based on SUBJECTS in the neighborhood,
  # testing whether subjects in disease_group are over-represented
  # nhood_counts = neighborhood x subject count matrix (i.e. milo@nhoodCounts). Should already include
  #                subjects absent from the current subset so totals reflect the WHOLE dataset
  # subj_info = one row per subject with subject IDs as row names and a da_variable column (i.e. lib_sizes)
  # disease_group = level of da_variable we are testing for enrichment

  subj_info <- subj_info[colnames(nhood_counts), , drop = FALSE]
  is_dis <- subj_info[[da_variable]] == disease_group

  tot_cond <- sum(is_dis)
  tot_not_cond <- sum(!is_dis)

  # count subjects with at least one sequence in each neighborhood, with and without condition
  present <- nhood_counts > 0
  in_nhood_cond <- unname(Matrix::rowSums(present[, is_dis, drop = FALSE]))
  in_nhood_not_cond <- unname(Matrix::rowSums(present[, !is_dis, drop = FALSE]))

  fisher_res_list <- lapply(seq_len(nrow(nhood_counts)), function(i){

    # build contingency table to test for a CONDITION neighborhood
    #
    #                   nhood
    #                No    Yes
    #               ___________
    #            No|     |     |
    # condition    |_____|_____|
    #           Yes|     |     |
    #              |_____|_____|

    contingency_table <- matrix(c(tot_not_cond - in_nhood_not_cond[i], tot_cond - in_nhood_cond[i],
                                  in_nhood_not_cond[i], in_nhood_cond[i]), 2, 2)

    fisher <- fisher.test(contingency_table, alternative = 'greater')

    return(data.frame(PValue = fisher$p.value,
                      odds_ratio = unname(fisher$estimate)))
  })

  fisher_result <- do.call(rbind, fisher_res_list)
  fisher_result$in_nhood_in_condition <- in_nhood_cond
  fisher_result$in_nhood_not_in_condition <- in_nhood_not_cond
  fisher_result$not_in_nhood_in_condition <- tot_cond - in_nhood_cond
  fisher_result$not_in_nhood_not_in_condition <- tot_not_cond - in_nhood_not_cond
  fisher_result$total_in_condition <- tot_cond
  fisher_result$total_not_in_condition <- tot_not_cond

  return(fisher_result)
}

run_nhood_wilcox <- function(nhood_counts, subj_info, da_variable, disease_group){
  # perform a one-sided Wilcoxon test on each neighborhood comparing neighborhood counts
  # normalized by subject depth, testing whether frequencies are higher in disease_group
  # nhood_counts = neighborhood x subject count matrix (i.e. milo@nhoodCounts). Should already include
  #                subjects absent from the current subset so totals reflect the WHOLE dataset
  # subj_info = one row per subject with subject IDs as row names, a da_variable column, and
  #             a depth column (i.e. lib_sizes)
  # disease_group = level of da_variable we are testing for enrichment

  subj_info <- subj_info[colnames(nhood_counts), , drop = FALSE]
  is_dis <- subj_info[[da_variable]] == disease_group

  # normalize each subject's neighborhood counts by their depth
  nhood_freqs <- sweep(as.matrix(nhood_counts), 2, subj_info$depth, '/')

  p_vals <- apply(nhood_freqs, 1, function(freqs){
    wilcox.test(freqs[!is_dis], freqs[is_dis], alternative = 'less', exact = FALSE)$p.value
  })

  return(list(wilcox_result = data.frame(PValue = unname(p_vals)),
              nhood_freqs = nhood_freqs))
}

get_spatial_fdr <- function(milo, pvalues, reduced_dim, weighting){
  # apply the same spatial FDR correction testNhoods uses to a new set of neighborhood p-values
  # pvalues = p-values in the same order as the neighborhoods (columns) in milo@nhoods

  graphSpatialFDR(x.nhoods = nhoods(milo),
                  graph = graph(milo),
                  weighting = weighting,
                  k = milo@.k,
                  pvalues = pvalues,
                  indices = nhoodIndex(milo),
                  distances = nhoodDistances(milo),
                  reduced.dimensions = reducedDim(milo, reduced_dim))
}

assign_min_p_nhood <- function(nhoods, nhood_results, fdr_col, keep_cols = c(), prefix = '', score_name = 'FDR'){
  # assign each sequence to the neighborhood with the lowest fdr_col value out of all
  # the neighborhoods it belongs to. Sequences in no neighborhood are returned with NA.
  # nhoods = sequence x neighborhood membership matrix (i.e. milo@nhoods). Columns must be
  #          in the same order as the rows of nhood_results
  # nhood_results = neighborhood-level test results containing nhood_id, fdr_col and keep_cols
  # keep_cols = named vector of additional columns to carry over from the assigned neighborhood,
  #             names are used for the output columns (i.e. c(PValue = 'fisher_PValue'))
  # prefix = prefix for output column names so results from multiple tests can be joined
  # score_name = output column name suffix for the fdr_col value (i.e. 'PValue' if assigning on raw p-values)

  memberships <- Matrix::summary(as(nhoods, 'CsparseMatrix'))
  memberships <- memberships[memberships$x != 0, ]

  # order by FDR, then neighborhood order, so ties go to the first neighborhood
  min_p <- data.frame(id_col = row.names(nhoods)[memberships$i],
                      nhood_idx = memberships$j) %>%
    dplyr::mutate(fdr = nhood_results[[fdr_col]][nhood_idx]) %>%
    dplyr::filter(!is.na(fdr)) %>%
    dplyr::arrange(fdr, nhood_idx) %>%
    dplyr::distinct(id_col, .keep_all = TRUE)

  min_p_df <- data.frame(id_col = min_p$id_col,
                         min_nhood_id = nhood_results$nhood_id[min_p$nhood_idx],
                         min_nhood_score = min_p$fdr)
  names(min_p_df)[names(min_p_df) == 'min_nhood_score'] <- paste0('min_nhood_', score_name)

  for (col in names(keep_cols)){
    min_p_df[[paste0('min_nhood_', col)]] <- nhood_results[[keep_cols[[col]]]][min_p$nhood_idx]
  }

  colnames(min_p_df)[-1] <- paste0(prefix, colnames(min_p_df)[-1])

  # add back sequences in no neighborhood
  min_p_df <- data.frame(id_col = row.names(nhoods)) %>%
    dplyr::left_join(min_p_df, by = 'id_col')

  return(min_p_df)
}

evaluate_results <- function(seq_table, p_val_col, auc_variable, name){
  # seq_table = table containing sequence-level information including auc_variable.
  #             Assumes p-values are already linked to sequences in this table
  # p_val_col = the column you want to use for p values to make AUC thresholds
  # auc_variable = variable to use for getting positives in the AUC curve
  # name = a name to specify for saving figures and tables (i.e. the name of the test)

  # get AUPRC baseline - fraction of positive events
  auprc_baseline <- mean(seq_table[[auc_variable]], na.rm = T)

  # sequences with no nhood get a min p of 1
  total_cells <- nrow(seq_table)
  invalid_cells <- sum(is.na(seq_table[[p_val_col]]))
  valid_cells <- total_cells - invalid_cells
  seq_table[is.na(seq_table[[p_val_col]]), p_val_col] <- 1

  # avoid floating point errors
  seq_table[[p_val_col]] <- round(seq_table[[p_val_col]], 6)
  auc_thresholds <- sort(unique(seq_table[[p_val_col]]))

  # add to the largest to make sure the entire curve is captured
  tot_thresh <- length(auc_thresholds)
  auc_thresholds[tot_thresh] <- auc_thresholds[tot_thresh] + 1e-3

  auc_data <- lapply(auc_thresholds, function(thresh){

    # get cells with min nhood p below threshold
    da.cell.list <- seq_table[[p_val_col]] < thresh

    true_pos <- sum(da.cell.list == T & seq_table[[auc_variable]] == T)
    false_neg <- sum(da.cell.list == F & seq_table[[auc_variable]] == T)
    true_neg <- sum(da.cell.list == F & seq_table[[auc_variable]] == F)
    false_pos <- sum(da.cell.list == T & seq_table[[auc_variable]] == F)

    return(data.frame('TPR' = true_pos / (true_pos + false_neg),
                      'FPR' = 1 - (true_neg / (true_neg + false_pos)),
                      'Precision' = true_pos / (false_pos + true_pos),
                      'FDR' = false_pos / (false_pos + true_pos),
                      'TP' = true_pos,
                      'FP' = false_pos,
                      'TN' = true_neg,
                      'FN' = false_neg
                      ))

  })

  auc_df <- do.call(rbind, auc_data)
  auc_df$threshold <- auc_thresholds

  # estimate the first precision point - it should always be NA b/c no false or true positives below the first threshold
  if (is.na(auc_df[1,'Precision']) & nrow(auc_df) > 1){
    auc_df[1,'Precision'] <- auc_df[2,'Precision']
  }

  write.table(auc_df,
              file.path(OUTPUT_DIR, 'tables', paste0('evaluation_curve_vals_', name, '.tsv')),
              sep = '\t', row.names = F, quote = F)

  # get auroc
  auroc <- pracma::trapz(auc_df$FPR, auc_df$TPR)

  # get auprc
  auprc <- pracma::trapz(auc_df$TPR, auc_df$Precision)

  # get average precision - precision at each threshold weighted by the increase in recall
  ap <- sum(diff(c(0, auc_df$TPR)) * auc_df$Precision, na.rm = T)

  pretty_name <- stringr::str_replace_all(name, '_', ' ')

  auc_df %>%
    ggplot(aes(x = FPR, y = TPR)) +
    geom_abline(slope = 1, intercept = 0, color = 'gray') +
    geom_point() +
    geom_line() +
    labs(title = paste0(pretty_name, ' threshold ', round(min(auc_thresholds)), ' to ', round(max(auc_thresholds), 3)),
         subtitle = paste0(pretty_name, ' AUROC: ', round(auroc, 3), '; ',
                           prettyNum(valid_cells, big.mark = ",", scientific = FALSE), '/',
                           prettyNum(total_cells, big.mark = ",", scientific = FALSE), ' cells in DA neighborhoods')) +
    theme_minimal()

  ggsave(file.path(OUTPUT_DIR, 'figures', paste0('AUROC_', name, '.png')),
         device = 'png',
         width = 7,
         height = 6)

  # PRC
  auc_df %>%
    ggplot(aes(x = TPR, y = Precision)) +
    geom_hline(yintercept = auprc_baseline, color = 'red', linetype = 'dashed') +
    geom_point() +
    geom_line() +
    labs(title = paste0(pretty_name, ' threshold ', round(min(auc_thresholds)), ' to ', round(max(auc_thresholds), 3)),
         subtitle = paste0(pretty_name, ' AUPRC: ', round(auprc, 3), '; AP: ', round(ap, 3), '; ',
                           prettyNum(valid_cells, big.mark = ",", scientific = FALSE), '/',
                           prettyNum(total_cells, big.mark = ",", scientific = FALSE), ' cells in DA neighborhoods'),
         x = 'Recall') +
    theme_minimal() +
    scale_y_continuous(limits = c(0, 1))

  ggsave(file.path(OUTPUT_DIR, 'figures', paste0('AUPRC_', name, '.png')),
         device = 'png',
         width = 7,
         height = 6)

  return(list('AUROC' = auroc,
              'AUPRC' = auprc,
              'average_precision' = ap))
}

calc_FDR <- function(results_table, p_val_col, auc_variable, alpha){
  # results table = table containing some kind of test result (Fisher Exact, Wilcox, etc.).
  #                 If rows are sequences, sequence-level FDR is calculated.
  #                 If rows are clusters, cluster-level FDR is calculated.
  # p_val_col = the column you want to apply the p-value or FDR threshold to
  # auc_variable = variable to use for getting positives

  # get all significant
  results_filtered <- results_table %>%
    dplyr::filter(!is.na(!!sym(p_val_col))) %>%
    dplyr::filter(!!sym(p_val_col) < alpha)

  # mean = TP / (TP + FP). We want FP / (TP + FP), which is 1 - (TP/(TP+FP))
  FDR <- 1 - mean(results_filtered[[auc_variable]])

  return(FDR)

}

########################
### PREP ENVIRONMENT ###
########################

# prepare to take input parameters
parser <- ArgumentParser(description = "Data location and Milo algorithm hyperparameters.")

parser$add_argument('-d', '--data_loc', type = 'character', default = 'data',
                    help = 'File path for the embedding or RNA-Seq data location.')

parser$add_argument('-md', '--metadata_loc', type = 'character', default = 'metadata',
                    help = 'File path for the metadata location. Metadata and data files should have 1:1 matching sequence identifiers.')

parser$add_argument('-li', '--library_sizes', type = 'character', default = NULL,
                    help = 'File path for the library sizes file location.')

parser$add_argument('-o', '--output_dir', type = 'character', default = 'DAseq_output',
                    help = 'Specify an output directory location.')

parser$add_argument('-da', '--da_variable', type = 'character', default = 'status',
                    help = 'Stratification variable that should be used to determine for differential abundance. There should be two levels in this factor/categorical variable.')

parser$add_argument('-dg', '--disease_group', type = 'character', default = 'disease',
                    help = 'Level of the DA variable to test for enrichment in the Fisher Exact and Wilcoxon tests.')

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

################################################################################

# Parse the arguments
args <- parser$parse_args()

# specify which dataset we are analyzing
DATA_LOC <- args$data_loc
MD_LOC <- args$metadata_loc
LIB_SIZES_LOC <- args$library_sizes
MD_NAME <- stringr::str_split_i(basename(MD_LOC), '_md', 1)

OUTPUT_DIR <- args$output_dir

message(paste0('Data will be saved to ', OUTPUT_DIR, '.'))

# get the variable for DA calculations
DA_VAR <- args$da_variable
DISEASE_GP <- args$disease_group

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
  message(paste0('AUC variable ', AUC_VAR, ' will be used.'))
} else{
  message('AUC will not be calculated.')
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

# NOTE: not using GLMM currently, but could be implemented if needed
formula_string <- paste0('~ ', DA_VAR)

design_formula <- as.formula(formula_string)

message(paste0('Using formula: ', formula_string))

# get subjects from whole dataset who may be missing in subset i.e. ASC
if(!is.null(LIB_SIZES_LOC)){
  lib_sizes <- read.csv(LIB_SIZES_LOC, sep = '\t') %>% as.data.frame()
  
  # find missing subjects
  absent_subj <- setdiff(lib_sizes$subject_id, colnames(milo@nhoodCounts))
  
  if (length(absent_subj) > 0){
    message(paste0('Adding missing subjects ', paste(absent_subj, collapse = ', '), ' to Milo neighborhood counts.'))
  
    # add to counts
    new_cols <- Matrix(0, nrow = nrow(milo@nhoodCounts), ncol = length(absent_subj), sparse = TRUE)
    colnames(new_cols) <- absent_subj
    
    new_nhood_counts <- cbind(milo@nhoodCounts, new_cols)
    
    milo@nhoodCounts <- new_nhood_counts
  } else{
    message('No missing subjects detected.')
  }
  
  # make sure library size df is consistent with Milo object
  lib_sizes <- as.data.frame(lib_sizes)
  row.names(lib_sizes) <- lib_sizes$subject_id
  
  lib_sizes <- lib_sizes[colnames(milo@nhoodCounts), , drop=FALSE]
  
  # sanity check
  subj_match <- all(row.names(lib_sizes) == colnames(milo@nhoodCounts))
  if (subj_match == F){
    warning('Subjects in Milo neighborhood count matrix are not in line with library size summary. Regression results will not be accurate.')
  }
  
} else{
  message('Calculating library sizes from metadata...')
  lib_sizes <- md %>%
    dplyr::group_by(subject_id, !!sym(DA_VAR)) %>%
    dplyr::summarize(depth = n(), .groups = "drop_last") %>%
    as.data.frame()
  
  row.names(lib_sizes) <- lib_sizes$subject_id
}

print('Library sizes to be used:')
print(lib_sizes)

# input custom cell.sizes for either the dataset as is OR 
# the entire dataset, even if in ASC mode
cell.sizes <- lib_sizes$depth
names(cell.sizes) <- row.names(lib_sizes)

# make design_df based on lib sizes
design_df <- lib_sizes
design_df$sample_id <- design_df$subject_id
design_df <- design_df %>% dplyr::select(-depth)

## Reorder rownames to match columns of nhoodCounts(milo) - should already match though
design_df <- design_df[colnames(nhoodCounts(milo)), , drop=FALSE]

design_df$sample_id <- as.factor(design_df$sample_id)
design_df$subject_id <- as.factor(design_df$subject_id)
design_df[,DA_VAR] <- as.factor(design_df[,DA_VAR])

print('Design:')
print(design_df)

print(table(data.frame(colData(milo))$subject_id))

da_results <- testNhoods(milo, 
                         cell.sizes = cell.sizes,
                         norm.method = 'logMS',
                         design = design_formula, 
                         design.df = design_df,
                         reduced.dim = 'embedding',
                         fdr.weighting = 'neighbour-distance')

# }

######################################
# FISHER EXACT AND WILCOXON TESTS    #
######################################
# uses nhood counts and library sizes, which include all subjects in the WHOLE dataset
if (!DISEASE_GP %in% lib_sizes[[DA_VAR]]){
  stop(paste0('Disease group ', DISEASE_GP, ' not found in ', DA_VAR, ' column.'))
}

if (any(is.na(lib_sizes[colnames(milo@nhoodCounts), 'depth']))){
  stop('Subjects in Milo neighborhood count matrix are missing from library sizes. Fisher and Wilcoxon tests cannot be run.')
}

if (!all(as.numeric(row.names(milo@nhoodCounts)) == da_results$Nhood)){
  warning('Neighborhood count matrix rows do not match GLM results. Fisher and Wilcoxon results will not be accurate.')
}

message(paste0('Running Fisher Exact tests for ', DISEASE_GP, ' enrichment...'))
fisher_res <- run_nhood_fisher(milo@nhoodCounts, lib_sizes, DA_VAR, DISEASE_GP)

da_results$fisher_PValue <- fisher_res$PValue
da_results$fisher_odds_ratio <- fisher_res$odds_ratio
da_results$fisher_SpatialFDR <- get_spatial_fdr(milo, fisher_res$PValue, 'embedding', 'neighbour-distance')
da_results$fisher_BH <- p.adjust(fisher_res$PValue, method = 'BH')

message(paste0('Running one-sided Wilcoxon tests for ', DISEASE_GP, ' enrichment...'))
wilcox_res <- run_nhood_wilcox(milo@nhoodCounts, lib_sizes, DA_VAR, DISEASE_GP)

da_results$wilcox_PValue <- wilcox_res$wilcox_result$PValue
da_results$wilcox_SpatialFDR <- get_spatial_fdr(milo, wilcox_res$wilcox_result$PValue, 'embedding', 'neighbour-distance')
da_results$wilcox_BH <- p.adjust(wilcox_res$wilcox_result$PValue, method = 'BH')

# save full test info
fisher_res$nhood_id <- as.character(unlist(milo@nhoodIndex))
write.table(fisher_res,
            file.path(OUTPUT_DIR, 'tables', 'fisher_results.tsv'),
            sep = '\t', row.names = F, quote = F)

nhood_freqs <- data.frame(nhood_id = as.character(unlist(milo@nhoodIndex)),
                          wilcox_res$nhood_freqs, check.names = F)
write.table(nhood_freqs,
            file.path(OUTPUT_DIR, 'tables', 'nhood_subject_freqs.tsv'),
            sep = '\t', row.names = F, quote = F)

print('Top DA results:')
da_results %>%
  arrange(SpatialFDR) %>%
  head() %>%
  print()

########### ADDED DIAGNOSTIC VISUALS ########### 
for (res in c('PValue', 'SpatialFDR', 'fisher_PValue', 'fisher_SpatialFDR',
              'fisher_BH', 'wilcox_PValue', 'wilcox_SpatialFDR', 'wilcox_BH')){

  da_results %>%
    ggplot(aes(x = !!sym(res))) + 
    geom_histogram(color = 'white', binwidth = 0.01) + 
    theme_bw() +
    labs(title = paste0('Milo ', res, ' Distribution')) +
    coord_cartesian(xlim = c(0, 1))

  ggsave(file.path(OUTPUT_DIR, 'figures', paste0(res, '_hist.png')),
        device = 'png', width = 8, height = 6, units = 'in')
}


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


nhoods_match <- all(colnames(milo@nhoods) == da_results$nhood_id)
if (!nhoods_match){
  message('WARNING: neighborhood matrix does not match results table. AUC values will not be accurate.')
}


# match each cell with the lowest p-value of all the neighborhoods it occupies
glm_min_p <- assign_min_p_nhood(milo@nhoods, da_results, 'SpatialFDR', 
                                keep_cols = c(PValue = 'PValue', logFC = 'logFC'))

fisher_min_p <- assign_min_p_nhood(milo@nhoods, da_results, 'fisher_SpatialFDR', 
                                   keep_cols = c(PValue = 'fisher_PValue', BH = 'fisher_BH', 
                                                 odds_ratio = 'fisher_odds_ratio'),
                                   prefix = 'fisher_')

# the Fisher test is very conservative, so also assign on raw p-values for AUCs
fisher_raw_min_p <- assign_min_p_nhood(milo@nhoods, da_results, 'fisher_PValue', 
                                       keep_cols = c(SpatialFDR = 'fisher_SpatialFDR', BH = 'fisher_BH', 
                                                     odds_ratio = 'fisher_odds_ratio'),
                                       prefix = 'fisher_raw_', score_name = 'PValue')

wilcox_min_p <- assign_min_p_nhood(milo@nhoods, da_results, 'wilcox_SpatialFDR', 
                                   keep_cols = c(PValue = 'wilcox_PValue', BH = 'wilcox_BH'),
                                   prefix = 'wilcox_')

min_p_nhoods_df <- glm_min_p %>%
  dplyr::left_join(fisher_min_p, by = 'id_col') %>%
  dplyr::left_join(fisher_raw_min_p, by = 'id_col') %>%
  dplyr::left_join(wilcox_min_p, by = 'id_col')

if (AUC_VAR != FALSE){
  min_p_nhoods_df <- min_p_nhoods_df %>%
    dplyr::left_join(md_reduced[c('id_col', AUC_VAR)], by = 'id_col')
}

write.table(min_p_nhoods_df, 
            file.path(OUTPUT_DIR, 'tables', 'min_p_nhoods.tsv'),
            sep = '\t', row.names = F, quote = F)

# create a version with ALL results included (so we can reference all the neighborhoods)
min_p_nhoods_df_merge <- da_results %>%
  full_join(min_p_nhoods_df, by = join_by(nhood_id == min_nhood_id), relationship = 'one-to-many', na_matches = 'never')

write.table(min_p_nhoods_df_merge, 
            file.path(OUTPUT_DIR, 'tables', paste0(MD_NAME, '_seq_results.tsv')), 
            sep = '\t', row.names = F, quote = F)

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
  # add AUC var info
  if (!AUC_VAR %in% colnames(min_p_nhoods_df)){
    min_p_nhoods_df <- min_p_nhoods_df %>%
      dplyr::left_join(md_reduced[c('id_col', AUC_VAR)])
  }

  message('Making AUC curves...')
  glm_eval <- evaluate_results(min_p_nhoods_df, 'min_nhood_FDR', AUC_VAR, 'GLM')
  fisher_eval <- evaluate_results(min_p_nhoods_df, 'fisher_raw_min_nhood_PValue', AUC_VAR, 'Fisher')
  wilcox_eval <- evaluate_results(min_p_nhoods_df, 'wilcox_min_nhood_FDR', AUC_VAR, 'One-Sided_Wilcoxon')

  # FDR
  FDR <- calc_FDR(min_p_nhoods_df, 'min_nhood_FDR', AUC_VAR, 0.05)
  fisher_FDR <- calc_FDR(min_p_nhoods_df, 'fisher_min_nhood_FDR', AUC_VAR, 0.05)
  wilcox_FDR <- calc_FDR(min_p_nhoods_df, 'wilcox_min_nhood_FDR', AUC_VAR, 0.05)

  # change the sequences with no nhood to a min p of 1 for Jaccard calculations
  min_p_nhoods_df[is.na(min_p_nhoods_df$min_nhood_FDR), 'min_nhood_FDR'] <- 1
  min_p_nhoods_df$min_nhood_FDR <- round(min_p_nhoods_df$min_nhood_FDR, 6)
  
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
stat_table <- data.frame('tool' = c('Milo', 'Milo + Fisher', 'Milo + One-sided Wilcoxon'),
                         'total_seqs' = c(ncol(milo)),
                         'total_subj' = ncol(milo@nhoodCounts),
                         'time (min)' = as.numeric(time_taken, units = "mins"),
                         'subjects' = paste(names(table(milo@colData$subject_id)), collapse = ', '),
                         'depths' = paste(table(milo@colData$subject_id), collapse = ', '),
                         check.names = F)

if (AUC_VAR != FALSE){
  
  purity_stats <- subj_nhood_cts %>%
    dplyr::filter(hit_seqs > 0) %>%
    dplyr::select(-c('subject_id', 'cells_per_subj', 'pct_subj')) %>%
    distinct()

  # document "purity" of clusters with simulated sequences visually
  make_purity_plot(purity_stats, 'nhood_id', 'pct_hits', 'cells_per_nhood', AUC_VAR)

  stat_table$num_hit_clusters <- nrow(purity_stats)
  stat_table$avg_pct_hits <- mean(purity_stats$pct_hits)
  stat_table$tot_hits <- c(sum(milo@colData[[AUC_VAR]], na.rm = T)) 
  stat_table$pct_hits <- c(mean(milo@colData[[AUC_VAR]], na.rm = T) * 100)
  stat_table$Jaccard_0.005 = c(jaccard_005, NA, NA)
  stat_table$Jaccard_0.05 = c(jaccard_05, NA, NA)
  stat_table$Jaccard_0.1 = c(jaccard_1, NA, NA)
  stat_table$Jaccard_max = c(Jaccard_max, NA, NA)
  stat_table$Jaccard_max_p = c(Jaccard_max_p, NA, NA)
  stat_table$AUROC <- c(glm_eval$AUROC, fisher_eval$AUROC, wilcox_eval$AUROC)
  stat_table$AUPRC <- c(glm_eval$AUPRC, fisher_eval$AUPRC, wilcox_eval$AUPRC)
  stat_table$average_precision <- c(glm_eval$average_precision, fisher_eval$average_precision, wilcox_eval$average_precision)
  stat_table$FDR <- c(FDR, fisher_FDR, wilcox_FDR)
  
  stat_table <- stat_table[c('tool', 'total_seqs', 'total_subj', 'tot_hits', 'pct_hits',
                             'num_hit_clusters', 'avg_pct_hits',
                             'AUROC', 'AUPRC', 'average_precision', 'FDR', 
                             'Jaccard_0.005', 'Jaccard_0.05',
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

# make viz for all cells in neighborhoods significant at alpha 0.05
sig_nhoods <- da_results %>%
  dplyr::filter(SpatialFDR < 0.05) %>%
  dplyr::pull(nhood_id)

sig_nhood_cells <- milo@nhoods[,as.character(sig_nhoods)]

if(length(sig_nhoods) == 1){
  sig_nhood_idx <- sig_nhood_cells > 0
  sig_nhood_cell_ids <- names(sig_nhood_cells[sig_nhood_idx])
} else{
  sig_nhood_idx <- rowSums(sig_nhood_cells) > 0
  sig_nhood_cell_ids <- row.names(sig_nhood_cells[sig_nhood_idx,])
}

umap_coords <- umap_coords %>%
  dplyr::mutate(da_cell = if_else(id_col %in% sig_nhood_cell_ids, TRUE, FALSE))

if (nrow(umap_coords) < 400000){
  message('Visualize significant neighborhoods.')
  make_UMAP_viz('da_cell', 'in significant neighborhood \n (spatial FDR < 0.05)', custom_pal = c('TRUE' = "red", 'FALSE' = "gray"))
}

message(paste0('Ending run: ', Sys.time()))

sessionInfo()
