#!/usr/bin/env Rscript
message(paste0('Starting run: ', Sys.time()))

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(Matrix)
  library(pracma)
  library(reticulate)
})

# Evaluates the output of a convergence tool. Every tool provides, for each group
# (an ASC in asc_mode, or 'all' otherwise):
#   1. a cluster x subject count matrix (cluster IDs already prefixed with the group ID)
#   2. a sequence-to-cluster membership table (sequence_id, subject_id, group_id, cluster_id,
#      and the AUC variable if available). Sequences can belong to multiple clusters (i.e. Milo
#      neighborhoods), or none (cluster_id = NA)
#   3. a tool info table (group_id, time_min)
#   4. Milo only: the Milo object, used for the negative binomial GLM and spatial FDR
# along with the library sizes from preprocessing, which give the depth and DA variable
# of every subject in the WHOLE dataset.

########################
### HELPER FUNCTIONS ###
########################

read_count_matrix <- function(count_files, subjects){
  # read the cluster x subject count matrices from all groups and combine into one matrix
  # with a column for every subject in subjects. Subjects absent from a group get 0 counts.

  count_list <- lapply(count_files, function(f){
    read.delim(f, check.names = F, colClasses = c(cluster_id = 'character'))
  })

  counts <- dplyr::bind_rows(count_list)

  if (any(duplicated(counts$cluster_id))){
    stop('Duplicate cluster IDs found across count matrices. Cluster IDs must be prefixed with their group ID.')
  }

  unknown_subj <- setdiff(colnames(counts), c('cluster_id', subjects))
  if (length(unknown_subj) > 0){
    stop(paste0('Subjects in count matrices are missing from library sizes: ', paste(unknown_subj, collapse = ', ')))
  }

  absent_subj <- setdiff(subjects, colnames(counts))
  if (length(absent_subj) > 0){
    message(paste0('Adding missing subjects ', paste(absent_subj, collapse = ', '), ' to cluster count matrix.'))
    counts[absent_subj] <- 0
  }

  count_mat <- as.matrix(counts[subjects])
  count_mat[is.na(count_mat)] <- 0
  row.names(count_mat) <- counts$cluster_id

  return(count_mat)
}

run_fisher <- function(count_mat, subj_info, da_variable, disease_group){
  # perform a one-sided Fisher Exact Test on each cluster based on SUBJECTS in the cluster,
  # testing whether subjects in disease_group are over-represented
  # count_mat = cluster x subject count matrix including all subjects in the WHOLE dataset
  # subj_info = one row per subject with subject IDs as row names and a da_variable column (i.e. lib_sizes)
  # disease_group = level of da_variable we are testing for enrichment

  subj_info <- subj_info[colnames(count_mat), , drop = FALSE]
  is_dis <- subj_info[[da_variable]] == disease_group

  tot_cond <- sum(is_dis)
  tot_not_cond <- sum(!is_dis)

  # count subjects with at least one sequence in each cluster, with and without condition
  present <- count_mat > 0
  in_cluster_cond <- unname(rowSums(present[, is_dis, drop = FALSE]))
  in_cluster_not_cond <- unname(rowSums(present[, !is_dis, drop = FALSE]))

  # there are only a few possible contingency tables, so only test each one once
  tables <- unique(data.frame(in_cluster_cond, in_cluster_not_cond))

  tables_res <- lapply(seq_len(nrow(tables)), function(i){

    # build contingency table to test for a CONDITION cluster
    #
    #                  cluster
    #                No    Yes
    #               ___________
    #            No|     |     |
    # condition    |_____|_____|
    #           Yes|     |     |
    #              |_____|_____|

    contingency_table <- matrix(c(tot_not_cond - tables$in_cluster_not_cond[i], tot_cond - tables$in_cluster_cond[i],
                                  tables$in_cluster_not_cond[i], tables$in_cluster_cond[i]), 2, 2)

    fisher <- fisher.test(contingency_table, alternative = 'greater')

    return(data.frame(p_value = fisher$p.value,
                      odds_ratio = unname(fisher$estimate)))
  })

  tables <- cbind(tables, do.call(rbind, tables_res))

  fisher_result <- data.frame(cluster_id = row.names(count_mat),
                              in_cluster_cond = in_cluster_cond,
                              in_cluster_not_cond = in_cluster_not_cond) %>%
    dplyr::left_join(tables, by = c('in_cluster_cond', 'in_cluster_not_cond')) %>%
    dplyr::transmute(cluster_id,
                     p_value,
                     odds_ratio,
                     in_cluster_in_condition = in_cluster_cond,
                     in_cluster_not_in_condition = in_cluster_not_cond,
                     not_in_cluster_in_condition = tot_cond - in_cluster_cond,
                     not_in_cluster_not_in_condition = tot_not_cond - in_cluster_not_cond,
                     total_in_condition = tot_cond,
                     total_not_in_condition = tot_not_cond)

  return(fisher_result)
}

get_cluster_freqs <- function(count_mat, subj_info){
  # normalize each subject's cluster counts by their depth
  subj_info <- subj_info[colnames(count_mat), , drop = FALSE]
  sweep(count_mat, 2, subj_info$depth, '/')
}

run_wilcox <- function(cluster_freqs, subj_info, da_variable, disease_group){
  # perform a one-sided Wilcoxon test on each cluster comparing cluster counts
  # normalized by subject depth, testing whether frequencies are higher in disease_group
  # cluster_freqs = cluster x subject normalized frequency matrix including all subjects in the WHOLE dataset
  # subj_info = one row per subject with subject IDs as row names and a da_variable column (i.e. lib_sizes)
  # disease_group = level of da_variable we are testing for enrichment

  subj_info <- subj_info[colnames(cluster_freqs), , drop = FALSE]
  is_dis <- subj_info[[da_variable]] == disease_group

  # many clusters (i.e. singletons) share the same frequencies, so only test each pattern once
  freq_keys <- apply(cluster_freqs, 1, paste, collapse = ',')
  unique_rows <- !duplicated(freq_keys)

  p_vals <- apply(cluster_freqs[unique_rows, , drop = FALSE], 1, function(freqs){
    wilcox.test(freqs[!is_dis], freqs[is_dis], alternative = 'less', exact = FALSE)$p.value
  })
  names(p_vals) <- freq_keys[unique_rows]

  return(data.frame(cluster_id = row.names(cluster_freqs),
                    p_value = unname(p_vals[freq_keys])))
}

run_beta_binomial <- function(count_mat, subj_info, da_variable, disease_group, control_group){
  # run the dejabu beta-binomial test on the cluster x subject count matrix, testing whether
  # clusters are enriched in disease_group compared to control_group
  # control_group = level of da_variable used as the reference

  dejabu_counts <- import('dejabu.Counts')
  dejabu_enrichment <- import('dejabu.Enrichment')

  # pass data frames so the cluster and subject names come through as the pandas index/columns
  clusters <- dejabu_counts$ClusterCounts(r_to_py(as.data.frame(count_mat)),
                                          r_to_py(subj_info),
                                          group_col = da_variable,
                                          subject_col = 'subject_id')

  contrast <- clusters$contrast(disease_group, control_group)
  bb_res <- dejabu_enrichment$BetaBinomialTest()$run(contrast)

  res_df <- py_to_r(bb_res$toFrame())
  res_df$cluster_id <- row.names(res_df)

  return(res_df)
}

run_overdispersed_z <- function(count_mat, subj_info, da_variable, disease_group, control_group){
  # run the dejabu overdispersed z-test on the cluster x subject count matrix, testing whether
  # clusters are enriched in disease_group compared to control_group
  # control_group = level of da_variable used as the reference

  dejabu_counts <- import('dejabu.Counts')
  dejabu_enrichment <- import('dejabu.Enrichment')

  # pass data frames so the cluster and subject names come through as the pandas index/columns
  clusters <- dejabu_counts$ClusterCounts(r_to_py(as.data.frame(count_mat)),
                                          r_to_py(subj_info),
                                          group_col = da_variable,
                                          subject_col = 'subject_id')

  contrast <- clusters$contrast(disease_group, control_group)
  oz_res <- dejabu_enrichment$OverdispersedZTest()$run(contrast)

  res_df <- py_to_r(oz_res$toFrame())
  res_df$cluster_id <- row.names(res_df)

  return(res_df)
}

run_pluginz <- function(count_mat, subj_info, da_variable, disease_group, control_group){
  # run the dejabu beta-binomial test on the cluster x subject count matrix, testing whether
  # clusters are enriched in disease_group compared to control_group
  # control_group = level of da_variable used as the reference

  dejabu_counts <- import('dejabu.Counts')
  dejabu_enrichment <- import('dejabu.Enrichment')

  # pass data frames so the cluster and subject names come through as the pandas index/columns
  clusters <- dejabu_counts$ClusterCounts(r_to_py(as.data.frame(count_mat)),
                                          r_to_py(subj_info),
                                          group_col = da_variable,
                                          subject_col = 'subject_id')

  contrast <- clusters$contrast(disease_group, control_group)
  pz_res <- dejabu_enrichment$PluginZTest(estimator = 'subject', combination = 'averaged')$run(contrast)

  res_df <- py_to_r(pz_res$toFrame())
  res_df$cluster_id <- row.names(res_df)

  return(res_df)
}

run_milo_glm <- function(milo, cluster_ids, subj_info, da_variable){
  # run the Milo negative binomial GLM on one group's Milo object, using the library
  # sizes of the WHOLE dataset (subjects absent from the group are added with 0 counts)
  # cluster_ids = cluster IDs of the neighborhoods, in the same order as nhoodCounts(milo)

  absent_subj <- setdiff(row.names(subj_info), colnames(nhoodCounts(milo)))

  if (length(absent_subj) > 0){
    message(paste0('Adding missing subjects ', paste(absent_subj, collapse = ', '), ' to Milo neighborhood counts.'))

    new_cols <- Matrix(0, nrow = nrow(nhoodCounts(milo)), ncol = length(absent_subj), sparse = TRUE)
    colnames(new_cols) <- absent_subj

    milo@nhoodCounts <- cbind(milo@nhoodCounts, new_cols)
  }

  subj_info <- subj_info[colnames(nhoodCounts(milo)), , drop = FALSE]

  # input custom cell.sizes for the entire dataset, even if in ASC mode
  cell.sizes <- subj_info$depth
  names(cell.sizes) <- row.names(subj_info)

  design_df <- data.frame(sample_id = as.factor(row.names(subj_info)),
                          row.names = row.names(subj_info))
  design_df[[da_variable]] <- as.factor(subj_info[[da_variable]])

  design_formula <- as.formula(paste0('~ ', da_variable))
  message(paste0('Using formula: ~ ', da_variable))

  da_results <- testNhoods(milo,
                           cell.sizes = cell.sizes,
                           norm.method = 'logMS',
                           design = design_formula,
                           design.df = design_df,
                           reduced.dim = 'embedding',
                           fdr.weighting = 'neighbour-distance')

  if (!all(da_results$Nhood == seq_len(nrow(nhoodCounts(milo))))){
    stop('Milo GLM results are not in the same order as the neighborhood counts.')
  }

  return(data.frame(cluster_id = cluster_ids,
                    p_value = da_results$PValue,
                    spatial_fdr = da_results$SpatialFDR,
                    logFC = da_results$logFC,
                    logCPM = da_results$logCPM))
}

get_spatial_fdr <- function(milo, pvalues){
  # apply the same spatial FDR correction testNhoods uses to a new set of neighborhood p-values
  # pvalues = p-values in the same order as the neighborhoods (columns) in milo@nhoods

  graphSpatialFDR(x.nhoods = nhoods(milo),
                  graph = graph(milo),
                  weighting = 'neighbour-distance',
                  k = milo@.k,
                  pvalues = pvalues,
                  indices = nhoodIndex(milo),
                  distances = nhoodDistances(milo),
                  reduced.dimensions = reducedDim(milo, 'embedding'))
}

score_sequences <- function(memberships, test_results, score_col){
  # assign each sequence the cluster with the lowest score_col value out of all the clusters
  # it belongs to. Ties go to the cluster listed first in the membership table.
  # Sequences in no cluster (or only untested clusters) are returned with NA.

  assigned <- memberships %>%
    dplyr::filter(!is.na(cluster_id)) %>%
    dplyr::inner_join(test_results[c('cluster_id', score_col, 'p_adj')], by = 'cluster_id') %>%
    dplyr::rename(score = !!sym(score_col)) %>%
    dplyr::filter(!is.na(score)) %>%
    dplyr::arrange(score, membership_order) %>%
    dplyr::distinct(sequence_id, .keep_all = TRUE) %>%
    dplyr::select(sequence_id, assigned_cluster_id = cluster_id, score, p_adj)

  seq_info <- memberships %>%
    dplyr::select(-cluster_id, -membership_order) %>%
    dplyr::distinct(sequence_id, .keep_all = TRUE)

  return(seq_info %>% dplyr::left_join(assigned, by = 'sequence_id'))
}

evaluate_curve <- function(scores, truth, log_digits = 8){
  # get the ROC and precision-recall curves for sequence-level scores, where a sequence is
  # called DA at a threshold if its score is below the threshold
  # scores = sequence-level scores (i.e. p-values); sequences with no score get a score of 1
  # truth = sequence-level ground truth (i.e. simulated)
  # log_digits = decimal places kept for log10(scores) when finding ties

  scores[is.na(scores)] <- 1
  score_range <- range(scores)

  # avoid floating point errors by rounding on the log scale. This is a relative tolerance,
  # so scores within floating point error of each other are tied while very small p-values
  # keep their full resolution (p-values of 0 become -Inf, which still sort first)
  scores <- round(log10(scores), log_digits)
  truth <- as.logical(truth)

  thresholds <- sort(unique(scores))

  # count positives and negatives at each unique score
  score_idx <- match(scores, thresholds)
  pos <- tabulate(score_idx[truth %in% TRUE], nbins = length(thresholds))
  neg <- tabulate(score_idx[truth %in% FALSE], nbins = length(thresholds))

  # sequences called DA at each threshold are those with scores strictly below it.
  # Add to the largest (on the log scale) to make sure the entire curve is captured
  n_thresh <- length(thresholds)
  thresholds[n_thresh] <- thresholds[n_thresh] + 1e-3

  true_pos <- c(0, cumsum(pos))[seq_len(n_thresh)]
  false_pos <- c(0, cumsum(neg))[seq_len(n_thresh)]
  true_pos[n_thresh] <- sum(pos)
  false_pos[n_thresh] <- sum(neg)

  false_neg <- sum(pos) - true_pos
  true_neg <- sum(neg) - false_pos

  auc_df <- data.frame('TPR' = true_pos / (true_pos + false_neg),
                       'FPR' = 1 - (true_neg / (true_neg + false_pos)),
                       'Precision' = true_pos / (false_pos + true_pos),
                       'FDR' = false_pos / (false_pos + true_pos),
                       'TP' = true_pos,
                       'FP' = false_pos,
                       'TN' = true_neg,
                       'FN' = false_neg,
                       'threshold' = 10^thresholds,
                       'log10_threshold' = thresholds)

  # estimate the first precision point - it should always be NA b/c no false or true positives below the first threshold
  if (is.na(auc_df[1,'Precision']) & nrow(auc_df) > 1){
    auc_df[1,'Precision'] <- auc_df[2,'Precision']
  }

  return(list('table' = auc_df,
              'AUROC' = pracma::trapz(auc_df$FPR, auc_df$TPR),
              'AUPRC' = pracma::trapz(auc_df$TPR, auc_df$Precision),
              # average precision - precision at each threshold weighted by the increase in recall
              'average_precision' = sum(diff(c(0, auc_df$TPR)) * auc_df$Precision, na.rm = T),
              'baseline' = mean(truth, na.rm = T),
              'score_range' = score_range))
}

plot_curves <- function(curve, name){

  # range of the observed scores (not the bumped final threshold)
  title <- paste0('Scores ', signif(curve$score_range[1], 3), ' to ', signif(curve$score_range[2], 3))

  curve$table %>%
    ggplot(aes(x = FPR, y = TPR)) +
    geom_abline(slope = 1, intercept = 0, color = 'gray') +
    geom_point() +
    geom_line() +
    labs(title = title,
         subtitle = paste0(name, ' AUROC: ', round(curve$AUROC, 3))) +
    theme_minimal()

  ggsave(file.path(OUTPUT_DIR, 'figures', paste0('AUROC_', name, '.png')),
         device = 'png', width = 7, height = 6)

  curve$table %>%
    ggplot(aes(x = TPR, y = Precision)) +
    geom_hline(yintercept = curve$baseline, color = 'red', linetype = 'dashed') +
    geom_point() +
    geom_line() +
    labs(title = title,
         subtitle = paste0(name, ' AUPRC: ', round(curve$AUPRC, 3), '; AP: ', round(curve$average_precision, 3)),
         x = 'Recall') +
    theme_minimal() +
    scale_y_continuous(limits = c(0, 1))

  ggsave(file.path(OUTPUT_DIR, 'figures', paste0('AUPRC_', name, '.png')),
         device = 'png', width = 7, height = 6)
}

calc_FDR <- function(seq_scores, auc_variable, alpha){
  # sequence-level FDR among sequences assigned to a cluster with an adjusted p-value below alpha

  sig <- seq_scores %>%
    dplyr::filter(!is.na(p_adj), p_adj < alpha)

  # mean = TP / (TP + FP). We want FP / (TP + FP), which is 1 - (TP/(TP+FP))
  return(1 - mean(as.logical(sig[[auc_variable]])))
}

get_purity_stats <- function(memberships, auc_variable){
  # "purity" of clusters containing at least one hit sequence

  memberships %>%
    dplyr::filter(!is.na(cluster_id)) %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarize(total_cluster_seqs = n(),
                     hit_seqs = sum(as.logical(!!sym(auc_variable)), na.rm = T),
                     .groups = 'drop') %>%
    dplyr::mutate(pct_hits = hit_seqs / total_cluster_seqs) %>%
    dplyr::filter(hit_seqs > 0)
}

make_purity_plot <- function(purity_data, auc_variable){

  ggplot(purity_data, aes(x = cluster_id, y = pct_hits)) +
    geom_col(fill = 'dodgerblue3') +
    geom_text(
      aes(label = total_cluster_seqs),
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

make_significant_cluster_plot <- function(fisher_res, seq_df, level, alpha, fill_var){
  # for each significant cluster, shows the number of subjects or sequences (level)
  # in the cluster colored by fill_var

  md_sig <- fisher_res %>%
    dplyr::filter(p_value <= alpha) %>%
    dplyr::select(cluster_id) %>%
    dplyr::inner_join(seq_df, by = 'cluster_id') %>%
    dplyr::select(all_of(c(level, 'cluster_id', fill_var))) %>%
    distinct()

  md_sig[[fill_var]] <- as.factor(md_sig[[fill_var]])

  md_sig %>%
    ggplot(aes(x=cluster_id, fill=!!sym(fill_var))) +
    geom_bar(stat="count",
             width=0.85) +
    labs(x="Convergent Clone ID") +
    theme_bw() +
    scale_fill_brewer(palette = "Dark2") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    geom_text(aes(label = after_stat(count)),
              stat = "count",
              position = position_stack(vjust = 0.5),
              color="gray16")
}

make_fisher_overview_plot <- function(fisher_res, seq_df, condition, alpha, max_x=6){

  # get seqs per cluster
  seq_count_df <- seq_df %>%
    dplyr::group_by(cluster_id) %>%
    summarise(seq_count = n(), .groups = "drop")

  df_plot <- fisher_res %>%
    dplyr::filter(p_value <= alpha) %>%
    dplyr::left_join(seq_count_df, by='cluster_id') %>%
    dplyr::mutate(log2_odds_ratio = log2(odds_ratio))

  # assign a value to the infinite or clusters
  df_plot$log2_odds_ratio[is.infinite(df_plot$log2_odds_ratio)] <- max_x

  df_plot %>%
    ggplot(aes(x=log2_odds_ratio, y=in_cluster_in_condition, color=p_value)) +
    geom_point(aes(size=seq_count), stroke=1, alpha = 0.6) +
    scale_color_gradient(low = "red4", high = "white") +
    geom_label_repel(label=df_plot$cluster_id, size = 2, nudge_y = 0.4, nudge_x = 0.2, color="gray6") +
    geom_vline(xintercept = max_x-1, linetype = "dashed") +
    labs(x=paste0(condition, " odds ratio (log2)"),
         y=paste0("Number of ", condition, " subjects per cluster"),
         size = "# sequences per cluster",
         color = paste0("p-value"),
         title=paste0("Convergent clusters for ", condition, " group (p<",alpha, ")"))
}

##############################
### SET UP THE ENVIRONMENT ###
##############################

parser <- ArgumentParser(description = "Evaluate the output of a convergence tool.")

parser$add_argument('-c', '--counts', type = 'character', required = TRUE,
                    help = 'Comma-separated list of cluster x subject count matrices (one per group).')

parser$add_argument('-m', '--memberships', type = 'character', required = TRUE,
                    help = 'Comma-separated list of sequence-to-cluster membership tables (one per group).')

parser$add_argument('-i', '--tool_info', type = 'character', default = NULL,
                    help = 'Comma-separated list of tool info tables (one per group) containing run times.')

parser$add_argument('-mo', '--milo_objects', type = 'character', default = NULL,
                    help = 'Milo only: comma-separated list of Milo object RDS files (one per group), named <group_id>_milo.rds.')

parser$add_argument('-li', '--library_sizes', type = 'character', required = TRUE,
                    help = 'Library sizes table with subject_id, da_variable and depth columns for every subject in the dataset.')

parser$add_argument('-t', '--tool_type', type = 'character', default = 'cluster',
                    help = 'Type of tool output: "cluster" (each sequence in at most one cluster) or "milo".')

parser$add_argument('-n', '--tool_name', type = 'character', default = 'Tool',
                    help = 'Tool name used to label results in run_stats.tsv (i.e. "CDR3 Similarity").')

parser$add_argument('-o', '--output_dir', type = 'character', default = '.',
                    help = 'Output directory.')

parser$add_argument('-da', '--da_variable', type = 'character', default = 'status',
                    help = 'Stratification variable that should be used to determine differential abundance. There should be two levels in this factor/categorical variable.')

parser$add_argument('-dg', '--disease_group', type = 'character', default = 'disease',
                    help = 'Level of the DA variable to test for enrichment.')

parser$add_argument('-cg', '--control_group', type = 'character', default = 'control',
                    help = 'Level of the DA variable used as the reference for the beta-binomial test.')

parser$add_argument('-a', '--auc_var', type = 'character', default = 'FALSE',
                    help = 'Specify which column should be used for generating AUC curve (i.e. "simulated" or "binder"). Column type should be logical. If no AUC variable, set to FALSE.')

parser$add_argument('-al', '--alpha', type = 'double', default = 0.05,
                    help = 'Adjusted p-value threshold used to calculate the FDR.')

args <- parser$parse_args()

split_files <- function(x) if (is.null(x) || x == '') character() else strsplit(x, ',')[[1]]

COUNT_LOCS <- split_files(args$counts)
MEMBERSHIP_LOCS <- split_files(args$memberships)
TOOL_INFO_LOCS <- split_files(args$tool_info)
MILO_LOCS <- split_files(args$milo_objects)
LIB_SIZES_LOC <- args$library_sizes
TOOL_TYPE <- args$tool_type
TOOL_NAME <- args$tool_name
OUTPUT_DIR <- args$output_dir
DA_VAR <- args$da_variable
DISEASE_GP <- args$disease_group
CONTROL_GP <- args$control_group
ALPHA <- args$alpha

# for reticulate:
use_python("/opt/conda/envs/evaluate/bin/python", required = TRUE)

AUC_VAR <- args$auc_var
if (toupper(AUC_VAR) == 'FALSE'){
  AUC_VAR <- FALSE
  message('AUC will not be calculated.')
} else{
  message(paste0('AUC variable ', AUC_VAR, ' will be used.'))
}

if (!TOOL_TYPE %in% c('cluster', 'milo')){
  stop('Tool type must be either "cluster" or "milo".')
}

if (TOOL_TYPE == 'milo'){
  suppressPackageStartupMessages({
    library(miloR)
    library(SingleCellExperiment)
  })

  if (length(MILO_LOCS) == 0){
    stop('Milo objects are required to evaluate Milo results.')
  }
}

for (d in c(OUTPUT_DIR, file.path(OUTPUT_DIR, 'figures'), file.path(OUTPUT_DIR, 'tables'))){
  if (!dir.exists(d)){
    dir.create(d)
  }
}

########################
### READ IN THE DATA ###
########################

message(paste0('Loading library sizes: ', LIB_SIZES_LOC))
lib_sizes <- read.delim(LIB_SIZES_LOC, check.names = F, colClasses = c(subject_id = 'character')) %>%
  dplyr::arrange(subject_id)
row.names(lib_sizes) <- lib_sizes$subject_id

if (!DISEASE_GP %in% lib_sizes[[DA_VAR]]){
  stop(paste0('Disease group ', DISEASE_GP, ' not found in ', DA_VAR, ' column.'))
}

if (!CONTROL_GP %in% lib_sizes[[DA_VAR]]){
  stop(paste0('Control group ', CONTROL_GP, ' not found in ', DA_VAR, ' column.'))
}

print('Library sizes to be used:')
print(lib_sizes)

message('Loading cluster counts...')
count_mat <- read_count_matrix(COUNT_LOCS, lib_sizes$subject_id)
message(paste0(nrow(count_mat), ' clusters loaded.'))

message('Loading cluster memberships...')
memberships <- lapply(MEMBERSHIP_LOCS, function(f){
  read.delim(f, check.names = F,
             colClasses = c(sequence_id = 'character', subject_id = 'character',
                            group_id = 'character', cluster_id = 'character'))
}) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(membership_order = dplyr::row_number())

if (AUC_VAR != FALSE){
  if (!AUC_VAR %in% colnames(memberships)){
    stop(paste0('AUC variable ', AUC_VAR, ' not found in membership tables.'))
  }
  memberships[[AUC_VAR]] <- as.logical(memberships[[AUC_VAR]])
}

missing_clusters <- setdiff(na.omit(memberships$cluster_id), row.names(count_mat))
if (length(missing_clusters) > 0){
  stop(paste0(length(missing_clusters), ' clusters in membership tables are missing from count matrices.'))
}

# map each cluster to its group
cluster_groups <- memberships %>%
  dplyr::filter(!is.na(cluster_id)) %>%
  dplyr::distinct(cluster_id, group_id)

if (length(TOOL_INFO_LOCS) > 0){
  tool_info <- lapply(TOOL_INFO_LOCS, read.delim) %>% dplyr::bind_rows()
  time_taken <- sum(tool_info$time_min)
} else{
  time_taken <- NA
}

#############
### TESTS ###
#############

# tests to run and the label used in run_stats.tsv
if (TOOL_TYPE == 'milo'){
  tests <- c('GLM' = TOOL_NAME,
             'Fisher' = paste0(TOOL_NAME, ' + Fisher'),
             'Wilcoxon' = paste0(TOOL_NAME, ' + One-sided Wilcoxon'),
             'BetaBinomial' = paste0(TOOL_NAME, ' + Beta Binomial'),
             'OverdispersedZ' = paste0(TOOL_NAME, ' + Overdispersed Z-Test'),
             'PluginZ' = paste0(TOOL_NAME, ' + Plug-in Z-Test'))
} else{
  tests <- c('Fisher' = paste0(TOOL_NAME, ' + Fisher'),
             'Wilcoxon' = paste0(TOOL_NAME, ' + Wilcoxon'),
             'BetaBinomial' = paste0(TOOL_NAME, ' + Beta Binomial'),
             'OverdispersedZ' = paste0(TOOL_NAME, ' + Overdispersed Z-Test'),
             'PluginZ' = paste0(TOOL_NAME, ' + Plug-in Z-Test'))
}

test_results <- list()

message(paste0('Running Fisher Exact tests for ', DISEASE_GP, ' enrichment...'))
test_results[['Fisher']] <- run_fisher(count_mat, lib_sizes, DA_VAR, DISEASE_GP)

message(paste0('Running one-sided Wilcoxon tests for ', DISEASE_GP, ' enrichment...'))
cluster_freqs <- get_cluster_freqs(count_mat, lib_sizes)
test_results[['Wilcoxon']] <- run_wilcox(cluster_freqs, lib_sizes, DA_VAR, DISEASE_GP)

message(paste0('Running Beta Binomial tests for ', DISEASE_GP, ' enrichment...'))
test_results[['BetaBinomial']] <- run_beta_binomial(count_mat, lib_sizes, DA_VAR, DISEASE_GP, CONTROL_GP)

message(paste0('Running Overdispersed z-tests for ', DISEASE_GP, ' enrichment...'))
test_results[['OverdispersedZ']] <- run_overdispersed_z(count_mat, lib_sizes, DA_VAR, DISEASE_GP, CONTROL_GP)

message(paste0('Running Plug-in z-tests for ', DISEASE_GP, ' enrichment...'))
test_results[['PluginZ']] <- run_pluginz(count_mat, lib_sizes, DA_VAR, DISEASE_GP, CONTROL_GP)

write.table(data.frame(cluster_id = row.names(cluster_freqs), cluster_freqs, check.names = F),
            file.path(OUTPUT_DIR, 'tables', 'cluster_subject_freqs.tsv'),
            sep = '\t', row.names = F, quote = F)

if (TOOL_TYPE == 'milo'){

  # the GLM and spatial FDR must be computed separately for each group's graph
  glm_res_list <- list()
  spatial_fdr_list <- list()
  milo_objs <- list()

  for (f in MILO_LOCS){
    group <- sub('_milo\\.rds$', '', basename(f))
    message(paste0('Running Milo GLM for group ', group, '...'))

    milo <- readRDS(f)
    group_clusters <- paste0(group, '_', unlist(nhoodIndex(milo)))

    if (!all(group_clusters %in% row.names(count_mat))){
      stop(paste0('Neighborhoods in Milo object ', f, ' do not match the cluster count matrix.'))
    }

    glm_res_list[[group]] <- run_milo_glm(milo, group_clusters, lib_sizes, DA_VAR)

    # spatial FDR for the tests on this group's graph
    fisher_p <- test_results[['Fisher']]$p_value[match(group_clusters, test_results[['Fisher']]$cluster_id)]
    wilcox_p <- test_results[['Wilcoxon']]$p_value[match(group_clusters, test_results[['Wilcoxon']]$cluster_id)]
    betabinomial_p <- test_results[['BetaBinomial']]$p_value[match(group_clusters, test_results[['BetaBinomial']]$cluster_id)]
    overdispersedz_p <- test_results[['OverdispersedZ']]$p_value[match(group_clusters, test_results[['OverdispersedZ']]$cluster_id)]
    pluginz_p <- test_results[['PluginZ']]$p_value[match(group_clusters, test_results[['PluginZ']]$cluster_id)]

    spatial_fdr_list[[group]] <- data.frame(cluster_id = group_clusters,
                                            fisher_spatial_fdr = get_spatial_fdr(milo, fisher_p),
                                            wilcox_spatial_fdr = get_spatial_fdr(milo, wilcox_p),
                                            betabinomial_spatial_fdr = get_spatial_fdr(milo, betabinomial_p),
                                            overdispersedz_spatial_fdr = get_spatial_fdr(milo, overdispersedz_p),
                                            pluginz_spatial_fdr = get_spatial_fdr(milo, pluginz_p))

    milo_objs[[group]] <- milo
  }

  test_results[['GLM']] <- do.call(rbind, glm_res_list)

  spatial_fdr <- do.call(rbind, spatial_fdr_list)

  test_results[['Fisher']] <- test_results[['Fisher']] %>%
    dplyr::left_join(spatial_fdr %>% dplyr::select(cluster_id, spatial_fdr = fisher_spatial_fdr), by = 'cluster_id')

  test_results[['Wilcoxon']] <- test_results[['Wilcoxon']] %>%
    dplyr::left_join(spatial_fdr %>% dplyr::select(cluster_id, spatial_fdr = wilcox_spatial_fdr), by = 'cluster_id')

  test_results[['BetaBinomial']] <- test_results[['BetaBinomial']] %>%
    dplyr::left_join(spatial_fdr %>% dplyr::select(cluster_id, spatial_fdr = betabinomial_spatial_fdr), by = 'cluster_id')
  
  test_results[['OverdispersedZ']] <- test_results[['OverdispersedZ']] %>%
    dplyr::left_join(spatial_fdr %>% dplyr::select(cluster_id, spatial_fdr = overdispersedz_spatial_fdr), by = 'cluster_id')

  test_results[['PluginZ']] <- test_results[['PluginZ']] %>%
    dplyr::left_join(spatial_fdr %>% dplyr::select(cluster_id, spatial_fdr = pluginz_spatial_fdr), by = 'cluster_id')
}

# BH correction across ALL clusters, including across groups
for (test in names(test_results)){
  test_results[[test]]$p_adj <- p.adjust(test_results[[test]]$p_value, method = 'BH')
}

all_test_results <- lapply(names(test_results), function(test){
  test_results[[test]] %>%
    dplyr::mutate(test = test, .after = cluster_id)
}) %>%
  dplyr::bind_rows() %>%
  dplyr::left_join(cluster_groups, by = 'cluster_id') %>%
  dplyr::relocate(group_id, .after = cluster_id)

write.table(all_test_results,
            file.path(OUTPUT_DIR, 'tables', 'cluster_test_results.tsv'),
            sep = '\t', row.names = F, quote = F)

for (test in names(test_results)){
  test_results[[test]] %>%
    ggplot(aes(x = p_value)) +
    geom_histogram(color = 'white', binwidth = 0.01) +
    theme_bw() +
    labs(title = paste0(tests[[test]], ' P-Value Distribution')) +
    coord_cartesian(xlim = c(0, 1))

  ggsave(file.path(OUTPUT_DIR, 'figures', paste0('pvalue_hist_', test, '.png')),
         device = 'png', width = 8, height = 6, units = 'in')
}

#######################
### SEQUENCE SCORES ###
#######################

# Milo sequences are assigned to their neighborhood with the lowest spatial FDR.
# For other tools, each sequence is in one cluster and gets that cluster's p-value
score_col <- ifelse(TOOL_TYPE == 'milo', 'spatial_fdr', 'p_value')

seq_scores <- lapply(names(test_results), function(test){
  score_sequences(memberships, test_results[[test]], score_col)
})
names(seq_scores) <- names(test_results)

seq_scores_wide <- seq_scores[[1]] %>%
  dplyr::select(-c(assigned_cluster_id, score, p_adj))

for (test in names(seq_scores)){
  test_cols <- seq_scores[[test]] %>%
    dplyr::select(sequence_id, assigned_cluster_id, score, p_adj)
  colnames(test_cols)[-1] <- paste0(test, '_', colnames(test_cols)[-1])

  seq_scores_wide <- seq_scores_wide %>%
    dplyr::left_join(test_cols, by = 'sequence_id')
}

write.table(seq_scores_wide,
            file.path(OUTPUT_DIR, 'tables', 'sequence_scores.tsv'),
            sep = '\t', row.names = F, quote = F)

##################
### EVALUATION ###
##################

stat_table <- data.frame('tool' = unname(tests),
                         'total_seqs' = dplyr::n_distinct(memberships$sequence_id),
                         'total_subj' = nrow(lib_sizes),
                         'time (min)' = time_taken,
                         'subjects' = paste(lib_sizes$subject_id, collapse = ', '),
                         'depths' = paste(lib_sizes$depth, collapse = ', '),
                         check.names = F)

if (AUC_VAR != FALSE){
  message('Making AUC curves...')

  for (i in seq_along(tests)){
    test <- names(tests)[i]

    curve <- evaluate_curve(seq_scores[[test]]$score, seq_scores[[test]][[AUC_VAR]])

    write.table(curve$table,
                file.path(OUTPUT_DIR, 'tables', paste0('evaluation_curve_vals_', test, '.tsv')),
                sep = '\t', row.names = F, quote = F)

    plot_curves(curve, test)

    stat_table[i, 'AUROC'] <- curve$AUROC
    stat_table[i, 'AUPRC'] <- curve$AUPRC
    stat_table[i, 'average_precision'] <- curve$average_precision
    stat_table[i, 'FDR'] <- calc_FDR(seq_scores[[test]], AUC_VAR, ALPHA)
  }

  # document "purity" of clusters with hit sequences
  purity_stats <- get_purity_stats(memberships, AUC_VAR)

  write.table(purity_stats,
              file.path(OUTPUT_DIR, 'tables', 'purity_stats.tsv'),
              sep = '\t', row.names = F, quote = F)

  make_purity_plot(purity_stats, AUC_VAR)

  seq_info <- memberships %>% dplyr::distinct(sequence_id, .keep_all = TRUE)

  stat_table$tot_hits <- sum(seq_info[[AUC_VAR]], na.rm = T)
  stat_table$pct_hits <- mean(seq_info[[AUC_VAR]], na.rm = T) * 100
  stat_table$num_hit_clusters <- nrow(purity_stats)
  stat_table$avg_pct_hits <- mean(purity_stats$pct_hits)

  stat_table <- stat_table[c('tool', 'total_seqs', 'total_subj', 'tot_hits', 'pct_hits',
                             'num_hit_clusters', 'avg_pct_hits',
                             'AUROC', 'AUPRC', 'average_precision', 'FDR',
                             'time (min)', 'subjects', 'depths')]
}

write.table(stat_table,
            file.path(OUTPUT_DIR, 'tables', 'run_stats.tsv'),
            sep = '\t', row.names = F, quote = F)

################
### DATA VIZ ###
################

if (TOOL_TYPE == 'cluster'){

  seq_df <- memberships %>%
    dplyr::filter(!is.na(cluster_id)) %>%
    dplyr::left_join(lib_sizes[c('subject_id', DA_VAR)], by = 'subject_id')

  if (AUC_VAR != FALSE){
    make_significant_cluster_plot(test_results[['Fisher']], seq_df, 'sequence_id', 0.1, AUC_VAR)

    ggsave(file.path(OUTPUT_DIR, 'figures', 'hits_by_seq_id.png'),
           device="png", width=5, height=4, units="in")
  }

  make_significant_cluster_plot(test_results[['Fisher']], seq_df, 'subject_id', 0.1, DA_VAR)

  ggsave(file.path(OUTPUT_DIR, 'figures', paste0(DA_VAR, '_results_by_subj_id.png')),
         device="png", width=5, height=4, units="in")

  make_significant_cluster_plot(test_results[['Fisher']], seq_df, 'sequence_id', 0.1, DA_VAR)

  ggsave(file.path(OUTPUT_DIR, 'figures', paste0(DA_VAR, '_results_by_seq_id.png')),
         device="png", width=5, height=4, units="in")

  make_fisher_overview_plot(test_results[['Fisher']], seq_df, DISEASE_GP, 0.1, max_x = 6)

  ggsave(file.path(OUTPUT_DIR, 'figures', 'fisher_overview_disease.png'),
         device="png", width=8, height=8, units="in")

} else{

  # Milo neighborhood graphs colored by GLM logFC
  for (group in names(milo_objs)){
    milo <- buildNhoodGraph(milo_objs[[group]])

    group_res <- test_results[['GLM']] %>%
      dplyr::filter(cluster_id %in% paste0(group, '_', unlist(nhoodIndex(milo)))) %>%
      dplyr::transmute(Nhood = dplyr::row_number(), logFC, PValue = p_value, SpatialFDR = spatial_fdr)

    plotNhoodGraphDA(milo, group_res, alpha=0.05)

    ggsave(file.path(OUTPUT_DIR, 'figures', paste0(group, '_final_UMAP.png')),
           device = 'png',  width = 12, height = 6, units = 'in')
  }
}

message(paste0('Ending run: ', Sys.time()))

sessionInfo()
