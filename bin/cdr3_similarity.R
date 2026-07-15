#!/usr/bin/env Rscript
message(paste0('Starting run: ', Sys.time()))

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(forcats)
  library(RColorBrewer)
  library(patchwork)
  library(airr)
  library(alakazam)
  library(shazam)
  library(scoper)
  library(data.table)
  library(ggrepel)
  library(pracma)
  library(pbapply)
  library(argparse)
})

########################
### HELPER FUNCTIONS ###
########################

# write a function to perform the Fisher Exact Test for a specific cluster based on SUBJECTS in cluster
fisher_test_cluster <- function(df, subj_summary, input_convergent_clone_id, condition, condition_col = 'status', clone_id_col = 'convergent_clone_id', count_col = 'subject_id'){
  # df: input dataframe only containing clusters relevant for Fisher test set
  # (i.e. with 2 or more subjects in condition of interest)
  # subj_summary: summary of all subjects and statuses made BEFORE filtering to df
  # input_convergent_clone_id: clone we are performing Fisher test on
  # condition: condition we are testing for enrichment
  # count_col: the count column, i.e. which column are we getting our counts from? Could be subject_id, sequence_id
  
  # first, get all the subjects or sequences in a cluster
  in_cluster <- df %>%
    dplyr::filter(!!sym(clone_id_col) == input_convergent_clone_id) %>%
    dplyr::pull(count_col) %>%
    unique()
  
  # establish healthy and diseased groups in the entire test_smaller1 group
  # using subj summary
  subj_cond <- subj_summary %>% 
    dplyr::filter(!!sym(condition_col) == condition) %>% 
    dplyr::pull(count_col) %>%
    unique()
  
  tot_cond <- length(subj_cond)
  
  subj_not_cond <- subj_summary %>% 
    dplyr::filter(!!sym(condition_col) != condition) %>% 
    dplyr::pull(count_col) %>%
    unique()
  
  tot_not_cond <- length(subj_not_cond)
  
  # count those in cluster with condition
  in_cluster_cond <- length(intersect(in_cluster, subj_cond))
  
  # count those in cluster without condition
  in_cluster_not_cond <- length(intersect(in_cluster, subj_not_cond))
  
  # count those NOT in cluster with condition
  not_in_cluster_cond <- tot_cond - in_cluster_cond
  
  # count those NOT in cluster without condition
  not_in_cluster_not_cond <- tot_not_cond - in_cluster_not_cond
  
  # do not do the test if only one
  # should not be the case anyway because we pre-filtered
  # NOTE: changed to if 0 here - not possible, but if it happens there has been some mistake
  if (length(in_cluster) < 1){
    return(list(fisher_test_result = NA, 
                subjects_in_cluster = length(in_cluster),
                in_cluster_in_condition = in_cluster_cond,
                in_cluster_not_in_condition = in_cluster_not_cond))
  } else{
    
    # build contingency table to test for a CONDITON cluster
    #
    #                  cluster
    #                No    Yes
    #               ___________
    #            No|     |     |
    # condition    |_____|_____|
    #           Yes|     |     |
    #              |_____|_____|
    
    contingency_table <- matrix(c(not_in_cluster_not_cond, not_in_cluster_cond, in_cluster_not_cond, in_cluster_cond), 2, 2)
    
    # do fisher test
    return(list(fisher_test_result = fisher.test(contingency_table, alternative="greater"), 
                num_in_cluster = length(in_cluster),
                in_cluster_in_condition = in_cluster_cond,
                in_cluster_not_in_condition = in_cluster_not_cond,
                not_in_cluster_cond = not_in_cluster_cond,
                not_in_cluster_not_cond = not_in_cluster_not_cond,
                tot_cond = tot_cond,
                tot_not_cond = tot_not_cond))
  }
  
}

get_fisher_exact_table <- function(hier_clone_df, condition, condition_col = 'status', clone_id_col = 'convergent_clone_id', count_col = 'subject_id', filter = TRUE){
  # go from a hierarchical clones output dataframe
  # then get the clones worth doing fisher's exact on
  # do the fisher's exact test on every clone to test for healthy or diseased patients
  # depending on condition
  # count col establishes whether fisher testing is done at subject or sequence level
  
  if (filter){
    cat(paste0("Getting clones with at least 2 unique ", count_col, " in ", condition, " group..."), end="\n")
    # get the convergent clones with at least 2 subjects in the disease and/or 
    # 2 subjects in the healthy group
    convergent_clones_testable <- filter_hier_clones(hier_clone_df, condition, condition_col, clone_id_col, count_col)
    cat(paste0(length(convergent_clones_testable), " clones found passing filtering conditions for ", condition, " group."), end="\n")
  } else{
    
    convergent_clones_testable <- unique(hier_clone_df[[clone_id_col]])
    cat(paste0(length(convergent_clones_testable), " clones will be tested for ", condition, " group."), end="\n")
    
  }
  
  
  cat("Preparing data for Fisher's Exact test...", end="\n")
  
  # get the total subject information summarized BEFORE filtering
  # in case subjects will get lost
  subj_summary <- hier_clone_df %>%
    dplyr::select(!!sym(count_col), !!sym(condition_col)) %>%
    distinct() 
  
  # reduce the table to prepare for fisher and do tests faster
  hier_clone_df_fisher <- hier_clone_df %>%
    dplyr::filter(!!sym(clone_id_col) %in% convergent_clones_testable)
  
  cat("Completing Fisher's Exact tests...", end="\n")
  # do all the fisher tests
  fisher_results_all <- pbapply::pblapply(convergent_clones_testable, function(clone_id){
    
    # do test
    fisher_results <- fisher_test_cluster(hier_clone_df_fisher, subj_summary, clone_id, condition, condition_col, clone_id_col, count_col)

    results_df <- data.frame(convergent_clone_id = clone_id,
                             cluster_type = condition,
                             count_column = count_col,
                             p_value = NA,
                             odds_ratio = NA,
                             num_in_cluster = fisher_results[['num_in_cluster']],
                             in_cluster_in_condition = fisher_results[['in_cluster_in_condition']],
                             in_cluster_not_in_condition = fisher_results[['in_cluster_not_in_condition']],
                             not_in_cluster_in_condition = fisher_results[['not_in_cluster_cond']],
                             not_in_cluster_not_in_condition = fisher_results[['not_in_cluster_not_cond']],
                             total_in_condition = fisher_results[['tot_cond']],
                             total_not_in_condition = fisher_results[['tot_not_cond']])
    
    # check for NA (not enough info) but should be filtered out
    
    # if (fisher_results$num_in_cluster > 1){
      
    # pull out the fisher test results looking for a disease and a healthy cluster
    fisher <- fisher_results$fisher_test_result
    
    results_df$p_value <- fisher$p.value
    
    results_df$odds_ratio <- fisher$estimate
      
    # }
    
    return(results_df)
    
  })
  
  fisher_results_all <- do.call(rbind, fisher_results_all)
  fisher_results_all$fdr <- p.adjust(fisher_results_all$p_value, method="fdr")
  
  return(fisher_results_all)
}

get_combined_fisher_exact_table <- function(hier_clone_df, condition_set, condition_col = 'status', clone_id_col = 'convergent_clone_id', count_col = 'subject_id', filter = TRUE){
  # hier_clone_df: hierarchical clones df
  # condition set: character vector containing all conditions to be tested
  
  condition_fisher_dfs <- lapply(condition_set, function(condition){
    
    get_fisher_exact_table(hier_clone_df, condition, condition_col, clone_id_col, count_col, filter)
    
  })
  
  fisher_results_all_cond <- do.call(rbind, condition_fisher_dfs)
  
  fisher_results_all_cond$convergent_clone_id <- as.character(fisher_results_all_cond$convergent_clone_id)
  
  return(fisher_results_all_cond)
  
}

summarize_clusters <- function(fisher_table, df_hier_clones, clone_id_col, count_col, var_of_interest){
  # get a table with info about the significant results coming from the fisher exact test table
  
  subj_info <- df_hier_clones %>%
    dplyr::group_by(!!sym(clone_id_col), !!sym(count_col)) %>%
    dplyr::summarise(count_per_cluster = n())
  
  if (var_of_interest == F){
    hit_info <- df_hier_clones %>%
      dplyr::group_by(!!sym(clone_id_col)) %>%
      dplyr::summarise(hit_seqs = NA)
  } else{
    hit_info <- df_hier_clones %>%
      dplyr::group_by(!!sym(clone_id_col)) %>%
      dplyr::summarise(hit_seqs = sum(!!sym(var_of_interest) == TRUE))
  }

  cluster_cts <- df_hier_clones %>%
    dplyr::group_by(!!sym(clone_id_col)) %>%
    dplyr::summarise(total_cluster_seqs = n())
  
  all_df <- cluster_cts %>%
    dplyr::left_join(subj_info, by = clone_id_col) %>%
    dplyr::mutate(pct_per_cluster = count_per_cluster / total_cluster_seqs) %>%
    dplyr::left_join(hit_info, by = clone_id_col) %>%
    dplyr::mutate(pct_hits = hit_seqs / total_cluster_seqs) %>%
    dplyr::right_join(fisher_table, by = clone_id_col, relationship = "many-to-many")
  
  return(all_df)
  
}

make_significant_cluster_plot <- function(fisher_res, df_hier_clones, level, alpha, clone_id_col, fill_var){
  # for each type of cluster, shows the number of subjects from each study
  # in the cluster 
  
  # get sig clusters
  md_sig <- fisher_res %>%
    dplyr::filter(p_value <= alpha) %>%
    dplyr::left_join(df_hier_clones, by = clone_id_col)
  
  # adjust for level - sequence or subject IDs
  md_sig <- md_sig %>%
    dplyr::select(all_of(c(level, clone_id_col, fill_var, clone_id_col, 'cluster_type'))) %>%
    distinct()
  
  p <- md_sig %>%
    ggplot(aes(x=!!sym(clone_id_col), fill=!!sym(fill_var))) +
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
  
  if(n_distinct(md_sig$cluster_type) > 1){
    p + facet_wrap(vars(cluster_type), scales="free") 
  } else{
    p
  }
  
}

make_fisher_overview_plot <- function(fisher_table, df_hier_clones, level, condition, alpha, clone_id_col, max_x=6, current_fold = ''){
  # level is the level at which fisher tests were done - i.e. "subject" or "sequence"
  
  # get seqs per cluster
  seq_count_df <- df_hier_clones %>%
    dplyr::group_by(!!sym(clone_id_col)) %>%
    summarise(seq_count = n())
  
  # add seqs per cluster to fisher exact table
  df_plot <- fisher_table %>%
    dplyr::filter(cluster_type == condition) %>%
    dplyr::filter(p_value <= alpha) %>%
    dplyr::left_join(seq_count_df, by=clone_id_col) %>%
    dplyr::mutate(log2_odds_ratio = log2(odds_ratio))
  
  # assign a value to the infinite or clusters
  df_plot$log2_odds_ratio[is.infinite(df_plot$log2_odds_ratio)] <- max_x
  
  df_plot %>%
    ggplot(aes(x=log2_odds_ratio, y=in_cluster_in_condition, color=p_value)) +
    geom_point(aes(size=seq_count), stroke=1, alpha = 0.6) +
    scale_color_gradient(low = "red4", high = "white") +
    geom_label_repel(label=df_plot[[clone_id_col]], size = 2, nudge_y = 0.4, nudge_x = 0.2, color="gray6") +
    geom_vline(xintercept = max_x-1, linetype = "dashed") +
    labs(x=paste0(condition, " odds ratio (log2)"),
         y=paste0("Number of ", condition, " ", level, "s per cluster"),
         size = "# sequences per cluster",
         color = paste0("p-value"),
         title=paste0("Convergent clusters for ", condition, " group (p<",alpha, ") ", current_fold))
  
}

do_wilcox_test <- function(results_df, da_variable, disease_group, cluster_col){
  # perform Wilcoxon test on each cluster comparing counts in each group normalized
  # by subject depth.
  # results_df should contain information about the cluster ID (cluster_col), da_variable,
  # and subject IDs

  # Subject sequencing depths
  subject_depths <- results_df %>%
    dplyr::group_by(subject_id, !!sym(da_variable)) %>%
    dplyr::summarize(subj_depth = n())

  # Cluster frequencies per subject
  cluster_subject_freqs <- results_df %>%
    dplyr::group_by(!!sym(cluster_col), subject_id) %>%
    dplyr::summarize(cluster_sequences = n()) %>%
    dplyr::left_join(subject_depths, by = 'subject_id') %>%
    dplyr::mutate(normalized_freq = cluster_sequences / subj_depth) %>%
    tidyr::pivot_wider(id_cols = cluster_col, 
                      names_from = 'subject_id', 
                      values_from = 'normalized_freq', 
                      values_fill = 0) %>%
    as.data.frame(check.names = F)

  write.table(cluster_subject_freqs, 
              file.path(OUTPUT_DIR, 'tables', 'cluster_subject_freqs.tsv'), 
              sep = '\t', row.names = F, quote = F)

  row.names(cluster_subject_freqs) <- as.character(cluster_subject_freqs[[cluster_col]])
  cluster_subject_freqs <- cluster_subject_freqs %>%
                              dplyr::select(-!!sym(cluster_col))

  ctrl <- subject_depths %>%
            dplyr::filter(!!sym(da_variable) != disease_group) %>%
            dplyr::pull(subject_id)

  dis <- subject_depths %>%
            dplyr::filter(!!sym(da_variable) == disease_group) %>%
            dplyr::pull(subject_id)

  wilcox_res_list <- lapply(row.names(cluster_subject_freqs), function(clust){
    x <- as.numeric(cluster_subject_freqs[clust,ctrl])
    y <- as.numeric(cluster_subject_freqs[clust,dis])
    p_val <- wilcox.test(x, y, alternative = 'less', exact = FALSE)$p.value
    return(data.frame(convergent_clone_id = clust,
                      p_value = p_val))
  })

  wilcox_result <- do.call(rbind, wilcox_res_list)
  wilcox_result$p_adj <- p.adjust(wilcox_result$p_value, method = 'BH')

  write.table(wilcox_result, 
              file.path(OUTPUT_DIR, 'tables', 'wilcox_res.tsv'), 
              sep = '\t', row.names = F, quote = F)

  return(wilcox_result)          
}


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

make_auc_curve <- function(results_table, seq_table, p_val_col, cluster_id_col, auc_variable, name){
  # results table = table containing some kind of test result (Fisher Exact, Wilcox, etc.)
  # seq_table = table containing sequence-level information including auc_variable and cluster_id_col
  #             used to identify the correct category for each sequence
  # p_val_col = the column you want to use for p values to make AUC thresholds
  # cluster_id_col = IDs specifying clusters that were tested
  # auc_variable = variable to use for getting positives in the AUC curve
  # name = a name to specify for saving figures and tables (i.e. the name of the test)

  auc_thresholds <- sort(unique(results_table[[p_val_col]]))
  
  # add to the largest to make sure the entire curve is captured
  tot_thresh <- length(auc_thresholds)
  auc_thresholds[tot_thresh] <- auc_thresholds[tot_thresh] + 1e-3
  
  auc_data <- lapply(auc_thresholds, function(thresh){
    
    # get whether the cells are DA or not at the given threshold
    
    # get significant clusters
    sig_clusters <- results_table %>%
      dplyr::filter(!!sym(p_val_col) < thresh) %>%
      dplyr::pull(cluster_id_col) %>%
      unique()
    
    da_result <- seq_table[c(auc_variable, cluster_id_col)] %>%
      dplyr::mutate(DA_cell = ifelse(!!sym(cluster_id_col) %in% sig_clusters, TRUE, FALSE))
    
    da_result[[auc_variable]] <- as.logical(da_result[[auc_variable]])
    
    true_pos <- sum(da_result[[auc_variable]] == T & da_result$DA_cell == T)
    
    false_neg <- sum(da_result[[auc_variable]] == T & da_result$DA_cell == F)
    
    true_neg <- sum(da_result[[auc_variable]] == F & da_result$DA_cell == F)
    
    false_pos <- sum(da_result[[auc_variable]] == F & da_result$DA_cell == T)
    
    TPR <- true_pos / (true_pos + false_neg)
    FPR <- 1 - (true_neg / (true_neg + false_pos))
    
    return(data.frame('TPR' = TPR,
                      'FPR' = FPR))
    
  })
  
  auc_df <- do.call(rbind, auc_data)
  auc_df[[p_val_col]] <- auc_thresholds
  
  write.table(auc_df, 
              file.path(OUTPUT_DIR, 'tables', paste0('auc_curve_vals_', name, '.tsv')), 
              sep = '\t', row.names = F, quote = F)
  
  # get auroc
  auroc <- pracma::trapz(auc_df$FPR, auc_df$TPR)
  
  auc_df %>%
    ggplot(aes(x = FPR, y = TPR)) +
    geom_point() +
    geom_line() +
    labs(title = paste0('Alpha Threshold ', round(min(auc_thresholds)), ' to ', round(max(auc_thresholds), 3)),
         subtitle = paste0(name, ' AUC: ', round(auroc, 3))) + 
    theme_minimal()
  
  ggsave(file.path(OUTPUT_DIR, 'figures', paste0('AUC_curve_', name, '.png')),
         device = 'png',
         width = 7,
         height = 6)

  return(auroc)
}

##############################
### SET UP THE ENVIRONMENT ###
##############################

parser <- ArgumentParser(description = "Data location and Mal-ID algorithm hyperparameters.")

parser$add_argument('-md', '--metadata_loc', type = 'character', default = 'metadata',
                    help = 'File path for the metadata location.')

parser$add_argument('-o', '--output_dir', type = 'character', default = 'MalID_internal_output',
                    help = 'Specify an output directory location.')

parser$add_argument('-da', '--da_variable', type = 'character', default = 'status',
                    help = 'Stratification variable that should be used to determine differential abundance. There should be two levels in this factor/categorical variable.')

parser$add_argument('-dg', '--disease_group', type = 'character', default = 'disease',
                    help = 'The disease category.')

parser$add_argument('-t', '--cluster_threshold', type = 'double', default = '0.15',
                    help = 'The distance threshold for forming clusters with hierarchical clones.')

parser$add_argument('-l', '--linkage_method', type = 'character', default = 'single',
                    help = 'The linkage method to be used in forming clusters.')

parser$add_argument('-a', '--auc_var', type = 'character', default = FALSE,
                    help = 'Specify which column should be used for generating AUC curve (i.e. "simulated" or "binder"). Column type should be logical. If no AUC variable, set to FALSE.')

parser$add_argument('-v', '--vdj_info', type = 'logical', default = TRUE,
                    help = 'Is v call and j call information included in the metadata? Can apply to expression or embedding data.')

parser$add_argument('-sc', '--single_cell', type = 'logical', default = FALSE,
                    help = 'Input TRUE if V(D)J info is present and contains paired heavy and light chain info.')

parser$add_argument('-r', '--remove_dups', type = 'logical', default = FALSE,
                    help = 'Will remove duplicate embeddings within an individual if TRUE.')


# Parse the arguments
args <- parser$parse_args()

MD_LOC <- args$metadata_loc
OUTPUT_DIR <- args$output_dir

DA_VAR <- args$da_variable

DISEASE_GP <- args$disease_group

THRESH <- args$cluster_threshold
LINKAGE <- args$linkage_method

VDJ <- args$vdj_info
SINGLE_CELL <- args$single_cell
AUC_VAR <- args$auc_var
REMOVE_DUPS <- args$remove_dups

if (AUC_VAR != FALSE){
  message(paste0('AUC variable ', AUC_VAR, ' will be used.'))
} else{
  message('AUC will not be calculated.')
}

if (REMOVE_DUPS){
  message('Duplicate embeddings within a subject will be collapsed.')
}

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

########################
### READ IN THE DATA ###
########################

# metadata
message(paste0('Loading metadata: ', MD_LOC))

tryCatch(
  
  {
    md <- airr::read_rearrangement(MD_LOC)
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

# change to a generic id column - look for sequence id first
if ('sequence_id' %in% colnames(md)){
  
  ID_COL_NAME <- 'sequence_id'
  names(md)[names(md) == 'sequence_id'] <- 'id_col'
  
} else if ('cell_id' %in% colnames(md)){
  
  ID_COL_NAME <- 'cell_id'
  names(md)[names(md) == 'cell_id'] <- 'id_col'
  
} else {
  
  stop('No cell_id or sequence_id column found in airr data.')
  
}

# add gene and allele info
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

if (AUC_VAR != FALSE){
  # make sure simulated is recognized
  md[[AUC_VAR]] <- as.logical(md[[AUC_VAR]])
}

message(paste0(dplyr::n_distinct(md$subject_id), ' unique subjects and ',
               dplyr::n_distinct(md$sample_id), ' unique samples found.'))

# remove any NA junctions which will mess up our analysis
num_NA_junc <- sum(is.na(md$junction))
message(paste0('WARNING: Removing ', num_NA_junc, ' sequences with NA in junction column.'))
md <- md %>%
  dplyr::filter(!is.na(junction))

# remove dups if necessary
if (REMOVE_DUPS){
  old_seq_num <- nrow(md)
  
  md <- md %>%
    distinct(v_gene, j_gene, cdr3_aa, subject_id, .keep_all = TRUE)
  
  new_seq_num <- nrow(md)
  
  seqs_removed <- old_seq_num - new_seq_num
  message(paste0('Duplicates removed. ', seqs_removed, ' sequences removed. New total: ', new_seq_num))
}

# measure how long the Mal-ID process itself takes
start_time <- Sys.time()

# STRATEGY IF WE HAVE A HIGHER SAMPLE SIZE
###############################################################################

###############
### CLUSTER ###
###############

# need to make clone IDs
# for initial test, use all the sequences
if (SINGLE_CELL){
  convergent_clones <- scoper::hierarchicalClones(md,
                                                  threshold=THRESH,
                                                  method="aa",
                                                  linkage=LINKAGE,
                                                  normalize="len",
                                                  junction="junction",
                                                  v_call="v_call", 
                                                  j_call="j_call",
                                                  clone="convergent_clone_id",
                                                  fields=NULL,
                                                  cell_id="id_col",
                                                  locus="locus",
                                                  only_heavy=FALSE,
                                                  split_light=FALSE,
                                                  first=FALSE,
                                                  cdr3=FALSE, 
                                                  mod3=FALSE,
                                                  max_n=0, 
                                                  nproc=16,
                                                  verbose=T, log=NULL,
                                                  summarize_clones=FALSE)
  
} else{
  convergent_clones <- scoper::hierarchicalClones(md,
                                                  threshold=THRESH,
                                                  method="aa",
                                                  linkage=LINKAGE,
                                                  normalize="len",
                                                  junction="junction",
                                                  v_call="v_call", 
                                                  j_call="j_call",
                                                  clone="convergent_clone_id",
                                                  fields=NULL,
                                                  cell_id=NULL,
                                                  locus="locus",
                                                  only_heavy=TRUE,
                                                  split_light=FALSE,
                                                  first=FALSE,
                                                  cdr3=FALSE, 
                                                  mod3=FALSE,
                                                  max_n=0, 
                                                  nproc=16,
                                                  verbose=T, log=NULL,
                                                  summarize_clones=FALSE)
}


convergent_clones$clone_id_full <- paste(convergent_clones$convergent_clone_id, convergent_clones$subject_id, sep ='_')

#########################
### FISHER EXACT TEST ###
#########################
# Fisher's exact test

fisher_table <- get_combined_fisher_exact_table(hier_clone_df = convergent_clones,
                                                condition_set = c(DISEASE_GP),
                                                condition_col = DA_VAR,
                                                clone_id_col = 'convergent_clone_id',
                                                count_col = 'subject_id',
                                                filter = FALSE)
  

  
write.table(fisher_table, file.path(OUTPUT_DIR, 'tables', 'fisher_table.tsv'), 
            sep="\t", quote = F, row.names = F)

fisher_table %>%
  ggplot(aes(x = p_value)) + 
  geom_histogram(color = 'white', binwidth = 0.01) + 
  theme_bw() +
  labs(title = 'CDR3 Similarity Fisher P-Value Distribution') +
  coord_cartesian(xlim = c(0, 1))

ggsave(file.path(OUTPUT_DIR, 'figures', 'pvalue_hist_fisher.png'),
       device = 'png', width = 8, height = 6, units = 'in')
  
# get summary info for clusters
summary_fisher <- summarize_clusters(fisher_table, 
                                     convergent_clones, 
                                     'convergent_clone_id', 'subject_id', AUC_VAR)

summary_fisher <- summary_fisher %>%
  rename(
    p_value_fisher = p_value,
    fdr_fisher = fdr,
    odds_ratio_fisher = odds_ratio
)
  
# make plots
if (AUC_VAR != FALSE){
  
  make_significant_cluster_plot(fisher_table, convergent_clones, 
                                'id_col', 0.1, 'convergent_clone_id', AUC_VAR)
  
  ggsave(file.path(OUTPUT_DIR, 'figures', 'hits_by_seq_id.png'), 
         device="png", width=5, height=4, units="in")
}

make_significant_cluster_plot(fisher_table, convergent_clones, 
                              'subject_id', 0.1, 'convergent_clone_id', DA_VAR)

ggsave(file.path(OUTPUT_DIR, 'figures', paste0(DA_VAR, '_results_by_subj_id.png')), 
       device="png", width=5, height=4, units="in")

make_significant_cluster_plot(fisher_table, convergent_clones, 
                              'id_col', 0.1, 'convergent_clone_id', DA_VAR)

ggsave(file.path(OUTPUT_DIR, 'figures', paste0(DA_VAR, '_results_by_seq_id.png')), 
       device="png", width=5, height=4, units="in")

make_fisher_overview_plot(fisher_table, convergent_clones, 
                          'subject', DISEASE_GP, 0.1, 'convergent_clone_id', max_x = 6)

ggsave(file.path(OUTPUT_DIR, 'figures', 'fisher_overview_disease.png'), 
       device="png", width=8, height=8, units="in")


# AUC summary
cols_of_interest <- c('id_col', 'v_gene', 'j_gene', 'subject_id', 'convergent_clone_id')

if (AUC_VAR != FALSE){
  cols_of_interest <- c(cols_of_interest, AUC_VAR)
}


sum1 <- convergent_clones[c(cols_of_interest)]
sum2 <- fisher_table[c('convergent_clone_id', 'p_value', 'odds_ratio', 'fdr')] %>% distinct()
colnames(sum2) <- c('convergent_clone_id', 'p_value_fisher', 'odds_ratio_fisher', 'fdr_fisher')

sum <- dplyr::left_join(sum1, sum2, by = 'convergent_clone_id')

# write.table(sum, file.path(OUTPUT_DIR, 'tables', "seq_summary.tsv"), 
#             sep="\t", quote = F, row.names = F)

# get ending time after getting clusters & Fisher Test and making basic figures/tables
end_time <- Sys.time()
time_taken <- end_time - start_time

#####################
### WILCOXON TEST ###
#####################

# information required: 
# depth per person
# number of sequences from each person in each cluster
# group identity of each person
# use to get ratio of seqs in cluster / depth for each
# then do Wilcox test disease vs. control for each --> use apply loop

message('Completing Wilcoxon DA tests...')

wilcox_res <- do_wilcox_test(convergent_clones, DA_VAR, DISEASE_GP, 'convergent_clone_id')

wilcox_sum <- wilcox_res
colnames(wilcox_sum) <- c('convergent_clone_id', 'p_value_wilcox', 'fdr_wilcox')

sum <- dplyr::left_join(sum, wilcox_sum, by = 'convergent_clone_id')

write.table(sum, file.path(OUTPUT_DIR, 'tables', "seq_summary.tsv"), 
            sep="\t", quote = F, row.names = F)

wilcox_res %>%
  ggplot(aes(x = p_value)) + 
  geom_histogram(color = 'white', binwidth = 0.01) + 
  theme_bw() +
  labs(title = 'CDR3 Similarity Wilcoxon P-Value Distribution') +
  coord_cartesian(xlim = c(0, 1))

ggsave(file.path(OUTPUT_DIR, 'figures', 'pvalue_hist_wilcox.png'),
       device = 'png', width = 8, height = 6, units = 'in')

summary <- summary_fisher %>%
  dplyr::inner_join(wilcox_sum, by = 'convergent_clone_id')
  
write.table(summary, file.path(OUTPUT_DIR, 'tables', 'cluster_subj_summary.tsv'), 
            sep="\t", quote = F, row.names = F)
  

###################
### RUN SUMMARY ###
###################

# make a summary of stats
stat_table <- data.frame('tool' = c('CDR3 Similarity + Fisher', 'CDR3 Similarity + Wilcoxon'),
                         'total_seqs' = nrow(convergent_clones),
                         'total_subj' = length(unique(convergent_clones$subject_id)),
                         'time (min)' = as.numeric(time_taken, units = "mins"),
                         'subjects' = paste(names(table(convergent_clones$subject_id)), collapse = ', '),
                         'depths' = paste(table(convergent_clones$subject_id), collapse = ', '),
                         check.names = F)

#######
# AUC #
#######

if (AUC_VAR != FALSE){
  message('Making AUC curves...')
  # do AUC curve with the Fisher Exact results
  fisher_auc <- make_auc_curve(fisher_table, convergent_clones, 'p_value', 'convergent_clone_id', AUC_VAR, 'Fisher')

  # do AUC curve with the Wilcoxon results
  wilcox_auc <- make_auc_curve(wilcox_res, convergent_clones, 'p_value', 'convergent_clone_id', AUC_VAR, 'Wilcoxon')
  
  stat_table[1, 'AUC'] <- c(fisher_auc)
  stat_table[2, 'AUC'] <- c(wilcox_auc)
  
  ###########
  # JACCARD #
  ###########
  
  jaccard_df <- sum %>%
    dplyr::mutate(p_under_0.005 = p_value_fisher <= 0.005,
                  p_under_0.05 = p_value_fisher <= 0.05,
                  p_under_0.1 = p_value_fisher <= 0.1)
  
  # calc jaccard index
  jaccard_005 <- sum(jaccard_df[[AUC_VAR]] & jaccard_df$p_under_0.005, na.rm = T) / sum(jaccard_df[[AUC_VAR]] | jaccard_df$p_under_0.005, na.rm = T)
  jaccard_05 <- sum(jaccard_df[[AUC_VAR]] & jaccard_df$p_under_0.05, na.rm = T) / sum(jaccard_df[[AUC_VAR]] | jaccard_df$p_under_0.05, na.rm = T)
  jaccard_1 <- sum(jaccard_df[[AUC_VAR]] & jaccard_df$p_under_0.1, na.rm = T) / sum(jaccard_df[[AUC_VAR]] | jaccard_df$p_under_0.1, na.rm = T)
  
  # jaccard_thresholds <- seq(0, 1, 0.005)
  jaccard_thresholds <- sort(unique(jaccard_df$p_value_fisher))
  jaccard_thresholds <- jaccard_thresholds[!is.na(jaccard_thresholds)]
  
  # get Jaccard across a range
  jaccards <- sapply(jaccard_thresholds, function(thresh){
    j <- sum(jaccard_df[[AUC_VAR]] & jaccard_df$p_value_fisher <= thresh, na.rm = T) / sum(jaccard_df[[AUC_VAR]] | jaccard_df$p_value_fisher <= thresh, na.rm = T)
  })
  
  # get max Jaccard and its corresponding p-value
  Jaccard_max <- max(jaccards)
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

  purity_stats <- summary %>%
    dplyr::filter(hit_seqs > 0) %>%
    dplyr::select(-c('subject_id', 'count_per_cluster', 'pct_per_cluster', 'count_column')) %>%
    distinct()

  # document "purity" of clusters with simulated sequences visually
  make_purity_plot(purity_stats, 'convergent_clone_id', 'pct_hits', 'total_cluster_seqs', AUC_VAR)

  stat_table$num_hit_clusters <- nrow(purity_stats)
  stat_table$avg_pct_hits <- mean(purity_stats$pct_hits)
  stat_table$tot_hits <- c(sum(jaccard_df[[AUC_VAR]], na.rm = T)) 
  stat_table$pct_hits <- c(mean(jaccard_df[[AUC_VAR]], na.rm = T) * 100)
  stat_table$Jaccard_0.005 = c(jaccard_005, NA) # not doing the Jaccard stats for Wilcoxon, at least for now
  stat_table$Jaccard_0.05 = c(jaccard_05, NA)
  stat_table$Jaccard_0.1 = c(jaccard_1, NA)
  stat_table$Jaccard_max = c(Jaccard_max, NA)
  stat_table$Jaccard_max_p = c(Jaccard_max_p, NA)
  
  stat_table <- stat_table[c('tool', 'total_seqs', 'total_subj', 'tot_hits', 'pct_hits',
                             'num_hit_clusters', 'avg_pct_hits', 
                             'AUC', 'Jaccard_0.005', 'Jaccard_0.05',
                             'Jaccard_0.1', 'Jaccard_max', 'Jaccard_max_p',
                             'time (min)', 'subjects', 'depths')]
  
}

write.table(stat_table, 
            file.path(OUTPUT_DIR, 'tables', 'run_stats.tsv'), 
            sep = '\t', row.names = F, quote = F)

message(paste0('Ending run: ', Sys.time()))

sessionInfo()
