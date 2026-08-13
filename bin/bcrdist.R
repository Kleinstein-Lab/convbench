#!/usr/bin/env Rscript

message(paste0('Starting run: ', Sys.time()))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(airr)
  library(alakazam)
  library(ggrepel)
  library(argparse)
  library(reticulate)
  library(arrow)
  library(pracma)
  library(Biostrings)
  library(scales)
})

####################################
######### HELPER FUNCTIONS #########
####################################

### Translate nt into aa ###
# !!! HOW TO DEAL WITH N!!! 
# !!! HOW TO DEAL WITH INCOMPLETE CODONS !!!
# Anything not falling into the codons become empty character
# The reading frame is fixed from the 1st
# If gaps are not explicitly labelled, it will cause frameshift
# Which is how real sequence looks like
translate_wGaps = function(nt_vec){
  
  # Build dictionary (standard code + gap mapping)
  code = c(Biostrings::GENETIC_CODE) #, "..." = ".")
  
  # Vectorize over the character strings
  aa_vec = vapply(nt_vec, function(seq) {
    
    # Anything not ACTG becomes an empty character
    seq = gsub("[^ACTG]", "", seq)
    
    len = nchar(seq)

    # If sequence is empty
    if(len < 3){
        aa = ""
        return(aa)
    }
    
    # Transform deletions into gaps
    # seq = gsub('-','.',seq)
    
    # Split string into triplets
    codons = substring(seq, seq(1, len - 2, by = 3), seq(3, len, by = 3))
    
    # Look up in dictionary
    # unrecognized stuff returns NA
    aa = code[codons]
    
    # Convert NAs to IMGT gaps
    # aa[is.na(aa)] = "."
    # Convert NAs into empty strings
    # Mainly deal with non-multiple-of-3 length
    aa[is.na(aa)] = "" 
    
    # Paste back into a single character string
    paste0(aa, collapse = "")
  }, FUN.VALUE = character(1), USE.NAMES = FALSE)
  
  return(aa_vec)
}

### Plot the distribution of n_subj/seq of all clusters ###
plot_cluster_composition = function(cluster_df,metric=NULL){

    # cluster_df has to have n_sequences
    cluster_df = cluster_df %>% filter(n_sequences > 1)

    log10p1 = trans_new(
        name = "log10p1",
        transform = function(x) log10(x + 1),
        inverse = function(x) 10^x - 1
    )
    
    if(metric == 'n_subjects'){
        p = cluster_df %>% ggplot(aes(x = .data[[metric]])) +
            geom_histogram(color='white', binwidth=1) + 
            labs(x='# subjects per cluster', y='# clusters')
    } else if(metric == 'n_sequences'){
        p = cluster_df %>% ggplot(aes(x = .data[[metric]])) +
            geom_histogram(binwidth=1) +
            labs(x='# sequences per cluster', y='# clusters')
    } else if(metric == 'n_unique_sequences'){
        p = cluster_df %>% ggplot(aes(x = .data[[metric]])) +
            geom_histogram(binwidth=1) +
            labs(x='# unique sequences per cluster', y='# clusters')
    } else{
        message(paste0("No ",metric," found in the clone composition summary!"))
        return(NULL)
    }

    p = p + 
        theme_bw() +
        scale_y_continuous(
            trans = log10p1,
            breaks = c(0, 9, 99, 999), 
            labels = c("0", "9", "99", "999"))
    return(p)
}

### Copied from Mal-ID helper functions ###
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

get_fisher_exact_table <- function(hier_clone_df, condition, condition_col = 'status', clone_id_col = 'convergent_clone_id', count_col = 'subject_id', subj_summary = NULL, filter = TRUE){
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

  # in the case of ASCs, total subj information can also be input as a parameter
  if (is.null(subj_summary)){
    subj_summary <- hier_clone_df %>%
      dplyr::select(!!sym(count_col), !!sym(condition_col)) %>%
      distinct() 
  }
  
  # reduce the table to prepare for fisher and do tests faster
  hier_clone_df_fisher <- hier_clone_df %>%
    dplyr::filter(!!sym(clone_id_col) %in% convergent_clones_testable)
  
  cat("Completing Fisher's Exact tests...", end="\n")
  # do all the fisher tests
  fisher_results_all <- lapply(convergent_clones_testable, function(clone_id){
    
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

get_combined_fisher_exact_table <- function(hier_clone_df, condition_set, condition_col = 'status', clone_id_col = 'convergent_clone_id', count_col = 'subject_id', subj_summary = NULL, filter = TRUE){
  # hier_clone_df: hierarchical clones df
  # condition set: character vector containing all conditions to be tested
  
  condition_fisher_dfs <- lapply(condition_set, function(condition){
    
    get_fisher_exact_table(hier_clone_df, condition, condition_col, clone_id_col, count_col, subj_summary, filter)
    
  })
  
  fisher_results_all_cond <- do.call(rbind, condition_fisher_dfs)
  
  fisher_results_all_cond$convergent_clone_id <- as.character(fisher_results_all_cond$convergent_clone_id)
  
  return(fisher_results_all_cond)
  
}

summarize_clusters <- function(fisher_table, df_hier_clones, clone_id_col, count_col, alpha, var_of_interest){
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

make_fisher_overview_plot <- function(fisher_table, df_hier_clones, level, condition, current_fold, alpha, clone_id_col, max_x=6){
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
         title=paste0("Convergent clusters for ", condition, " group (p<",alpha, "), ", current_fold))
  
}

do_wilcox_test <- function(results_df, da_variable, disease_group, cluster_col, subject_depths = NULL){
  # perform Wilcoxon test on each cluster comparing counts in each group normalized
  # by subject depth.
  # results_df should contain information about the cluster ID (cluster_col), da_variable,
  # and subject IDs

  # Subject sequencing depths
  if(is.null(subject_depths)){
    subject_depths <- results_df %>%
      dplyr::group_by(subject_id, !!sym(da_variable)) %>%
      dplyr::summarize(depth = n())
  }

  # Cluster frequencies per subject
  cluster_subject_freqs <- results_df %>%
    dplyr::group_by(!!sym(cluster_col), subject_id) %>%
    dplyr::summarize(cluster_sequences = n()) %>%
    dplyr::left_join(subject_depths, by = 'subject_id') %>%
    dplyr::mutate(normalized_freq = cluster_sequences / depth) %>%
    tidyr::pivot_wider(id_cols = cluster_col, 
                      names_from = 'subject_id', 
                      values_from = 'normalized_freq', 
                      values_fill = 0) %>%
    as.data.frame(check.names = F)

  row.names(cluster_subject_freqs) <- as.character(cluster_subject_freqs[[cluster_col]])
  cluster_subject_freqs <- cluster_subject_freqs %>%
                              dplyr::select(-!!sym(cluster_col))

  # add individuals not in the ASC if needed
  absent_subj <- setdiff(subject_depths$subject_id, colnames(cluster_subject_freqs))
    
  if (length(absent_subj) > 0){
    message(paste0('Adding missing subjects ', paste(absent_subj, collapse = ', '), ' to cluster subject frequency table.'))
    
    # add to counts
    cluster_subject_freqs[, absent_subj] <- 0

  } else{
    message('No missing subjects detected.')
  }
  
  write.table(cluster_subject_freqs,
              file.path(OUTPUT_DIR, 'tables', 'cluster_subject_freqs.tsv'), 
              sep = '\t', row.names = T, quote = F)

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

parser <- ArgumentParser(description = "Data location and BCRdist hyperparameters.")

parser$add_argument('-he', '--helper_loc', type = 'character', default = 'bin/tcrdist3.py',
                    help = 'File path for the helper function location.')

parser$add_argument('-md', '--metadata_loc', type = 'character', default = 'metadata',
                    help = 'File path for the metadata location.')

parser$add_argument('-li', '--library_sizes', type = 'character', default = NULL,
                    help = 'File path for the library sizes file location.')

parser$add_argument('-o', '--output_dir', type = 'character', default = 'tcrdist3_fisher_output',
                    help = 'Specify an output directory location.')

parser$add_argument('-da', '--da_variable', type = 'character', default = 'status',
                    help = 'Stratification variable that should be used to determine for differential abundance. There should be two levels in this factor/categorical variable.')

parser$add_argument('-dg', '--disease_group', type = 'character', default = 'disease',
                    help = 'The disease category.')

parser$add_argument('-t', '--cluster_threshold', type = 'double', default = '60',
                    help = 'The distance threshold for forming clusters with hierarchical clones.')

# Not enabled currently
# parser$add_argument('-l', '--linkage_method', type = 'character', default = 'single',
#                    help = 'The linkage method to be used in forming clusters.')

# Not enabled currently
# parser$add_argument('-sc', '--single_cell', type = 'character', default = "FALSE",
#                    help = 'Input true if V(D)J info is present and contains paired heavy and light chain info.')

parser$add_argument('-a', '--auc_var', type = 'character', default = FALSE,
                    help = 'Specify which column should be used for generating AUC curve (i.e. "simulated" or "binder"). Column type should be logical. If no AUC variable, set to FALSE.')

# Added for running tcrdist
# e.g. "/home/zhaochenye/anaconda3/envs/env_scissor/bin/python"
parser$add_argument('-py', '--python_location', type = 'character', default = NULL,
                    help = 'Point reticulate to the correct python environment where the tcrdist3 is installed. Mutually exclusive with -e/--conda_env.')

# Added for running tcrdist
# e.g. "env_scissor"
parser$add_argument('-e', '--conda_env', type = 'character', default = NULL,
                    help = 'Point reticulate to the correct conda environment where the tcrdist3 is installed. Mutually exclusive with -py/--python_env.')

# Added for tcrdist (T) or bcrdist (F)
parser$add_argument('-i', '--infer_v2cdr', type = 'character', default = "TRUE",
                    help = 'Specify whether to infer CDR from V gene (T) or to directly input CDR sequences (F).')

# Deprecated
# Adding IMGT gaps to sequences
# e.g. '/home/zhaochenye/share/germlines/imgtdb_base_cvg/'
# parser$add_argument('-imgt', '--imgtdb', type = 'character', default = NULL,
#                    help = 'Only used when -i FALSE, specify where the IMGT database was installed to add gaps to CDR regions.')

# parser$add_argument('-og', '--organism', type = 'character', default = 'human',
#                    help = 'Only used when -i FALSE, specify the organism under IMGT database, should match the input data.')

# parser$add_argument('-lc', '--locus', type = 'character', default = 'Ig',
#                    help = 'Only used when -i FALSE, specify the locus under IMGT database, should match the input data.')

# Added for tcrdist3 parameters
parser$add_argument('-c', '--cpu', type = 'integer', default = 1,
                    help = 'Specify the number of cpus used to run tcrdist3')

# Added for tcrdist3 parameters
parser$add_argument('-cs', '--chunk_size', type = 'integer', default = 1000,
                    help = 'Specify the chunk size to break input data for running tcrdist3')

# Added for tcrdist3 parameters
parser$add_argument('-r', '--radius', type = 'double', default = 60,
                    help = 'Specify the radius (maximum distance) in the pairwise distance matrix.')

# Added for checking the pipeline
# parser$add_argument('-dm', '--save_distancemtx', type = 'character', default = NULL,
#                    help = 'If provided, the location where the pairwise distance matrix is saved (.pt). Note: This can be time-consuming.')

# Added for tcrdist3 parameters
parser$add_argument('-m', '--mode', type = 'character', default = 'greedy',
                    help = 'Specify the clustering method: "greedy" or "single".')


# Parse the arguments
args <- parser$parse_args()
### set the working directory ###
# TCRDIST_SCRIPT <- args$helper_loc # loaded after directed to python
TCRDIST_SCRIPT <- Sys.which('tcrdist3.py')

MD_LOC <- args$metadata_loc
LIB_SIZES_LOC <- args$library_sizes
MD_NAME <- stringr::str_split_i(basename(MD_LOC), '_md', 1)

OUTPUT_DIR <- args$output_dir

DA_VAR <- args$da_variable

DISEASE_GP <- args$disease_group

THRESH <- args$cluster_threshold

AUC_VAR <- args$auc_var

PYTHON_LOC <- args$python_location
CONDA_ENV <- args$conda_env
if(!is.null(PYTHON_LOC)){
  use_python(PYTHON_LOC)
} else if(!is.null(CONDA_ENV)){
  use_condaenv(CONDA_ENV)
} else{ 
  stop('Python or Conda environment should be provided where tcrdist3 is installed.')
}
source_python(TCRDIST_SCRIPT)

V2CDR <- as.logical(args$infer_v2cdr)

# IMGTDB <- args$imgtdb
# ORGANISM<- args$organism
# LOCUS <- args$locus
# if(!V2CDR & is.null(IMGTDB)){
#  stop('Cannot directly input CDR sequences without IMGT database location.')
# }

CPU <- args$cpu
CHUNK_SIZE <- args$chunk_size
RADIUS <- args$radius
MODE <- args$mode

# SAVE_DIST <- args$save_distancemtx
# if(!is.null(SAVE_DIST)){
#    SAVE_DIST_PATH = SAVE_DIST
#    SAVE_DIST = TRUE
# } else{
#    SAVE_DIST_PATH = ''
#    SAVE_DIST = FALSE
# }

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

#############################################
######## Preprocessings for tcrdist3 ########
#############################################

# md_ = md
# source('/home/zhaochenye/Projects/convergence/tcrdist3_cvg/igblast_funcs.R')
# 86533 - addGaps removes an extra C in the original junction column
# That C makes the junction_length = 58, not multiple of 3
# AddGaps can't deal with junction correction
# Or this C is a missing/un-labelled insertion
# 69877 - insertions at the beginning of junction

if(V2CDR){
  message('Infer CDR sequences from V gene!')
    
  # Add pmhc_aa column for code consistency
  md$pmhc_aa = ''
  if(!"cdr1_aa" %in% colnames(md)){
      md$cdr1_aa = ''
  }
  if(!"cdr2_aa" %in% colnames(md)){
      md$cdr2_aa = ''
  }
    
} else{
  message('Directly input CDR sequences!')
  
  ### Add IMGT gaps ###
  # imgt_dir = '/home/zhaochenye/share/germlines/imgtdb_base_cvg/'
  # md = addGaps(db = md,gapdb = IMGTDB,organism = ORGANISM,locus = LOCUS)
  
  # validation
  # cdr1_aa_ts = translate_wGaps(md$cdr1)
  # cdr1_aa_sa = substr(md$sequence_alignment_aa,start = 27,stop = 38)
  # cdr1_aa = md$cdr1_aa
  # cdr1 = data.frame(cdr1 = md$cdr1,cdr1_aa_ts,cdr1_aa_sa,cdr1_aa)
  # idx = which(substr(md$sequence_alignment_aa,start = 30,stop = 30) == 'X')
  # substr(md$sequence_alignment[idx],start = 79,stop = 114)

  ### CDR1/2/2.5 sequences ###
  # md$sequence_alignment_aa = translate_wGaps(md$sequence_alignment)
  # md$cdr1_aa = substr(md$sequence_alignment_aa,start = 27,stop = 38)
  # md$cdr2_aa = substr(md$sequence_alignment_aa,start = 56,stop = 65)
  # md$pmhc_aa = substr(md$sequence_alignment_aa,start = 81,stop = 86)

  ### Translate nt into aa ###
  if(!"cdr1_aa" %in% colnames(md)){
      if(!"cdr1" %in% colnames(md)){
          stop('Please have either cdr1_aa or cdr1 column in the input.')
      }
      md$cdr1_aa = translate_wGaps(md$cdr1)
  }
  if(!"cdr2_aa" %in% colnames(md)){
      if(!"cdr2" %in% colnames(md)){
          stop('Please have either cdr2_aa or cdr2 column in the input.')
      }
      md$cdr2_aa = translate_wGaps(md$cdr2)
  }
  # Add pmhc_aa column for code consistency
  md$pmhc_aa = ''

  ### Deletions spanning whole CDR1/2 ###
  na_cdr1 = which(is.na(md$cdr1_aa))
  if(length(na_cdr1)){
      md$cdr1_aa[na_cdr1] = ""
  }
  na_cdr2 = which(is.na(md$cdr2_aa))
  if(length(na_cdr2)){
      md$cdr2_aa[na_cdr2] = ""
  }
    
}

# standardize column names
colnames(md) <- tolower(colnames(md))

# create artificial sample_id copies from subject ID if not present
if (!'sample_id' %in% colnames(md)){
  md$sample_id <- md$subject_id
}

# change to a generic id column
if ('cell_id' %in% colnames(md)){
  
  ID_COL_NAME <- 'cell_id'
  names(md)[names(md) == 'cell_id'] <- 'id_col'
  
} else if ('sequence_id' %in% colnames(md)){
  
  ID_COL_NAME <- 'sequence_id'
  names(md)[names(md) == 'sequence_id'] <- 'id_col'
  
} else{
  
  stop('No cell_id or sequence_id column found in airr data.')
  
}

### Adding allele info to fit tcrdist input requirement ###
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
if(!'junction_aa' %in% colnames(md)){
    if(!'junction' %in% colnames(md)){
        stop('No junction or junction_aa found in airr data.')
    } else{
        message('No junction_aa found in airr data. Translate junction into aa sequences.')
        md$junction_aa = translate_wGaps(md$junction)
    }
}
num_NA_junc <- sum(is.na(md$junction) | is.na(md$junction_aa))
message(paste0('WARNING: Removing ', num_NA_junc, ' sequences with NA in junction or junction_aa column.'))
md <- md %>%
  dplyr::filter(!is.na(junction)) %>%
  dplyr::filter(!is.na(junction_aa)) 

### Rename the columns to fit tcrdist input format ###

# !!!!!! Maybe should group identical V_CDR3 pairs first !!!!!!
# !!!!!! count column is inactive now !!!!!!

seqs = select(.data = md,id_col,subject_id,v_allele,junction_aa,j_allele,cdr1_aa,cdr2_aa,pmhc_aa)
colnames(seqs) = c('id_col','subject_id','v_b_gene','cdr3_b_aa','j_b_gene','cdr1_b_aa','cdr2_b_aa','pmhc_b_aa')
seqs = seqs %>%
  mutate(seq_id = paste0('s',1:nrow(seqs)), # seq_id is necessary for linking correct sequence ID after tcrdist3
         count = 1)
# cat("Saving modified metadata for running tcrdist3...", end="\n")
# write.table(seqs, file.path(OUTPUT_DIR, 'intermediates', "tcrdist_input.tsv"), 
#            sep="\t", quote = F, row.names = F)

#############################################
############### Run tcrdist3 ################
#############################################

# measure how long the tcrdist3 process itself takes
start_time <- Sys.time()

seqs_new_py = tryCatch({
    tcrdist3_clusters(seqs,
                      cpus=CPU,
                      chunk_size=CHUNK_SIZE, # the minibatch size
                      radius=RADIUS, # introducing sparsity, any distance above radius will be 0
                      id_col='seq_id', # the order of seqs will change, a unique id is necessary for matching the seqs afterwards
                      mode = MODE, 
                      threshold=THRESH,
                      #save_dist_mtx=SAVE_DIST,save_path=paste0(SAVE_DIST_PATH,'/pwdist.pt'),
                      v2cdr=V2CDR) # whether to infer CDRs from V genes
    
    }, error = function(e) {
      message("PYTHON CRASHED! Here is the actual Python traceback:")
      print(reticulate::py_last_error())
      stop(e) # re-throw the error to halt execution
  })

seqs_new <- py_to_r(seqs_new_py)
seqs_new$subject_id = NULL # avoid duplication with md
convergent_clone = right_join(md,seqs_new,by = 'id_col')

convergent_clone$count <- unlist(convergent_clone$count)

print(convergent_clone[1:5,c('id_col', 'count')])

# Saving the clusters
# Only for testing, will comment out
cat(paste0("Saving clustering results (convergent_clone)..."), end="\n")
write.table(convergent_clone, file.path(OUTPUT_DIR, 'tables',"bcrdist_clusters.tsv"), 
              sep="\t", quote = F, row.names = F)

# Characterizing the clusters
clone_comp = convergent_clone %>% 
    group_by(convergent_clone_id) %>%
    summarise(
        n_subjects = n_distinct(subject_id),
        n_sequences = n(),
        #n_unique_sequences=n_distinct(v_gene, j_gene, junction_aa)
    )
write.table(clone_comp, file.path(OUTPUT_DIR, 'tables',"cluster_summary.tsv"), sep="\t", quote = F, row.names = F)

# Plot
clone_comp = clone_comp %>% filter(n_sequences > 1)
p1 = plot_cluster_composition(clone_comp,metric='n_subjects')
p2 = plot_cluster_composition(clone_comp,metric='n_sequences')
# p3 = plot_cluster_composition(clone_comp,metric='n_unique_sequences')

if(!is.null(p1)) ggsave(file.path(OUTPUT_DIR, 'figures','n_subjects_distribution.png'), p1, device="png", width=5, height=4, units="in")
if(!is.null(p2)) ggsave(file.path(OUTPUT_DIR, 'figures','n_sequences_distribution.png'), p2, device="png", width=5, height=4, units="in")
# if(!is.null(p3)) ggsave(file.path(OUTPUT_DIR, 'figures', paste0(current_fold, '_n_unique_sequences_distribution.png')), p3, device="png", width=5, height=4, units="in")
#}

#########################
### FISHER EXACT TEST ###
#########################
# Fisher's exact test

if(!is.null(LIB_SIZES_LOC)){
  lib_sizes <- read.csv(LIB_SIZES_LOC, sep = '\t') %>% as.data.frame()
} else{
  lib_sizes <- NULL
}

fisher_table <- get_combined_fisher_exact_table(hier_clone_df = convergent_clone,
                                                  condition_set = c(DISEASE_GP),
                                                  condition_col = DA_VAR,
                                                  clone_id_col = 'convergent_clone_id',
                                                  count_col = 'subject_id',
                                                  subj_summary = lib_sizes,
                                                  filter = FALSE)
  
write.table(fisher_table, file.path(OUTPUT_DIR, 'tables', "fisher_table.tsv"), sep="\t", quote = F, row.names = F)

# get summary info for clusters
summary_fisher <- summarize_clusters(fisher_table, convergent_clone, 'convergent_clone_id', 'subject_id', 0.1, AUC_VAR)

summary_fisher <- summary_fisher %>%
  dplyr::rename(
    p_value_fisher = p_value,
    fdr_fisher = fdr,
    odds_ratio_fisher = odds_ratio
)

# write.table(summary_fisher, file.path(OUTPUT_DIR, 'tables', "fisher_summary.tsv"), sep="\t", quote = F, row.names = F)
  
# make plots
if(AUC_VAR != FALSE){
    make_significant_cluster_plot(fisher_table, convergent_clone,'id_col', 0.1, 'convergent_clone_id', AUC_VAR)
    ggsave(file.path(OUTPUT_DIR, 'figures', 'results_by_seq_id.png'), device="png", width=5, height=4, units="in")
}

  
make_significant_cluster_plot(fisher_table, convergent_clone, 'subject_id', 0.1, 'convergent_clone_id', 'status')
ggsave(file.path(OUTPUT_DIR, 'figures', 'results_by_subj_id.png'), device="png", width=5, height=4, units="in")
  
make_fisher_overview_plot(fisher_table, convergent_clone,'subject', DISEASE_GP, "1fold", 0.1, 'convergent_clone_id', max_x = 6)

ggsave(file.path(OUTPUT_DIR, 'figures', 'fisher_overview_disease.png'), device="png", width=8, height=8, units="in")
  
# make_fisher_overview_plot(fisher_table, convergent_clones[[current_fold]], 'subject', 'control', current_fold, 0.1, 'convergent_clone_id', max_x = 6)
# ggsave(file.path(OUTPUT_DIR, 'figures', 'fisher_results', 'fisher_overview_control.png'), device="png", width=8, height=8, units="in")

# AUC summary
cols_of_interest <- c('id_col', 'v_gene', 'j_gene', 'subject_id', 'convergent_clone_id')

if (AUC_VAR != FALSE){
  cols_of_interest <- c(cols_of_interest, AUC_VAR)
}

sum1 <- convergent_clone[c(cols_of_interest)]
sum2 <- fisher_table[c('convergent_clone_id', 'p_value', 'odds_ratio', 'fdr')] %>% distinct()
colnames(sum2) <- c('convergent_clone_id', 'p_value_fisher', 'odds_ratio_fisher', 'fdr_fisher')

sum <- dplyr::left_join(sum1, sum2, by = 'convergent_clone_id')
  
# write.table(sum, file.path(OUTPUT_DIR, 'tables', paste0(MD_NAME, "_seq_summary.tsv")), sep="\t", quote = F, row.names = F)

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

wilcox_res <- do_wilcox_test(convergent_clone, DA_VAR, DISEASE_GP, 'convergent_clone_id', subject_depths = lib_sizes)

wilcox_sum <- wilcox_res
colnames(wilcox_sum) <- c('convergent_clone_id', 'p_value_wilcox', 'fdr_wilcox')

sum <- dplyr::left_join(sum, wilcox_sum, by = 'convergent_clone_id')

write.table(sum, file.path(OUTPUT_DIR, 'tables', paste0(MD_NAME, "_seq_summary.tsv")), 
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
stat_table <- data.frame('tool' = c('BCRdist + Fisher', 'BCRdist + Wilcoxon'),
                         'total_seqs' = nrow(sum),
                         'total_subj' = length(unique(sum$subject_id)),
                         'time (min)' = as.numeric(time_taken, units = "mins"),
                         'subjects' = paste(names(table(sum$subject_id)), collapse = ', '),
                         'depths' = paste(table(sum$subject_id), collapse = ', '),
                         check.names = F)

#######
# AUC #
#######

if (AUC_VAR != FALSE){

  message('Making AUC curves...')
  # do AUC curve with the Fisher Exact results
  fisher_auc <- make_auc_curve(fisher_table, convergent_clone, 'p_value', 'convergent_clone_id', AUC_VAR, 'Fisher')

  # do AUC curve with the Wilcoxon results
  wilcox_auc <- make_auc_curve(wilcox_res, convergent_clone, 'p_value', 'convergent_clone_id', AUC_VAR, 'Wilcoxon')
  
  stat_table[1, 'AUC'] <- c(fisher_auc)
  stat_table[2, 'AUC'] <- c(wilcox_auc)
    
  ###########
  # JACCARD #
  ###########
  
  # sum <- read.csv(file.path(OUTPUT_DIR, 'tables', paste0(MD_NAME, "_seq_summary.tsv")), sep = '\t')
  
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
  
  ggsave(file.path(OUTPUT_DIR, 'figures', 'jaccard_plot.png'),
          plot = jaccard_plot,
          device = 'png',
          width = 7,
          height = 5)
  
  ###################
  ### RUN SUMMARY ###
  ###################

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
