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

summarize_clusters <- function(fisher_table, df_hier_clones, clone_id_col, count_col, alpha){
  # get a table with info about the significant results coming from the fisher exact test table
  
  subj_info <- df_hier_clones %>%
    dplyr::group_by(!!sym(clone_id_col), !!sym(count_col)) %>%
    dplyr::summarise(count_per_cluster = n())
  
  sim_info <- df_hier_clones %>%
    dplyr::group_by(!!sym(clone_id_col)) %>%
    dplyr::summarise(simulated_per_cluster = sum(simulated == TRUE))
  
  cluster_cts <- df_hier_clones %>%
    dplyr::group_by(!!sym(clone_id_col)) %>%
    dplyr::summarise(total_cluster_seqs = n())
  
  all_df <- cluster_cts %>%
    dplyr::left_join(subj_info, by = clone_id_col) %>%
    dplyr::mutate(pct_per_cluster = count_per_cluster / total_cluster_seqs) %>%
    dplyr::left_join(sim_info, by = clone_id_col) %>%
    dplyr::mutate(pct_sim = simulated_per_cluster / total_cluster_seqs) %>%
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

##############################
### SET UP THE ENVIRONMENT ###
##############################

parser <- ArgumentParser(description = "Data location and Mal-ID algorithm hyperparameters.")

parser$add_argument('-md', '--metadata_loc', type = 'character', default = 'metadata',
                    help = 'File path for the metadata location.')

parser$add_argument('-o', '--output_dir', type = 'character', default = 'MalID_internal_output',
                    help = 'Specify an output directory location.')

parser$add_argument('-da', '--da_variable', type = 'character', default = 'status',
                    help = 'Stratification variable that should be used to determine for differential abundance. There should be two levels in this factor/categorical variable.')

parser$add_argument('-dg', '--disease_group', type = 'character', default = 'disease',
                    help = 'The disease category.')

parser$add_argument('-t', '--cluster_threshold', type = 'double', default = '0.15',
                    help = 'The distance threshold for forming clusters with hierarchical clones.')

parser$add_argument('-l', '--linkage_method', type = 'character', default = 'single',
                    help = 'The linkage method to be used in forming clusters.')

parser$add_argument('-v', '--vdj_info', type = 'logical', default = TRUE,
                    help = 'Is v call and j call information included in the metadata? Can apply to expression or embedding data.')

parser$add_argument('-sc', '--single_cell', type = 'logical', default = FALSE,
                    help = 'Input true if V(D)J info is present and contains paired heavy and light chain info.')

parser$add_argument('-r', '--remove_dups', type = 'logical', default = FALSE,
                    help = 'Will remove duplicate embeddings within an individual if TRUE.')

parser$add_argument('-si', '--simulated', type = 'logical', default = FALSE,
                    help = 'Specify whether input data is simulated or real.')

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
SIMULATED <- args$simulated
REMOVE_DUPS <- args$remove_dups

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

if (!'v_gene' %in% colnames(md)){
  # assume if v_gene not included, J probably isn't either
  md$v_gene <- alakazam::getGene(md$v_call, strip_d = F, omit_nl = F)
  md$v_allele <- alakazam::getAllele(md$v_call, strip_d = F, omit_nl = F)
  md$j_gene <- alakazam::getGene(md$j_call, strip_d = F, omit_nl = F)
  
}

if (SIMULATED){
  # make sure simulated is recognized
  md$simulated <- as.logical(md$simulated)
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
  labs(title = 'Mal-ID P-Value Distribution') +
  coord_cartesian(xlim = c(0, 1))

ggsave(file.path(OUTPUT_DIR, 'figures', 'pvalue_hist.png'),
       device = 'png', width = 8, height = 6, units = 'in')
  
# get summary info for clusters
summary <- summarize_clusters(fisher_table, 
                              convergent_clones, 
                             'convergent_clone_id', 'subject_id', 0.1)
  
write.table(summary, file.path(OUTPUT_DIR, 'tables', 'fisher_summary.tsv'), 
            sep="\t", quote = F, row.names = F)
  
# make plots
if (SIMULATED){
  
  make_significant_cluster_plot(fisher_table, convergent_clones, 
                                'id_col', 0.1, 'convergent_clone_id', 'simulated')
  
  ggsave(file.path(OUTPUT_DIR, 'figures', 'simulated_results_by_seq_id.png'), 
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

if (SIMULATED){
  cols_of_interest <- c(cols_of_interest, 'simulated')
}


sum1 <- convergent_clones[c(cols_of_interest)]
sum2 <- fisher_table[c('convergent_clone_id', 'p_value', 'odds_ratio', 'fdr')] %>% distinct()
sum <- dplyr::left_join(sum1, sum2, by = 'convergent_clone_id')

write.table(sum, file.path(OUTPUT_DIR, 'tables', "seq_summary.tsv"), 
            sep="\t", quote = F, row.names = F)

# get ending time after getting clusters & Fisher Test and making basic figures/tables
end_time <- Sys.time()
time_taken <- end_time - start_time

###################
### RUN SUMMARY ###
###################

# make a summary of stats
stat_table <- data.frame('tool' = c('Mal-ID Model 2'),
                         'total_seqs' = nrow(convergent_clones),
                         'total_subj' = length(unique(convergent_clones$subject_id)),
                         'time (min)' = as.numeric(time_taken, units = "mins"),
                         'subjects' = paste(names(table(convergent_clones$subject_id)), collapse = ', '),
                         'depths' = paste(table(convergent_clones$subject_id), collapse = ', '),
                         check.names = F)

#######
# AUC #
#######

if (SIMULATED){
    
  auc_thresholds <- sort(unique(fisher_table$p_value))
  # auc_thresholds <- quantile(fisher_table$p_value, seq(0, 1, 0.01), names=F)
  # auc_thresholds[1] <- auc_thresholds[1] - 1e-8
  
  # add to the largest to make sure the entire curve is captured
  tot_thresh <- length(auc_thresholds)
  auc_thresholds[tot_thresh] <- auc_thresholds[tot_thresh] + 1e-3
  
  auc_data <- lapply(auc_thresholds, function(thresh){
    
    # get whether the cells are DA or not at the given threshold
    
    # get significant clusters
    sig_clusters <- fisher_table %>%
      dplyr::filter(p_value < thresh) %>%
      dplyr::pull(convergent_clone_id) %>%
      unique()
    
    da_result <- convergent_clones[c('simulated', 'convergent_clone_id')] %>%
      dplyr::mutate(DA_cell = ifelse(convergent_clone_id %in% sig_clusters, TRUE, FALSE))
    
    da_result$simulated <- as.logical(da_result$simulated)
    
    true_pos <- sum(da_result$simulated == T & da_result$DA_cell == T)
    
    false_neg <- sum(da_result$simulated == T & da_result$DA_cell == F)
    
    true_neg <- sum(da_result$simulated == F & da_result$DA_cell == F)
    
    false_pos <- sum(da_result$simulated == F & da_result$DA_cell == T)
    
    TPR <- true_pos / (true_pos + false_neg)
    FPR <- 1 - (true_neg / (true_neg + false_pos))
    
    return(data.frame('TPR' = TPR,
                      'FPR' = FPR))
    
  })
  
  auc_df <- do.call(rbind, auc_data)
  auc_df$p_value <- auc_thresholds
  
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
         subtitle = paste0('AUC: ', round(auroc, 3))) + 
    theme_minimal()
  
  ggsave(file.path(OUTPUT_DIR, 'figures', 'AUC_curve.png'),
         device = 'png',
         width = 7,
         height = 6)
  
  stat_table$AUC <- c(auroc)
  
  ###########
  # JACCARD #
  ###########
  
  jaccard_df <- sum %>%
    dplyr::mutate(p_under_0.005 = p_value <= 0.005,
                  p_under_0.05 = p_value <= 0.05,
                  p_under_0.1 = p_value <= 0.1)
  
  # calc jaccard index
  jaccard_005 <- sum(jaccard_df$simulated & jaccard_df$p_under_0.005, na.rm = T) / sum(jaccard_df$simulated | jaccard_df$p_under_0.005, na.rm = T)
  jaccard_05 <- sum(jaccard_df$simulated & jaccard_df$p_under_0.05, na.rm = T) / sum(jaccard_df$simulated | jaccard_df$p_under_0.05, na.rm = T)
  jaccard_1 <- sum(jaccard_df$simulated & jaccard_df$p_under_0.1, na.rm = T) / sum(jaccard_df$simulated | jaccard_df$p_under_0.1, na.rm = T)
  
  # jaccard_thresholds <- seq(0, 1, 0.005)
  jaccard_thresholds <- sort(unique(jaccard_df$p_value))
  jaccard_thresholds <- jaccard_thresholds[!is.na(jaccard_thresholds)]
  
  # get Jaccard across a range
  jaccards <- sapply(jaccard_thresholds, function(thresh){
    j <- sum(jaccard_df$simulated & jaccard_df$p_value <= thresh, na.rm = T) / sum(jaccard_df$simulated | jaccard_df$p_value <= thresh, na.rm = T)
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
    
  stat_table$pct_simulated <- c(mean(jaccard_df$simulated, na.rm = T) * 100)
  stat_table$Jaccard_0.005 = jaccard_005
  stat_table$Jaccard_0.05 = jaccard_05
  stat_table$Jaccard_0.1 = jaccard_1
  stat_table$Jaccard_max = Jaccard_max
  stat_table$Jaccard_max_p = Jaccard_max_p
  
  stat_table <- stat_table[c('tool', 'total_seqs', 'total_subj', 'pct_simulated',
                             'AUC', 'Jaccard_0.005', 'Jaccard_0.05',
                             'Jaccard_0.1', 'Jaccard_max', 'Jaccard_max_p',
                             'time (min)', 'subjects', 'depths')]
  
}

write.table(stat_table, 
            file.path(OUTPUT_DIR, 'tables', 'run_stats.tsv'), 
            sep = '\t', row.names = F, quote = F)

message(paste0('Ending run: ', Sys.time()))

sessionInfo()
