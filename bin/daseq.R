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

    results_df <- data.frame(cluster_type = condition,
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

    results_df[[clone_id_col]] <- clone_id
    
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
  
  fisher_results_all_cond[[clone_id_col]] <- as.character(fisher_results_all_cond[[clone_id_col]])
  
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

make_auc_curve <- function(seq_table, p_val_col, cluster_id_col, auc_variable, name){
  # seq_table = table containing sequence-level information including auc_variable and cluster_id_col
  #             used to identify the correct category for each sequence. Also assumes p-values are already
  #             linked to sequences in this table
  # p_val_col = the column you want to use for p values to make AUC thresholds
  # cluster_id_col = IDs specifying clusters that were tested
  # auc_variable = variable to use for getting positives in the AUC curve
  # name = a name to specify for saving figures and tables (i.e. the name of the test)
  total_cells <- nrow(seq_table)
  unassigned_cells <- sum(is.na(seq_table[[p_val_col]]))
  assigned_cells <- total_cells - unassigned_cells
  seq_table[is.na(seq_table[[p_val_col]]), p_val_col] <- 1
  
  p_auc_thresholds <- sort(unique(seq_table[[p_val_col]]))
  # p_auc_thresholds <- quantile(X.cells$wilcox.adj.BH, seq(0, 1, 0.01), names=F)
  # auc_thresholds[1] <- auc_thresholds[1] - 1e-8
  
  # add to the largest to make sure the entire curve is captured
  tot_p_thresh <- length(p_auc_thresholds)
  p_auc_thresholds[tot_p_thresh] <- p_auc_thresholds[tot_p_thresh] + 1e-3
  
  p_auc_data <- lapply(p_auc_thresholds, function(thresh){
    
    DA_cells <- seq_table %>%
      dplyr::filter(!!sym(p_val_col) < thresh)
    
    non_DA_cells <- seq_table %>%
      dplyr::filter(!!sym(p_val_col) >= thresh)
    
    true_pos <- sum(DA_cells[[auc_variable]] == TRUE)
    false_neg <- sum(non_DA_cells[[auc_variable]] == TRUE)
    true_neg <- sum(non_DA_cells[[auc_variable]] == FALSE)
    false_pos <- sum(DA_cells[[auc_variable]] == FALSE)
    
    return(data.frame('TPR' = true_pos / (true_pos + false_neg),
                      'FPR' = 1 - (true_neg / (true_neg + false_pos))))
    
  })
  
  p_auc_df <- do.call(rbind, p_auc_data)
  p_auc_df$da_score_threshold <- p_auc_thresholds
  
  write.table(p_auc_df, 
              file.path(OUTPUT_DIR, 'tables', paste0('p_auc_curve_vals_', name, '.tsv')), 
              sep = '\t', row.names = F, quote = F)
  
  # get auroc
  p_auroc <- pracma::trapz(p_auc_df$FPR, p_auc_df$TPR)
  
  p_auc_df %>%
    ggplot(aes(x = FPR, y = TPR)) +
    geom_point() +
    geom_line() +
    labs(x = 'FPR',
         y = 'TPR',
         title = paste0(name, ' threshold ', round(min(p_auc_thresholds)), ' to ', round(max(p_auc_thresholds), 3)),
         subtitle = paste0(name, ' AUC: ', round(p_auroc, 3), '; ', 
                           prettyNum(assigned_cells, big.mark = ",", scientific = FALSE), '/', 
                           prettyNum(total_cells, big.mark = ",", scientific = FALSE), ' cells in DA clusters')) + 
    theme_minimal()
  
  ggsave(file.path(OUTPUT_DIR, 'figures', paste0('AUC_curve_', name, '.png')),
         device = 'png',
         width = 7,
         height = 6)

  return(p_auroc)

}

########################
### PREP ENVIRONMENT ###
########################

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

parser$add_argument('-dg', '--disease_group', type = 'character', default = 'disease',
                    help = 'The disease category.')

parser$add_argument('-m', '--k_min', type = 'integer', default = 20,
                    help = 'Minimum number of neighbors to use in KNN algorithm.')

parser$add_argument('-t', '--k_step', type = 'integer', default = 50,
                    help = 'Step value for KNN neighbor values.')

parser$add_argument('-x', '--k_max', type = 'integer', default = 450,
                    help = 'Maximum number of neighbors to use in KNN algorithm.')

parser$add_argument('-re', '--resolution', type = 'character', default = '0.01',
                    help = 'Clustering resolution. Use a larger value to get more regions. Note: will be cast as numeric in script.')

parser$add_argument('-mi', '--min_cell', type = 'character', default = 'NULL',
                    help = 'Minimum cells required to form a cluster. If NULL, it will be set to the minimum value in the k-vector.')

parser$add_argument('-a', '--auc_variable', type = 'character', default = FALSE,
                    help = 'Specify which column should be used for generating AUC curve (i.e. "simulated" or "binder"). Column type should be logical. If no AUC variable, set to FALSE.')

parser$add_argument('-v', '--vdj_info', type = 'logical', default = TRUE,
                    help = 'Is v call and j call information included in the metadata? Can apply to expression or embedding data.')

# TODO: can change to be more granular/option to plot at gene, family etc. level
# right now defaults to v_call and j_call columns and removes allele info
parser$add_argument('-sc', '--single_cell', type = 'logical', default = FALSE,
                    help = 'Input true if V(D)J info is present and contains paired heavy and light chain info.')

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

DISEASE_GP <- args$disease_group

KVEC <- seq(args$k_min, args$k_max, args$k_step)
RESOLUTION <- as.numeric(args$resolution)
MIN_CELL <- args$min_cell

if (MIN_CELL == 'NULL'){
    MIN_CELL <- NULL
} else{
  MIN_CELL <- as.numeric(MIN_CELL)
}

message(paste0('K nearest neighbor values: ', paste0(KVEC, collapse = ' ')))

VDJ <- args$vdj_info
SINGLE_CELL <- args$single_cell
AUC_VAR <- args$auc_variable

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
if (AUC_VAR != FALSE){
  umap_embeddings_coords <- umap_embeddings_coords %>%
    dplyr::left_join(md[c('id_col', AUC_VAR)], by = 'id_col')
  
  make_umap_viz(AUC_VAR, 'Hits', custom_pal = c('TRUE' = "red", 'FALSE' = "gray"))
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
    resolution = RESOLUTION,
    plot.embedding = X.embed,
    min.cell = MIN_CELL
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

if (AUC_VAR != FALSE){
  X.cells <- X.cells %>%
    dplyr::left_join(md[c(AUC_VAR, 'id_col')], by = 'id_col')
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
  ggplot(aes(x = pval.wilcoxon)) + 
  geom_histogram(color = 'white', binwidth = 0.01) + 
  theme_bw() +
  labs(title = 'DA-Seq Wilcoxon P-Value Distribution') +
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

#####################
### FISHER TEST ###
#####################

fisher_table <- get_combined_fisher_exact_table(hier_clone_df = X.cells %>% dplyr::filter(da.region.label != '0'),
                                                condition_set = c(DISEASE_GP),
                                                condition_col = DA_VAR,
                                                clone_id_col = 'da.region.label',
                                                count_col = 'subject_id',
                                                filter = FALSE)

write.table(fisher_table, file.path(OUTPUT_DIR, 'tables', 'fisher_table.tsv'), 
            sep="\t", quote = F, row.names = F)

fisher_sum <- fisher_table[c('da.region.label', 'p_value', 'fdr')]
colnames(fisher_sum) <- c('da.region.label', 'p_value_fisher', 'fdr_fisher')

X.cells <- dplyr::left_join(X.cells, fisher_sum, by = 'da.region.label')

# cell-level info
# restore ID column for writing
names(X.cells)[names(X.cells) == 'id_col'] <- ID_COL_NAME

write.table(X.cells, 
            file.path(OUTPUT_DIR, 'tables', 'da_seqs.tsv'), 
            sep='\t', quote = F, row.names = F)

names(X.cells)[names(X.cells) == ID_COL_NAME] <- 'id_col'

# make a summary of stats
stat_table <- data.frame('tool' = c('DAseq', 'DAseq + Fisher'),
                         'total_seqs' = nrow(X.cells),
                         'total_subj' = length(unique(X.cells$subject_id)),
                         'time (min)' = as.numeric(time_taken, units = "mins"),
                         'subjects' = paste(names(table(X.cells$subject_id)), collapse = ', '),
                         'depths' = paste(table(X.cells$subject_id), collapse = ', '),
                         check.names = F)

#######
# AUC #
#######

if (AUC_VAR != FALSE & AUC_VAR %in% colnames(X.cells)){ 
  message('Making AUC curves...')
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
    
    true_pos <- sum(DA_cells[[AUC_VAR]] == TRUE)
    false_neg <- sum(non_DA_cells[[AUC_VAR]] == TRUE)
    true_neg <- sum(non_DA_cells[[AUC_VAR]] == FALSE)
    false_pos <- sum(DA_cells[[AUC_VAR]] == FALSE)
    
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
  # do AUC curve with the Fisher Exact results
  fisher_auc <- make_auc_curve(X.cells, 'p_value_fisher', 'da.region.label', 'simulated', 'Fisher')

  # do AUC curve with the Wilcoxon results
  wilcox_auc <- make_auc_curve(X.cells, 'wilcox.adj.BH', 'da.region.label', 'simulated', 'Wilcoxon')
  
  stat_table[1, 'AUC'] <- c(wilcox_auc)
  stat_table[2, 'AUC'] <- c(fisher_auc)
  
  
  ###########
  # JACCARD #
  ###########
  jaccard_df <- X.cells
  
  jaccard_df <- jaccard_df %>%
    dplyr::mutate(p_under_0.005 = wilcox.adj.BH <= 0.005,
                  p_under_0.05 = wilcox.adj.BH <= 0.05,
                  p_under_0.1 = wilcox.adj.BH <= 0.1)
  
  # calc jaccard index
  jaccard_005 <- sum(jaccard_df[[AUC_VAR]] & jaccard_df$p_under_0.005, na.rm = T) / sum(jaccard_df[[AUC_VAR]] | jaccard_df$p_under_0.005, na.rm = T)
  jaccard_05 <- sum(jaccard_df[[AUC_VAR]] & jaccard_df$p_under_0.05, na.rm = T) / sum(jaccard_df[[AUC_VAR]] | jaccard_df$p_under_0.05, na.rm = T)
  jaccard_1 <- sum(jaccard_df[[AUC_VAR]] & jaccard_df$p_under_0.1, na.rm = T) / sum(jaccard_df[[AUC_VAR]] | jaccard_df$p_under_0.1, na.rm = T)
  
  # get Jaccard across a range
  # jaccard_thresholds <- quantile(jaccard_df$wilcox.adj.BH, seq(0, 1, 0.01), names=F, na.rm = T)
  jaccard_thresholds <- sort(unique(jaccard_df$wilcox.adj.BH))
  jaccard_thresholds <- jaccard_thresholds[!is.na(jaccard_thresholds)]
  
  # auc_thresholds[1] <- auc_thresholds[1] - 1e-8
  # jaccard_thresholds[101] <- jaccard_thresholds[101] + 1e-8
  
  jaccards <- sapply(jaccard_thresholds, function(thresh){
    j <- sum(jaccard_df[[AUC_VAR]] & jaccard_df$wilcox.adj.BH <= thresh, na.rm = T) / sum(jaccard_df[[AUC_VAR]] | jaccard_df$wilcox.adj.BH <= thresh, na.rm = T)
  })
  
  # get max Jaccard and its corresponding p-value
  Jaccard_max <- max(jaccards, na.rm = T)
  Jaccard_max_p <- jaccard_thresholds[which.max(jaccards)]
  
  jaccard_plot_df <- data.frame('Adjusted Wilcoxon P-Value Threshold' = jaccard_thresholds,
                                'Jaccard Similarity Index' = jaccards,
                                check.names = F)
  
  write.table(jaccard_plot_df, 
              file.path(OUTPUT_DIR, 'tables', 'jaccard_plot_vals.tsv'), 
              sep = '\t', row.names = F, quote = F)
  
  jaccard_plot <- jaccard_plot_df %>%
    ggplot(aes(x = !!sym('Adjusted Wilcoxon P-Value Threshold'), y = !!sym('Jaccard Similarity Index'))) +
    geom_point() +
    geom_line() +
    theme_bw() +
    labs(title = 'Jaccard Similarity Across Adjusted Wilcoxon P Thresholds',
         subtitle = paste0('Max Jaccard: ', round(Jaccard_max, 3), 
                           ' at adjusted P-value ', round(Jaccard_max_p, 3)))
  
  ggsave(filename = file.path(OUTPUT_DIR, 'figures', 'jaccard_plot.png'),
         plot = jaccard_plot,
         device = 'png',
         width = 7,
         height = 5)
  
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

if (AUC_VAR != FALSE & AUC_VAR %in% colnames(X.cells)){
  # also add what percentage of the cluster is simulated sequences
  hit_cts <- X.cells %>%
    dplyr::select('da.region.label', AUC_VAR) %>%
    dplyr::group_by(da.region.label) %>%
    dplyr::summarise(hit_seqs = sum(!!sym(AUC_VAR) == TRUE)) %>%
    dplyr::filter(da.region.label != 0)
  
  hit_cts$da.region.label <- as.character(hit_cts$da.region.label)
  
  subj_cts <- subj_cts %>%
    dplyr::left_join(hit_cts, by = 'da.region.label') %>%
    dplyr::mutate(pct_hits = hit_seqs / seqs_per_cluster)
}

# add the statistical info
subj_cts <- subj_cts %>%
  dplyr::left_join(DA.stat, by = 'da.region.label') %>%
  dplyr::left_join(fisher_sum, by = 'da.region.label')

for (test_choice in c(tests, 'p_value_fisher', 'fdr_fisher')){
  
  subj_cts[paste0(test_choice, '_significant')] <- subj_cts[test_choice] < 0.05
  
}

subj_cts <- subj_cts %>%
  dplyr::arrange(da.region.label, desc(pct_subj))

write.table(subj_cts, 
            file.path(OUTPUT_DIR, 'tables', 'cluster_cts.tsv'), 
            sep='\t', quote = F, row.names = F)

if (AUC_VAR != FALSE & AUC_VAR %in% colnames(X.cells)){
  
  purity_stats <- subj_cts %>%
    dplyr::filter(hit_seqs > 0) %>%
    dplyr::select(-c('subject_id', 'seqs_per_subj', 'pct_subj')) %>%
    distinct()

  # document "purity" of clusters with simulated sequences visually
  make_purity_plot(purity_stats, 'da.region.label', 'pct_hits', 'seqs_per_cluster', AUC_VAR)

  stat_table$num_hit_clusters <- nrow(purity_stats)
  stat_table$avg_pct_hits <- mean(purity_stats$pct_hits)
  stat_table$tot_hits <- c(sum(X.cells[[AUC_VAR]], na.rm = T)) 
  stat_table$pct_hits <- c(mean(X.cells[[AUC_VAR]], na.rm = T) * 100)
  stat_table$Jaccard_0.005 = c(jaccard_005, NA)
  stat_table$Jaccard_0.05 = c(jaccard_05, NA)
  stat_table$Jaccard_0.1 = c(jaccard_1, NA)
  stat_table$Jaccard_max = c(Jaccard_max, NA)
  stat_table$Jaccard_max_p = c(Jaccard_max_p, NA)
  
  stat_table <- stat_table[c('tool', 'total_seqs', 'total_subj', 'tot_hits', 'pct_hits',
                             'num_hit_clusters', 'avg_pct_hits', 'AUC', 'Jaccard_0.005', 'Jaccard_0.05',
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
