#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
    library(stringr)
    library(ggplot2)
    library(pracma)
    library(argparse)
})

##########################
### set up environment ###
##########################

parser <- ArgumentParser(description = "Data location and Mal-ID algorithm hyperparameters.")

parser$add_argument('-md', '--metadata_locs', type = 'character', default = 'metadata',
                    help = 'File paths for the metadata locations.')

parser$add_argument('-a', '--auc_var', type = 'character', default = FALSE,
                    help = 'Specify which column should be used for generating AUC curve (i.e. "simulated" or "binder"). Column type should be logical. If no AUC variable, set to FALSE.')

parser$add_argument('-t', '--tool', type = 'character', default = FALSE,
                    help = 'Specify which tool was used to generate the results tables.')


# Parse the arguments
args <- parser$parse_args()

MD_LOCS <- args$metadata_locs

AUC_VAR <- args$auc_var

if (AUC_VAR == 'FALSE'){
    AUC_VAR <- as.logical(AUC_VAR)
    message('AUC will not be calculated.')
}

TOOL <- args$tool

# get all the metadata files (1 per AUC)
MD_LOCS <- strsplit(MD_LOCS, ',')[[1]]

if (!dir.exists('tables')){
    dir.create('tables')
}

if (!dir.exists('figures')){
    dir.create('figures')
}

########################
### helper functions ###
########################

evaluation_curve <- function(res, p_col, binder_col, tool = '', simplify_p = F){
  
  # get AUPRC baseline - fraction of positive events
  auprc_baseline <- mean(res[[binder_col]], na.rm = T)

  invalid_seqs <- sum(is.na(res[[p_col]]))
  total_seqs <- nrow(res)
  valid_seqs <- total_seqs - invalid_seqs
  
  # change the sequences with no nhood to a min p of 1
  res[is.na(res[[p_col]]), p_col] <- 1

  # avoid floating point errors
  res[[p_col]] <- round(res[[p_col]], 6)
  
  if (simplify_p){
    auc_thresholds_init <- sort(unique(res[[p_col]]))
    
    auc_thresholds <- c(auc_thresholds_init[1])
    
    for (i in 2:length(auc_thresholds_init)){
      final_auc_thresholds <- length(auc_thresholds)
      if (auc_thresholds_init[i] - auc_thresholds[final_auc_thresholds] > 1e-5){ # has to be different enough from last
        auc_thresholds <- c(auc_thresholds, auc_thresholds_init[i])
      }
    }
  } else{
    auc_thresholds <- sort(unique(res[[p_col]]))
  }
  
  # add to the largest to make sure the entire curve is captured
  tot_thresh <- length(auc_thresholds)
  auc_thresholds[tot_thresh] <- auc_thresholds[tot_thresh] + 1e-3
  
  # get binder info (or whatever the AUC var is)
  binder <- res[[binder_col]] == T
  nonbinder <- res[[binder_col]] == F
  
  auc_data <- lapply(auc_thresholds, function(thresh){
    
    # get cells with min nhood p below threshold
    da.cell.list <- res[[p_col]] < thresh
    
    true_pos <- sum(da.cell.list == T & binder)
    false_neg <- sum(da.cell.list == F & binder)
    true_neg <- sum(da.cell.list == F & nonbinder)
    false_pos <- sum(da.cell.list == T & nonbinder)
    
    # auroc
    TPR <- true_pos / (true_pos + false_neg)
    FPR <- 1 - (true_neg / (true_neg + false_pos))

    # FDR
    FDR <- false_pos / (false_pos + true_pos)

    # AUPRC
    precision <- true_pos / (false_pos + true_pos)
    
    return(data.frame('TPR' = TPR,
                      'FPR' = FPR,
                      'Precision' = precision,
                      'FDR' = FDR,
                      'TP' = true_pos,
                      'FP' = false_pos,
                      'TN' = true_neg,
                      'FN' = false_neg))
    
    
  })
  
  auc_df <- do.call(rbind, auc_data)
  auc_df$p_value <- auc_thresholds

  # estimate the first precision point - it should always be NA b/c no false or true positives below the first threshold
  if (is.na(auc_df[1,'Precision']) & nrow(auc_df) > 1){
    auc_df[1,'Precision'] <- auc_df[2,'Precision']
  }
  
  # get auroc
  auroc <- pracma::trapz(auc_df$FPR, auc_df$TPR)
  title_auroc <- paste0(tool, '\nAUROC: ', round(auroc, 2))

  # get auprc
  auprc <- pracma::trapz(auc_df$TPR, auc_df$Precision)
  title_auprc <- paste0(tool, '\nAUPRC: ', round(auprc, 2))

  if (invalid_seqs > 0){
    subtitle <- paste0(prettyNum(sum(valid_seqs), big.mark = ",", scientific = FALSE), '/', 
                       prettyNum(total_seqs, big.mark = ",", scientific = FALSE), ' seqs clustered')
  } else{
    subtitle <- ''
  }
  
  p_auroc <- auc_df %>%
    ggplot(aes(x = FPR, y = TPR)) +
    geom_abline(slope = 1, intercept = 0, color = 'gray') +
    geom_point(size = 1.25) +
    geom_line() +
    labs(title = title_auroc,
         subtitle = subtitle,
         x = 'False Positive Rate',
         y = 'True Positive Rate') + 
    theme_minimal(base_size = 18) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 12, margin = margin(t = -5)))

  p_auprc <- auc_df %>%
    ggplot(aes(x = TPR, y = Precision)) +
    geom_hline(yintercept = auprc_baseline, color = 'red', linetype = 'dashed') +
    geom_point(size = 1.25) +
    geom_line() +
    labs(title = title_auprc,
         subtitle = subtitle,
         x = 'Recall',
         y = 'Precision') + 
    theme_minimal(base_size = 18) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 12, margin = margin(t = -5))) +
    scale_y_continuous(limits = c(0, 1))
  
  return(list('auroc' = auroc, 'plot_auroc' = p_auroc, 
              'auprc' = auprc, 'plot_auprc' = p_auprc,
              'table' = auc_df))
  
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

####################
### load results ###
####################

# dependent on the input tool

##### CDR3 SIM #####
if (TOOL ==  'cdr3_similarity'){
    # gather ASC files and create combined results file
    message('Combining ASC results for CDR3 similarity...')

    res_files <- lapply(MD_LOCS, function(f){
      df <- readr::read_tsv(f, show_col_types = F)

    asc_id <- stringr::str_split_i(basename(f), '_seq_summary', 1)

      df <- df %>%
        dplyr::mutate(ASC = asc_id,
                      convergent_clone_id_full = paste0(ASC, '_', convergent_clone_id))

      return(df)
    })
    
    cdr3_sim <- do.call(rbind, res_files)
    
    cdr3_sim <- cdr3_sim %>%
      dplyr::rename(sequence_id = id_col)

    # get new corrected p-values
    cluster_results <- cdr3_sim[c('convergent_clone_id_full', 'p_value_fisher', 'p_value_wilcox')] %>% distinct()

    # get correction for aggregated result
    cluster_results$agg_fdr_fisher <- p.adjust(cluster_results$p_value_fisher, method = 'BH')
    cluster_results$agg_fdr_wilcox <- p.adjust(cluster_results$p_value_wilcox, method = 'BH')

    # add back
    cdr3_sim <- cdr3_sim %>%
      dplyr::left_join(cluster_results[c('convergent_clone_id_full', 'agg_fdr_fisher', 'agg_fdr_wilcox')], 
                       by = 'convergent_clone_id_full',
                       relationship = 'many-to-one')
    
    readr::write_tsv(cdr3_sim, file.path('tables', 'combined_res.tsv.gz'))

    if (AUC_VAR != FALSE){
        
        message('Calculating CDRH3 Similarity AUC curves...')

        p_cdr3_fisher_asc <- evaluation_curve(cdr3_sim, 'p_value_fisher', AUC_VAR, tool = 'CDRH3 Similarity + Fisher')

        ggsave(file.path('figures', 'CDR3SIM_ASC_FISHER_AUROC.png'),
            p_cdr3_fisher_asc$plot_auroc,
            device = 'png',
            width = 7,
            height = 6)
        
        ggsave(file.path('figures', 'CDR3SIM_ASC_FISHER_AUPRC.png'),
            p_cdr3_fisher_asc$plot_auprc,
            device = 'png',
            width = 7,
            height = 6)

        write.table(p_cdr3_fisher_asc$table, file.path('tables', 'CDR3SIM_ASC_FISHER_EVALUATION.tsv'), sep = '\t',
                    row.names = F, quote = F)

        ###

        p_cdr3_wilcox_asc <- evaluation_curve(cdr3_sim, 'p_value_wilcox', AUC_VAR, tool = 'CDRH3 Similarity + Wilcoxon')

        ggsave(file.path('figures', 'CDR3SIM_ASC_WILCOX_AUROC.png'),
            p_cdr3_wilcox_asc$plot_auroc,
            device = 'png',
            width = 7,
            height = 6)

        ggsave(file.path('figures', 'CDR3SIM_ASC_WILCOX_AUPRC.png'),
            p_cdr3_wilcox_asc$plot_auprc,
            device = 'png',
            width = 7,
            height = 6)

        write.table(p_cdr3_wilcox_asc$table, file.path('tables', 'CDR3SIM_ASC_WILCOX_EVALUATION.tsv'), sep = '\t',
                    row.names = F, quote = F)

        all_auc_res <- data.frame(tool = c('CDRH3 Similarity + Fisher',
                                           'CDRH3 Similarity + Wilcoxon'),
                                  AUROC = c(p_cdr3_fisher_asc$auroc,
                                            p_cdr3_wilcox_asc$auroc),
                                  AUPRC = c(p_cdr3_fisher_asc$auprc,
                                            p_cdr3_wilcox_asc$auprc),
                                  FDR = c(calc_FDR(cdr3_sim, 'agg_fdr_fisher', AUC_VAR, 0.05),
                                          calc_FDR(cdr3_sim, 'agg_fdr_wilcox', AUC_VAR, 0.05)))

        write_tsv(all_auc_res, file.path('tables', 'ASC_evaluation_summary.tsv'))
    }
}

##### DASEQ #####
if (TOOL == 'DA-seq'){
    # gather ASC files and create combined results file
    message('Combining ASC results for DA-seq')
    
    res_files <- lapply(MD_LOCS, function(f){
      df <- readr::read_tsv(f, show_col_types = F)

      asc_id <- stringr::str_split_i(basename(f), '_da_seqs', 1)

      df <- df %>%
        dplyr::mutate(ASC = asc_id,
                      da.region.label.full = paste0(ASC, '_', da.region.label))

      return(df)
    })
    
    daseq <- do.call(rbind, res_files)

    # get new corrected p-values
    cluster_results <- daseq[c('da.region.label.full', 'pval.wilcoxon', 'pval.ttest', 'p_value_fisher', 'p_value_wilcox_onesided')] %>% distinct()

    # get correction for aggregated result
    cluster_results$agg_fdr_fisher <- p.adjust(cluster_results$p_value_fisher, method = 'BH')
    cluster_results$agg_fdr_wilcox <- p.adjust(cluster_results$pval.wilcoxon, method = 'BH')
    cluster_results$agg_fdr_ttest <- p.adjust(cluster_results$pval.ttest, method = 'BH')
    cluster_results$agg_fdr_wilco_onesided <- p.adjust(cluster_results$p_value_wilcox_onesided, method = 'BH')

    # add back
    daseq <- daseq %>%
      dplyr::left_join(cluster_results[c('da.region.label.full', 'agg_fdr_fisher', 'agg_fdr_wilcox', 
                                        'agg_fdr_ttest', 'agg_fdr_wilco_onesided')], 
                       by = 'da.region.label.full',
                       relationship = 'many-to-one')
    
    readr::write_tsv(daseq, file.path('tables', 'combined_res.tsv.gz'))

    if (AUC_VAR != FALSE){
        message('Calculating DA-Seq AUC curves...')

        p_daseq <- evaluation_curve(daseq, 'wilcox.adj.BH', AUC_VAR, tool = 'DA-seq')

        ggsave(file.path('figures', 'DAseq_ASC_AUROC.png'),
               p_daseq$plot_auroc,
               device = 'png',
               width = 7,
               height = 6)
            
        ggsave(file.path('figures', 'DAseq_ASC_AUPRC.png'),
               p_daseq$plot_auprc,
               device = 'png',
               width = 7,
               height = 6)

        write.table(p_daseq$table, file.path('tables', 'DAseq_ASC_EVALUATION.tsv'), sep = '\t',
                    row.names = F, quote = F)

        ###

        p_daseq_onesided <- evaluation_curve(daseq, 'fdr_wilcox_onesided', AUC_VAR, tool = 'DA-seq + OS Wilcoxon')

        ggsave(file.path('figures', 'DAseq_ASC_OS_WILCOX_AUROC.png'),
            p_daseq_onesided$plot_auroc,
            device = 'png',
            width = 7,
            height = 6)

        ggsave(file.path('figures', 'DAseq_ASC_OS_WILCOX_AUPRC.png'),
            p_daseq_onesided$plot_auprc,
            device = 'png',
            width = 7,
            height = 6)

        write.table(p_daseq_onesided$table, file.path('tables', 'DAseq_ASC_OS_WILCOX_EVALUATION.tsv'), sep = '\t',
                    row.names = F, quote = F)

        ###

        p_daseq_fisher <- evaluation_curve(daseq, 'p_value_fisher', AUC_VAR, tool = 'DA-seq + Fisher')

        ggsave(file.path('figures', 'DAseq_ASC_FISHER_AUROC.png'),
            p_daseq_fisher$plot_auroc,
            device = 'png',
            width = 7,
            height = 6)

        ggsave(file.path('figures', 'DAseq_ASC_FISHER_AUPRC.png'),
            p_daseq_fisher$plot_auprc,
            device = 'png',
            width = 7,
            height = 6)

        write.table(p_daseq_fisher$table, file.path('tables', 'DAseq_ASC_FISHER_EVALUATION.tsv'), sep = '\t',
                    row.names = F, quote = F)

        all_auc_res <- data.frame(tool = c('DA-Seq',
                                           'DA-Seq + Fisher',
                                           'DA-Seq + One-sided Wilcoxon'),
                                  AUROC = c(p_daseq$auroc,
                                            p_daseq_fisher$auroc,
                                            p_daseq_onesided$auroc),
                                  AUPRC = c(p_daseq$auprc,
                                            p_daseq_fisher$auprc,
                                            p_daseq_onesided$auprc),
                                  FDR = c(calc_FDR(daseq, 'agg_fdr_wilcox', AUC_VAR, 0.05),
                                          calc_FDR(daseq, 'agg_fdr_fisher', AUC_VAR, 0.05),
                                          calc_FDR(daseq, 'agg_fdr_wilco_onesided', AUC_VAR, 0.05)))

        write_tsv(all_auc_res, file.path('tables', 'ASC_evaluation_summary.tsv'))
    }
}

##### MILO #####
if (TOOL == 'Milo'){
    # gather ASC files and create combined results file
    message('Combining ASC results for Milo')
    res_files <- lapply(MD_LOCS, function(f){

      df <- readr::read_tsv(f, show_col_types = F)

     asc_id <- stringr::str_split_i(basename(f), '_seq_results', 1)

      df <- df %>%
        dplyr::mutate(ASC = asc_id,
                      nhood_id = paste0(ASC, '_', nhood_id))

      return(df)
    })
    
    milo <- do.call(rbind, res_files)
    
    milo <- milo %>%
      dplyr::rename(sequence_id = id_col)

    # get new corrected p-values on ALL nhoods (not just ones in min_nhood_ID)
    cluster_results <- milo %>%
                        dplyr::filter(!is.na(nhood_id)) %>%
                        dplyr::select(c('nhood_id', 'PValue')) %>% 
                        distinct()

    # get correction for aggregated result
    cluster_results$agg_fdr <- p.adjust(cluster_results$PValue, method = 'BH')

    # add back
    milo <- milo %>%
      dplyr::left_join(cluster_results[c('nhood_id', 'agg_fdr')], 
                       by = 'nhood_id', relationship = 'many-to-one')
    
    readr::write_tsv(milo, file.path('tables', 'combined_res.tsv.gz'))

    if (AUC_VAR != FALSE){
        message('Calculating Milo AUC curve...')

        p_milo <- evaluation_curve(milo %>% dplyr::filter(!is.na(sequence_id)), 'min_nhood_FDR', 
                                  AUC_VAR, tool = 'Milo', simplify_p = T)

        ggsave(file.path('figures', 'Milo_ASC_AUROC.png'),
                p_milo$plot_auroc,
                device = 'png',
                width = 7,
                height = 6)

        ggsave(file.path('figures', 'Milo_ASC_AUPRC.png'),
                p_milo$plot_auprc,
                device = 'png',
                width = 7,
                height = 6)

        write.table(p_milo$table, file.path('tables', 'Milo_ASC_EVALUATION.tsv'), sep = '\t',
                    row.names = F, quote = F)
        
        all_auc_res <- data.frame(tool = c('Milo'),
                                  AUROC = c(p_milo$auroc),
                                  AUPRC = c(p_milo$auprc),
                                  FDR = c(calc_FDR(milo %>% dplyr::filter(!is.na(sequence_id)), 
                                                  'agg_fdr', AUC_VAR, 0.05)))

        write_tsv(all_auc_res, file.path('tables', 'ASC_evaluation_summary.tsv'))
    }
} 

##### BCRdist #####
if (TOOL == 'BCRdist'){
    # gather ASC files and create combined results file
    message('Combining ASC results for BCRdist...')
    
    res_files <- lapply(MD_LOCS, function(f){
      df <- readr::read_tsv(f, show_col_types = F)

      asc_id <- stringr::str_split_i(basename(f), '_seq_summary', 1)

      df <- df %>%
        dplyr::mutate(ASC = asc_id,
                      convergent_clone_id_full = paste0(ASC, '_', convergent_clone_id))

      return(df)
    })
    
    bcrdist <- do.call(rbind, res_files)
    
    bcrdist <- bcrdist %>%
      dplyr::rename(sequence_id = id_col)
    
    # get new corrected p-values
    cluster_results <- bcrdist[c('convergent_clone_id_full', 'p_value_fisher', 'p_value_wilcox')] %>% distinct()

    # get correction for aggregated result
    cluster_results$agg_fdr_fisher <- p.adjust(cluster_results$p_value_fisher, method = 'BH')
    cluster_results$agg_fdr_wilcox <- p.adjust(cluster_results$p_value_wilcox, method = 'BH')

    # add back
    bcrdist <- bcrdist %>%
      dplyr::left_join(cluster_results[c('convergent_clone_id_full', 'agg_fdr_fisher', 'agg_fdr_wilcox')], 
                       by = 'convergent_clone_id_full',
                       relationship = 'many-to-one')
    
    readr::write_tsv(bcrdist, file.path('tables', 'combined_res.tsv.gz'))

    if (AUC_VAR != FALSE){
        message('Calculating BCRdist AUC curves...')

        p_bcrdist_fisher_asc <- evaluation_curve(bcrdist, 'p_value_fisher', AUC_VAR, tool = 'BCRdist + Fisher')

        ggsave(file.path('figures', 'BCRdist_ASC_FISHER_AUROC.png'),
            p_bcrdist_fisher_asc$plot_auroc,
            device = 'png',
            width = 7,
            height = 6)

        ggsave(file.path('figures', 'BCRdist_ASC_FISHER_AUPRC.png'),
            p_bcrdist_fisher_asc$plot_auprc,
            device = 'png',
            width = 7,
            height = 6)

        write.table(p_bcrdist_fisher_asc$table, file.path('tables', 'BCRdist_ASC_FISHER_EVALUATION.tsv'), sep = '\t',
                    row.names = F, quote = F)

        ###

        p_bcrdist_wilcox_asc <- evaluation_curve(bcrdist, 'p_value_wilcox', AUC_VAR, tool = 'BCRdist + Wilcoxon')

        ggsave(file.path('figures', 'BCRdist_ASC_WILCOX_AUROC.png'),
            p_bcrdist_wilcox_asc$plot_auroc,
            device = 'png',
            width = 7,
            height = 6)

        ggsave(file.path('figures', 'BCRdist_ASC_WILCOX_AUPRC.png'),
            p_bcrdist_wilcox_asc$plot_auprc,
            device = 'png',
            width = 7,
            height = 6)

        write.table(p_bcrdist_wilcox_asc$table, file.path('tables', 'BCRdist_ASC_WILCOX_EVALUATION.tsv'), sep = '\t',
                    row.names = F, quote = F)

        all_auc_res <- data.frame(tool = c('BCRdist + Fisher',
                                           'BCRdist + Wilcoxon'),
                                  AUROC = c(p_bcrdist_fisher_asc$auroc,
                                            p_bcrdist_wilcox_asc$auroc),
                                  AUPRC = c(p_bcrdist_fisher_asc$auprc,
                                            p_bcrdist_wilcox_asc$auprc),
                                  FDR = c(calc_FDR(bcrdist, 'agg_fdr_fisher', AUC_VAR, 0.05),
                                          calc_FDR(bcrdist, 'agg_fdr_wilcox', AUC_VAR, 0.05)))

        write_tsv(all_auc_res, file.path('tables', 'ASC_evaluation_summary.tsv'))
    }
}

message('Analysis finished.')
