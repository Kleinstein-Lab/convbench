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

auc_curve <- function(res, p_col, binder_col, tool = '', simplify_p = F){
  
  invalid_seqs <- sum(is.na(res[[p_col]]))
  total_seqs <- nrow(res)
  valid_seqs <- total_seqs - invalid_seqs
  
  # change the sequences with no nhood to a min p of 1
  res[is.na(res[[p_col]]), p_col] <- 1
  
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
  
  # get binder info
  binder <- res[[binder_col]] == T
  nonbinder <- res[[binder_col]] == F
  
  auc_data <- lapply(auc_thresholds, function(thresh){
    
    # get cells with min nhood p below threshold
    da.cell.list <- res[[p_col]] < thresh
    
    true_pos <- sum(da.cell.list == T & binder)
    false_neg <- sum(da.cell.list == F & binder)
    true_neg <- sum(da.cell.list == F & nonbinder)
    false_pos <- sum(da.cell.list == T & nonbinder)
    
    TPR <- true_pos / (true_pos + false_neg)
    FPR <- 1 - (true_neg / (true_neg + false_pos))
    
    return(data.frame('TPR' = TPR,
                      'FPR' = FPR))
    
    
  })
  
  auc_df <- do.call(rbind, auc_data)
  auc_df$p_value <- auc_thresholds
  
  # get auroc
  auroc <- pracma::trapz(auc_df$FPR, auc_df$TPR)
  title <- paste0(tool, '\nAUC: ', round(auroc, 2))
  
  if (invalid_seqs > 0){
    subtitle <- paste0(prettyNum(sum(valid_seqs), big.mark = ",", scientific = FALSE), '/', 
                       prettyNum(total_seqs, big.mark = ",", scientific = FALSE), ' seqs clustered')
  } else{
    subtitle <- ''
  }
  
  p <- auc_df %>%
    ggplot(aes(x = FPR, y = TPR)) +
    geom_point(size = 1.25) +
    geom_line() +
    labs(title = title,
         subtitle = subtitle,
         x = 'False Positive Rate',
         y = 'True Positive Rate') + 
    theme_minimal(base_size = 18) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 12, margin = margin(t = -5)))
  
  return(list('auroc' = auroc, 'plot' = p, 'table' = auc_df))
  
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
    
    readr::write_tsv(cdr3_sim, file.path('tables', 'combined_res.tsv.gz'))

    if (AUC_VAR != FALSE){
        
        message('Calculating CDRH3 Similarity AUC curves...')

        p_cdr3_fisher_asc <- auc_curve(cdr3_sim, 'p_value_fisher', AUC_VAR, tool = 'CDRH3 Similarity + Fisher')

        ggsave(file.path('figures', 'CDR3SIM_ASC_FISHER_AUC_curve.png'),
            p_cdr3_fisher_asc$plot,
            device = 'png',
            width = 7,
            height = 6)

        write.table(p_cdr3_fisher_asc$table, file.path('tables', 'CDR3SIM_ASC_FISHER_AUC.tsv'), sep = '\t',
                    row.names = F, quote = F)

        ###

        p_cdr3_wilcox_asc <- auc_curve(cdr3_sim, 'p_value_wilcox', AUC_VAR, tool = 'CDRH3 Similarity + Wilcoxon')

        ggsave(file.path('figures', 'CDR3SIM_ASC_WILCOX_AUC_curve.png'),
            p_cdr3_wilcox_asc$plot,
            device = 'png',
            width = 7,
            height = 6)

        write.table(p_cdr3_wilcox_asc$table, file.path('tables', 'CDR3SIM_ASC_WILCOX_AUC.tsv'), sep = '\t',
                    row.names = F, quote = F)

        all_auc_res <- data.frame(tool = c('CDRH3 Similarity + Fisher',
                                        'CDRH3 Similarity + Wilcoxon'),
                                AUC = c(p_cdr3_fisher_asc$auroc,
                                        p_cdr3_wilcox_asc$auroc))

        write_tsv(all_auc_res, file.path('tables', 'AUC_summary.tsv'))
    }
}

##### DASEQ #####
if (TOOL == 'DA-seq'){
    # gather ASC files and create combined results file
    message('Combining ASC results for DA-seq')
    
    res_files <- lapply(MD_LOCS, function(f){
      df <- readr::read_tsv(f, show_col_types = F)

      asc_id <- stringr::str_split_i(basename(f), '_seq_summary', 1)

      df <- df %>%
        dplyr::mutate(ASC = asc_id,
                      da.region.label.full = paste0(ASC, '_', da.region.label))

      return(df)
    })
    
    daseq <- do.call(rbind, res_files)
    
    readr::write_tsv(daseq, file.path('tables', 'combined_res.tsv.gz'))

    if (AUC_VAR != FALSE){
        message('Calculating DA-Seq AUC curves...')

        p_daseq <- auc_curve(daseq, 'wilcox.adj.BH', AUC_VAR, tool = 'DA-seq')

        ggsave(file.path('figures', 'DAseq_ASC_AUC_curve.png'),
            p_daseq$plot,
            device = 'png',
            width = 7,
            height = 6)

        write.table(p_daseq$table, file.path('tables', 'DAseq_ASC_AUC.tsv'), sep = '\t',
                    row.names = F, quote = F)

        ###

        p_daseq_onesided <- auc_curve(daseq, 'fdr_wilcox_onesided', AUC_VAR, tool = 'DA-seq + OS Wilcoxon')

        ggsave(file.path('figures', 'DAseq_ASC_OS_WILCOX_AUC_curve.png'),
            p_daseq_onesided$plot,
            device = 'png',
            width = 7,
            height = 6)

        write.table(p_daseq_onesided$table, file.path('tables', 'DAseq_ASC_OS_WILCOX_AUC.tsv'), sep = '\t',
                    row.names = F, quote = F)

        ###

        p_daseq_fisher <- auc_curve(daseq, 'p_value_fisher', AUC_VAR, tool = 'DA-seq + Fisher')

        ggsave(file.path('figures', 'DAseq_ASC_FISHER_AUC_curve.png'),
            p_daseq_fisher$plot,
            device = 'png',
            width = 7,
            height = 6)

        write.table(p_daseq_fisher$table, file.path('tables', 'DAseq_ASC_FISHER_AUC.tsv'), sep = '\t',
                    row.names = F, quote = F)

        all_auc_res <- data.frame(tool = c('DA-Seq',
                                        'DA-Seq + Fisher',
                                        'DA-Seq + One-sided Wilcoxon'),
                                AUC = c(p_daseq$auroc,
                                        p_daseq_fisher$auroc,
                                        p_daseq_onesided$auroc))

        write_tsv(all_auc_res, file.path('tables', 'AUC_summary.tsv'))
    }
}

##### MILO #####
if (TOOL == 'Milo'){
    # gather ASC files and create combined results file
    message('Combining ASC results for Milo')
    res_files <- lapply(MD_LOCS, function(f){

      df <- readr::read_tsv(f, show_col_types = F)

     asc_id <- stringr::str_split_i(basename(f), '_seq_summary', 1)

      df <- df %>%
        dplyr::mutate(ASC = asc_id,
                      min_nhood_id = paste0(ASC, '_', min_nhood_id))

      return(df)
    })
    
    milo <- do.call(rbind, res_files)
    
    milo <- milo %>%
      dplyr::rename(sequence_id = id_col)
    
    readr::write_tsv(milo, file.path('tables', 'combined_res.tsv.gz'))

    if (AUC_VAR != FALSE){
        message('Calculating Milo AUC curve...')

        p_milo <- auc_curve(milo, 'min_nhood_FDR', 
                            AUC_VAR, tool = 'Milo', simplify_p = T)

        ggsave(file.path('figures', 'Milo_ASC_AUC_curve.png'),
                p_milo$plot,
                device = 'png',
                width = 7,
                height = 6)

        write.table(p_milo$table, file.path('tables', 'Milo_ASC_AUC.tsv'), sep = '\t',
                    row.names = F, quote = F)
        
        all_auc_res <- data.frame(tool = c('Milo'),
                                AUC = c(p_milo$auroc))

        write_tsv(all_auc_res, file.path('tables', 'AUC_summary.tsv'))
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
    
    readr::write_tsv(bcrdist, file.path('tables', 'combined_res.tsv.gz'))

    if (AUC_VAR != FALSE){
        message('Calculating BCRdist AUC curves...')

        p_bcrdist_fisher_asc <- auc_curve(bcrdist, 'p_value_fisher', AUC_VAR, tool = 'BCRdist + Fisher')

        ggsave(file.path('figures', 'BCRdist_ASC_FISHER_AUC_curve.png'),
            p_bcrdist_fisher_asc$plot,
            device = 'png',
            width = 7,
            height = 6)

        write.table(p_bcrdist_fisher_asc$table, file.path('tables', 'BCRdist_ASC_FISHER_AUC.tsv'), sep = '\t',
                    row.names = F, quote = F)

        ###

        p_bcrdist_wilcox_asc <- auc_curve(bcrdist, 'p_value_wilcox', AUC_VAR, tool = 'BCRdist + Wilcoxon')

        ggsave(file.path('figures', 'BCRdist_ASC_WILCOX_AUC_curve.png'),
            p_bcrdist_wilcox_asc$plot,
            device = 'png',
            width = 7,
            height = 6)

        write.table(p_bcrdist_wilcox_asc$table, file.path('tables', 'BCRdist_ASC_WILCOX_AUC.tsv'), sep = '\t',
                    row.names = F, quote = F)

        all_auc_res <- data.frame(tool = c('BCRdist + Fisher',
                                        'BCRdist + Wilcoxon'),
                                AUC = c(p_bcrdist_fisher_asc$auroc,
                                        p_bcrdist_wilcox_asc$auroc))

        write_tsv(all_auc_res, file.path('tables', 'AUC_summary.tsv'))
    }
}

message('Analysis finished.')
