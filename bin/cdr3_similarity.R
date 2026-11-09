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
  library(tidyr)
  library(argparse)
})


########################
### HELPER FUNCTIONS ###
########################

write_evaluation_inputs <- function(memberships, group_id, auc_variable, time_taken){
  # save the inputs for evaluate.R:
  #   <group_id>_cluster_counts.tsv = cluster x subject count matrix
  #   <group_id>_memberships.tsv = sequence-to-cluster memberships with ground truth (if available)
  #   <group_id>_tool_info.tsv = run time
  # memberships = one row per sequence with sequence_id, subject_id and cluster_id
  #               (cluster IDs already prefixed with group_id)

  memberships$group_id <- group_id

  membership_cols <- c('sequence_id', 'subject_id', 'group_id', 'cluster_id')
  if (auc_variable != FALSE){
    membership_cols <- c(membership_cols, auc_variable)
  }

  cluster_counts <- memberships %>%
    dplyr::filter(!is.na(cluster_id)) %>%
    dplyr::count(cluster_id, subject_id) %>%
    tidyr::pivot_wider(names_from = 'subject_id', values_from = 'n', values_fill = 0)

  write.table(cluster_counts, file.path(OUTPUT_DIR, 'tables', paste0(group_id, '_cluster_counts.tsv')),
              sep = '\t', quote = F, row.names = F)

  write.table(memberships[membership_cols], file.path(OUTPUT_DIR, 'tables', paste0(group_id, '_memberships.tsv')),
              sep = '\t', quote = F, row.names = F)

  write.table(data.frame(group_id = group_id, time_min = as.numeric(time_taken, units = 'mins')),
              file.path(OUTPUT_DIR, 'tables', paste0(group_id, '_tool_info.tsv')),
              sep = '\t', quote = F, row.names = F)
}

##############################
### SET UP THE ENVIRONMENT ###
##############################

parser <- ArgumentParser(description = "Data location and Mal-ID algorithm hyperparameters.")

parser$add_argument('-md', '--metadata_loc', type = 'character', default = 'metadata',
                    help = 'File path for the metadata location.')

parser$add_argument('-g', '--group_id', type = 'character', default = 'all',
                    help = 'Group ID (i.e. the ASC) used to prefix cluster IDs and output files.')

parser$add_argument('-o', '--output_dir', type = 'character', default = 'MalID_internal_output',
                    help = 'Specify an output directory location.')

parser$add_argument('-t', '--cluster_threshold', type = 'double', default = '0.15',
                    help = 'The distance threshold for forming clusters with hierarchical clones.')

parser$add_argument('-l', '--linkage_method', type = 'character', default = 'single',
                    help = 'The linkage method to be used in forming clusters.')

parser$add_argument('-c', '--cpu', type = 'integer', default = 16,
                    help = 'Specify the number of cpus used to do hierarchical clustering (default nproc = 16).')

parser$add_argument('-a', '--auc_var', type = 'character', default = FALSE,
                    help = 'Specify which column contains the ground truth for evaluation (i.e. "simulated" or "binder"). Column type should be logical. If no AUC variable, set to FALSE.')

parser$add_argument('-v', '--vdj_info', type = 'logical', default = TRUE,
                    help = 'Is v call and j call information included in the metadata? Can apply to expression or embedding data.')

parser$add_argument('-sc', '--single_cell', type = 'logical', default = FALSE,
                    help = 'Input TRUE if V(D)J info is present and contains paired heavy and light chain info.')

parser$add_argument('-r', '--remove_dups', type = 'logical', default = FALSE,
                    help = 'Will remove duplicate embeddings within an individual if TRUE.')

parser$add_argument('-f', '--first', type = 'logical', default = FALSE,
                    help = 'Specify whether to consider only the first occurrence of each clone.')


# Parse the arguments
args <- parser$parse_args()

MD_LOC <- args$metadata_loc
GROUP_ID <- args$group_id

OUTPUT_DIR <- args$output_dir

THRESH <- args$cluster_threshold
LINKAGE <- args$linkage_method
CPUS <- args$cpu

VDJ <- args$vdj_info
SINGLE_CELL <- args$single_cell
AUC_VAR <- args$auc_var
REMOVE_DUPS <- args$remove_dups
FIRST <- args$first

if (AUC_VAR != FALSE){
  message(paste0('AUC variable ', AUC_VAR, ' will be saved for evaluation.'))
} else{
  message('No AUC variable provided.')
}

if (REMOVE_DUPS){
  message('Duplicate embeddings within a subject will be collapsed.')
}

# create locations for results to be saved within output dir
if(!dir.exists(file.path(OUTPUT_DIR))){
  dir.create(file.path(OUTPUT_DIR))
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
message(paste0('Starting hierarchical clustering with ', CPUS, ' processors.'))
message(paste0('Linkage method: ', LINKAGE))
message(paste0('Threshold: ', THRESH))

start_time <- Sys.time()

# STRATEGY IF WE HAVE A HIGHER SAMPLE SIZE
###############################################################################

###############
### CLUSTER ###
###############

# required for method="aa" to produce valid AA translations
md <- md %>%
  dplyr::filter(
    nchar(junction) > 0,
    nchar(junction) %% 3 == 0
  )

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
                                                  first=FIRST,
                                                  cdr3=FALSE, 
                                                  mod3=FALSE,
                                                  max_n=0, 
                                                  nproc=CPUS,
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
                                                  first=FIRST,
                                                  cdr3=FALSE, 
                                                  mod3=FALSE,
                                                  max_n=0, 
                                                  nproc=CPUS,
                                                  verbose=T, log=NULL,
                                                  summarize_clones=FALSE)
}

message('Hierarchical clustering complete.')

end_time <- Sys.time()
time_taken <- end_time - start_time

##################################
### SAVE INPUTS FOR EVALUATION ###
##################################

convergent_clones <- convergent_clones %>%
  dplyr::mutate(sequence_id = id_col,
                cluster_id = paste0(GROUP_ID, '_', convergent_clone_id))

write_evaluation_inputs(convergent_clones, GROUP_ID, AUC_VAR, time_taken)

message(paste0(dplyr::n_distinct(convergent_clones$cluster_id), ' clusters saved for evaluation.'))

message(paste0('Ending run: ', Sys.time()))

sessionInfo()
