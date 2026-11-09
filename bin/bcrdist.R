#!/usr/bin/env Rscript

message(paste0('Starting run: ', Sys.time()))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(airr)
  library(alakazam)
  library(argparse)
  library(reticulate)
  library(arrow)
  library(tidyr)
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

parser <- ArgumentParser(description = "Data location and BCRdist hyperparameters.")

parser$add_argument('-he', '--helper_loc', type = 'character', default = 'bin/tcrdist3.py',
                    help = 'File path for the helper function location.')

parser$add_argument('-md', '--metadata_loc', type = 'character', default = 'metadata',
                    help = 'File path for the metadata location.')

parser$add_argument('-g', '--group_id', type = 'character', default = 'all',
                    help = 'Group ID (i.e. the ASC) used to prefix cluster IDs and output files.')

parser$add_argument('-o', '--output_dir', type = 'character', default = 'tcrdist3_output',
                    help = 'Specify an output directory location.')

parser$add_argument('-t', '--cluster_threshold', type = 'double', default = '60',
                    help = 'The distance threshold for forming clusters with hierarchical clones.')

# Not enabled currently
# parser$add_argument('-l', '--linkage_method', type = 'character', default = 'single',
#                    help = 'The linkage method to be used in forming clusters.')

# Not enabled currently
# parser$add_argument('-sc', '--single_cell', type = 'character', default = "FALSE",
#                    help = 'Input true if V(D)J info is present and contains paired heavy and light chain info.')

parser$add_argument('-a', '--auc_var', type = 'character', default = FALSE,
                    help = 'Specify which column contains the ground truth for evaluation (i.e. "simulated" or "binder"). Column type should be logical. If no AUC variable, set to FALSE.')

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
GROUP_ID <- args$group_id

OUTPUT_DIR <- args$output_dir

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

end_time <- Sys.time()
time_taken <- end_time - start_time

##################################
### SAVE INPUTS FOR EVALUATION ###
##################################

convergent_clone <- convergent_clone %>%
  dplyr::mutate(sequence_id = id_col,
                cluster_id = paste0(GROUP_ID, '_', convergent_clone_id))

write_evaluation_inputs(convergent_clone, GROUP_ID, AUC_VAR, time_taken)

message(paste0(dplyr::n_distinct(convergent_clone$cluster_id), ' clusters saved for evaluation.'))

message(paste0('Ending run: ', Sys.time()))

sessionInfo()
