#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(airr)
  library(alakazam)
  library(argparse)
})

##############################
### SET UP THE ENVIRONMENT ###
##############################

parser <- ArgumentParser(description = "Pre-process AIRR/embedding pairs: record library sizes and optionally split by ASC.")

parser$add_argument('-a', '--asc_guide', type = 'character', default = NULL,
                    help = 'File path to a tab-separated file containing IMGT to ASC allele translations. If not provided, data are not split by ASC.')

parser$add_argument('-d', '--data_loc', type = 'character', default = 'data',
                    help = 'File path for the embedding or RNA-Seq data location.')

parser$add_argument('-md', '--metadata_loc', type = 'character', default = 'metadata',
                    help = 'File path for the metadata location. Metadata and data files should have 1:1 matching sequence identifiers.')

parser$add_argument('-da', '--da_variable', type = 'character', default = 'status',
                    help = 'Stratification variable that should be used to determine differential abundance. There should be two levels in this factor/categorical variable.')

parser$add_argument('-e', '--use_embedding', type = 'logical', default = TRUE,
                    help = 'Specify whether embedding AND metadata (TRUE) or ONLY metadata (FALSE) should be accounted for in library sizes and ASC splitting.')

# Parse the arguments
args <- parser$parse_args()

ASC_GUIDE_LOC <- args$asc_guide
MD_LOC <- args$metadata_loc
DATA_LOC <- args$data_loc
USE_EMB <- as.logical(args$use_embedding)
DA_VAR <- args$da_variable

SPLIT_ASC <- !is.null(ASC_GUIDE_LOC)

if (SPLIT_ASC){
  message('ASC guide provided. Data will be split by ASC.')
} else{
  message('No ASC guide provided. Only library sizes will be recorded.')
}

#################
### LOAD DATA ###
#################

# metadata
message(paste0('Loading metadata: ', MD_LOC))

tryCatch(

  {
    md <- airr::read_rearrangement(MD_LOC)
  }, error = function(e){

    stop(e)

  }

)

if (!DA_VAR %in% colnames(md)){
  stop(paste0('DA variable ', DA_VAR, ' not found in metadata.'))
}

# load embeddings or expr data
if (USE_EMB){
  message(paste0('Loading data: ', DATA_LOC))

  tryCatch(

    {
      # only the sequence IDs are needed if the embedding is not being split
      if (SPLIT_ASC){
        data <- data.table::fread(DATA_LOC, sep = '\t', header = T)
      } else{
        data <- data.table::fread(DATA_LOC, sep = '\t', header = T, select = 'sequence_id')
      }
    }, error = function(e){

      stop(e)

    }
  )

  # get common sequences
  common_seqs <- intersect(md$sequence_id, data$sequence_id)

  message(paste0(length(common_seqs), ' of ', nrow(md), ' metadata sequences found in the data.'))

  md <- md %>%
    dplyr::filter(sequence_id %in% common_seqs) %>%
    as.data.frame(check.names = F)

  data <- data %>%
    dplyr::filter(sequence_id %in% common_seqs) %>%
    as.data.frame(check.names = F)
}

#########################
### GET LIBRARY SIZES ###
#########################

# save library size and da var information for future steps
lib_sizes <- md %>%
  dplyr::group_by(subject_id, !!sym(DA_VAR)) %>%
  dplyr::summarize(depth = n(), .groups = 'drop')

if (any(duplicated(lib_sizes$subject_id))){
  stop(paste0('Subjects found with more than one value of ', DA_VAR, ': ',
              paste(unique(lib_sizes$subject_id[duplicated(lib_sizes$subject_id)]), collapse = ', ')))
}

write_tsv(lib_sizes, 'library_sizes.tsv')

if (!SPLIT_ASC){
  message('Finished.')
  quit(save = 'no', status = 0)
}

###################
### ASSIGN ASCs ###
###################

# load the ASC guide
tryCatch(

  {
    asc_guide <- readr::read_tsv(ASC_GUIDE_LOC)
  }, error = function(e){

    stop(e)

  }

)

# expand parts where multiple names for alleles exist
asc_guide <- asc_guide %>%
  tidyr::separate_longer_delim(allele, delim = ";") %>%
  tidyr::separate_longer_delim(imgt_genes, delim = ";")

# get the ASC group
asc_guide$asc_group <- getGene(asc_guide$asc_allele, strip_d = F, omit_nl = F)

asc_allele_guide <- distinct(asc_guide[c('allele', 'asc_allele', 'asc_group')])

# add the v alleles
md$v_allele <- getAllele(md$v_call, strip_d = F, omit_nl = F)

# check if there are any V genes in the dataset that do not match up with ASCs
alleles <- unique(asc_guide$allele)
genes <- unique(asc_guide$imgt_genes)

not_in_table <- md$v_allele[!md$v_allele %in% alleles]

if (length(not_in_table) > 0){
    message('Alleles not present in the ASC guide:')
    print(unique(not_in_table))
} else{
    message('All alleles match to ASCs.')
}

# JOIN AT THE ALLELE LEVEL TO GET ASCs
md <- md %>%
  dplyr::left_join(asc_allele_guide, by = join_by(v_allele == allele))

message(paste0('Number of uncategorized sequences: ', sum(is.na(md$asc_group)), '.'))

#########################
### SPLIT UP THE DATA ###
#########################

ASCs <- unique(md$asc_group)

# combine ASCs under 1,000 sequences
message('Combining ASCs with under 1,000 sequences...')
small_ASCs <- table(md$asc_group)[table(md$asc_group) < 1000]
small_ASCs <- names(small_ASCs)

md[md$asc_group %in% small_ASCs, 'asc_group'] <- 'small-ASC'

# do not continue if it just makes one group
if (length(unique(md$asc_group)) == 1){
  stop('Execution halted: only 1 ASC is generated. Run without asc_mode enabled instead.')
}

ASC_split <- split(md, md$asc_group)

for (asc_choice in names(ASC_split)){

    message(paste0('Saving info for ASC ', asc_choice, '...'))

    asc_md <- ASC_split[[asc_choice]]
    asc_md$ASC <- asc_choice

    write_tsv(asc_md,
             file.path(paste0(asc_choice, '_md.tsv.gz')))

    if (USE_EMB){
    emb_filtered <- data %>%
        dplyr::filter(sequence_id %in% ASC_split[[asc_choice]]$sequence_id)

    write_tsv(emb_filtered,
                file.path(paste0(asc_choice, '_emb.tsv.gz')))
  } else{
    # write dummy file
    write_tsv(data.frame(sequence_id = character()),
                file.path(paste0(asc_choice, '_emb.tsv.gz')))
  }

  message(paste0(nrow(ASC_split[[asc_choice]]), ' sequences recorded for ', asc_choice, '.'))

}

message('Finished.')
