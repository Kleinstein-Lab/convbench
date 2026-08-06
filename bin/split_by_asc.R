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

parser <- ArgumentParser(description = "Data location and Mal-ID algorithm hyperparameters.")

parser$add_argument('-a', '--asc_guide', type = 'character', default = 'asc_guide.tsv',
                    help = 'File path to a tab-separated file containing IMGT to ASC allele translations.')

# parser$add_argument('-l', '--label', type = 'character', default = 'COVID',
#                     help = 'Metadata label for a given AIRR/embedding pair.')

parser$add_argument('-d', '--data_loc', type = 'character', default = 'data',
                    help = 'File path for the embedding or RNA-Seq data location.')

parser$add_argument('-md', '--metadata_loc', type = 'character', default = 'metadata',
                    help = 'File path for the metadata location. Metadata and data files should have 1:1 matching sequence identifiers.')

parser$add_argument('-e', '--use_embedding', type = 'logical', default = TRUE,
                    help = 'Specify whether embedding AND metadata (TRUE) or ONLY metadata (FALSE) should be accounted for in ASC splitting.')

# Parse the arguments
args <- parser$parse_args()

ASC_GUIDE_LOC <- args$asc_guide
MD_LOC <- args$metadata_loc
DATA_LOC <- args$data_loc
USE_EMB <- as.logical(args$use_embedding)
# META_LABEL <- args$label

# if (!dir.exists(file.path('ASCs'))){
#   dir.create(file.path('ASCs'))
# }

###################
### ASSIGN ASCs ###
###################

# load the ASC guide
tryCatch(
  
  {
    asc_guide <- airr::read_rearrangement(ASC_GUIDE_LOC)
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

# load the data
# metadata
message(paste0('Loading metadata: ', MD_LOC))

tryCatch(
  
  {
    md <- readr::read_tsv(MD_LOC)
  }, error = function(e){
    
    stop(e)
    
  }
  
)

# load embeddings or expr data
message(paste0('Loading data: ', DATA_LOC))

if (USE_EMB){

    tryCatch(
  
  {
    data <- data.table::fread(DATA_LOC, sep = '\t', header = T)
  }, error = function(e){
    
    stop(e)
    
  }
)

}

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
small_ASCs <- table(md$asc_group)[table(md$asc_group) < 200] 
small_ASCs <- names(small_ASCs)

# do not continue if it just makes one group
if (length(small_ASCs) == 1){
  stop('Execution halted: only 1 ASC is generated. Run without asc_mode enabled instead.')
}

md[md$asc_group %in% small_ASCs, 'asc_group'] <- 'small-ASC'

ASC_split <- split(md, md$asc_group)

# initialize manifest
# n_rows <- length(ASC_split)
# manifest <- data.frame(asc_id = character(n_rows),
#                        airr = character(n_rows),
#                        embedding = character(n_rows))

# row.names(manifest) <- names(ASC_split)

for (asc_choice in names(ASC_split)){
  
    message(paste0('Saving info for ASC ', asc_choice, '...'))
    # manifest[asc_choice, 'asc_id'] <- paste0(asc_choice, '_', META_LABEL)

    # if (!dir.exists(file.path('ASCs', asc_choice))){
    #   dir.create(file.path('ASCs', asc_choice))
    # }

    asc_md <- ASC_split[[asc_choice]]
    asc_md$ASC <- asc_choice
  
    write_tsv(asc_md,
             file.path(paste0(asc_choice, '_md.tsv.gz')))

    # manifest[asc_choice, 'airr'] <- file.path(asc_choice, paste0(asc_choice, '_', META_LABEL, '_md.tsv.gz'))

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
  
  # include a filename regardless
  # manifest[asc_choice, 'embedding'] <- file.path(asc_choice, paste0(asc_choice, '_', META_LABEL, '_emb.tsv.gz'))

  message(paste0(nrow(ASC_split[[asc_choice]]), ' sequences recorded for ', asc_choice, '.'))
  
}

# write_csv(manifest, file.path('ASCs', 'manifest.csv'))

message('Finished.')

