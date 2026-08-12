process SPLIT_BY_ASC{
    tag "${meta.id}"
    label 'process_medium'

    container "docker.io/cfsullivan16/split_by_asc:1.0.0dev"

    input:
    tuple val(meta), path(airr), path(embedding)
    path asc_guide

    output:
    tuple val(meta),
          path("*_md.tsv.gz"),
          path("*_emb.tsv.gz"),
          path("library_sizes.tsv")

    script:
    """
    split_by_asc.R \
    -a $asc_guide \
    -d $embedding \
    -md $airr \
    -da ${params.da_variable} \
    -e TRUE

    """
}