process SPLIT_BY_ASC{
    tag "${meta.id}"
    label 'process_medium'

    container "docker.io/cfsullivan16/split_by_asc:1.0.0dev"

    input:
    tuple val(meta), path(airr), path(embedding)
    path asc_guide

    output:
    tuple val(meta),
          path("ASCs/*")

    script:
    """
    split_by_asc.R \
    -a $asc_guide \
    -l ${meta.id} \
    -d $embedding \
    -md $airr \
    -e TRUE

    """
}