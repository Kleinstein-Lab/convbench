process PREPROCESS{
    tag "${meta.id}"
    label 'process_medium'

    container "docker.io/cfsullivan16/split_by_asc:1.0.0"

    input:
    tuple val(meta), path(airr), path(embedding)
    path asc_guide // optional: pass [] to skip ASC splitting

    output:
    tuple val(meta), path("library_sizes.tsv"), emit: library_sizes
    tuple val(meta), path("*_md.tsv.gz"), path("*_emb.tsv.gz"), emit: asc_files, optional: true

    script:
    def asc_guide_arg = asc_guide ? "-a ${asc_guide}" : ""
    """
    preprocess.R \
    ${asc_guide_arg} \
    -d $embedding \
    -md $airr \
    -da ${params.da_variable} \
    -e ${params.split_embeddings}

    """
}
