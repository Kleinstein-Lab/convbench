process GET_ASC_AUROC{
    tag "${meta_id}_${tool_id}"
    label 'process_medium'

    container "docker.io/cfsullivan16/asc_auroc:1.0.0dev"

    input:
    tuple val(tool_id), val(meta_id), path(seq_tables)

    output:
    path "tables/combined_res.tsv.gz", emit: combined_res
    path "tables/*.tsv", emit: auc_files
    path "figures/*.png", emit: figs

    script:
    """
    get_asc_auroc.R \
    -md ${seq_tables.join(',')} \
    -a ${params.auc_variable} \
    -t $tool_id

    """
}
