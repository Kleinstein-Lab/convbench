process GET_ASC_AUROC{
    tag "${meta_id}_${tool_id}"
    label 'process_medium'

    container "docker.io/cfsullivan16/asc_auroc:1.0.0dev"

    input:
    tuple val(tool_id), val(meta_id), path(seq_tables)

    output:
    path "tables/combined_res.tsv.gz", emit: combined_res
    path "tables/*.tsv", emit: auc_files, optional: true
    path "tables/purity_stats*.tsv", emit: purity_stats, optional: true
    path "figures/*.png", emit: figs, optional: true

    script:
    """
    get_asc_auroc.R \
    -md ${seq_tables.join(',')} \
    -a ${params.auc_variable} \
    -t $tool_id

    """
}
