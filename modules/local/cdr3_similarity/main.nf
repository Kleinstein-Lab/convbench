process CDR3_SIMILARITY{
    tag "${meta.id}"
    label 'process_medium'

    container "docker.io/cfsullivan16/cdr3similarity:1.0.0dev"

    input:
    tuple val(meta), path(airr), path(embedding)

    output:
    path "tables/run_stats.tsv", emit: run_stats
    path "tables/auc_curve_vals.tsv", emit: auc_vals, optional: true
    path "tables/jaccard_plot_vals.tsv", emit: jaccard_vals, optional: true
    path "tables/seq_summary.tsv", emit: seq_summary
    path "tables/fisher_summary.tsv", emit: fisher_summary
    path "tables/fisher_table.tsv", emit: fisher_table
    path "figures/*.png", emit: figs

    script:
    """
    cdr3_similarity.R \
    -md $airr \
    -o . \
    -da ${params.da_variable} \
    -dg "condition" \
    -t 0.15 \
    -l "single" \
    -a ${params.auc_variable} \
    -v ${params.vdj_info} \
    -sc ${params.single_cell} \
    -r ${params.remove_dups} \
    

    """
}
