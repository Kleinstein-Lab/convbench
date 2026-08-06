process CDR3_SIMILARITY{
    tag "${meta.id}"
    label 'process_medium'

    container "docker.io/cfsullivan16/cdr3similarity:1.0.0dev"

    input:
    tuple val(meta), path(airr), path(embedding)

    output:
    path "tables/run_stats.tsv", emit: run_stats
    path "tables/auc_curve_vals_*.tsv", emit: auc_vals, optional: true
    path "tables/jaccard_plot_vals.tsv", emit: jaccard_vals, optional: true
    path "tables/*_seq_summary.tsv", emit: seq_summary
    path "tables/cluster_subj_summary.tsv", emit: cluster_subj_summary
    path "tables/fisher_table.tsv", emit: fisher_table
    path "tables/cluster_subject_freqs.tsv", emit: cluster_subject_freqs
    path "tables/wilcox_res.tsv", emit: wilcox_res
    path "figures/*.png", emit: figs

    script:
    """
    cdr3_similarity.R \
    -md $airr \
    -o . \
    -da ${params.da_variable} \
    -dg ${params.disease_gp} \
    -t 0.15 \
    -l "single" \
    -a ${params.auc_variable} \
    -v ${params.vdj_info} \
    -sc ${params.single_cell} \
    -r ${params.remove_dups} \
    

    """
}

process CDR3_SIMILARITY_ASC{
    tag "${meta_id}_${asc_id}"
    label 'process_medium'

    container "docker.io/cfsullivan16/cdr3similarity:1.0.0dev"

    input:
    tuple val(meta_id), val(asc_id), path(airr), path(embedding)

    output:
    tuple val(meta_id), path("tables/*_seq_summary.tsv"), emit: auc_input
    path "tables/run_stats.tsv", emit: run_stats
    path "tables/jaccard_plot_vals.tsv", emit: jaccard_vals, optional: true
    path "tables/cluster_subj_summary.tsv", emit: cluster_subj_summary
    path "tables/fisher_table.tsv", emit: fisher_table
    path "tables/cluster_subject_freqs.tsv", emit: cluster_subject_freqs
    path "tables/wilcox_res.tsv", emit: wilcox_res
    path "figures/*.png", emit: figs

    script:
    """
    cdr3_similarity.R \
    -md $airr \
    -o . \
    -da ${params.da_variable} \
    -dg ${params.disease_gp} \
    -t 0.15 \
    -l "single" \
    -a ${params.auc_variable} \
    -v ${params.vdj_info} \
    -sc ${params.single_cell} \
    -r ${params.remove_dups} \

    """
}
