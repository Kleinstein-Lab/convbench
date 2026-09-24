process CDR3_SIMILARITY{
    tag "${meta.id}"
    label 'process_medium'

    container "docker.io/cfsullivan16/cdr3similarity:1.0.0dev"

    input:
    tuple val(meta), path(airr), path(embedding)

    output:
    path "tables/run_stats.tsv", emit: run_stats
    path "tables/evaluation_curve_vals_*.tsv", emit: auc_vals, optional: true
    path "tables/jaccard_plot_vals.tsv", emit: jaccard_vals, optional: true
    path "tables/*_seq_summary.tsv", emit: seq_summary
    path "tables/cluster_subj_summary.tsv", emit: cluster_subj_summary
    path "tables/fisher_table.tsv", emit: fisher_table
    path "tables/cluster_subject_freqs.tsv", emit: cluster_subject_freqs
    path "tables/wilcox_res.tsv", emit: wilcox_res
    path "figures/*.png", emit: figs

    script:
    def args = task.ext.args ? task.ext.args : ""
    """
    cdr3_similarity.R \
    -md ${airr} \
    -o . \
    -da ${params.da_variable} \
    -dg ${params.disease_gp} \
    -t ${meta.threshold} \
    -l ${meta.linkage} \
    -c ${params.cdr3_sim_nproc} \
    -a ${params.auc_variable} \
    -v ${params.vdj_info} \
    -sc ${params.single_cell} \
    -r ${params.remove_dups} \
    ${args}

    """
}

process CDR3_SIMILARITY_ASC{
    tag "${meta.id}_${meta.asc_id}"
    label 'process_medium'

    container "docker.io/cfsullivan16/cdr3similarity:1.0.0dev"

    input:
    tuple val(meta), path(airr), path(embedding)

    output:
    tuple val(meta.id), path("tables/*_seq_summary.tsv"), emit: auc_input
    path "tables/run_stats.tsv", emit: run_stats
    path "tables/jaccard_plot_vals.tsv", emit: jaccard_vals, optional: true
    path "tables/cluster_subj_summary.tsv", emit: cluster_subj_summary
    path "tables/fisher_table.tsv", emit: fisher_table
    path "tables/cluster_subject_freqs.tsv", emit: cluster_subject_freqs
    path "tables/wilcox_res.tsv", emit: wilcox_res
    path "figures/*.png", emit: figs

    script:
    def args = task.ext.args ? task.ext.args : ""
    """
    cdr3_similarity.R \
    -md ${airr} \
    -li ${meta.library_sizes} \
    -o . \
    -da ${params.da_variable} \
    -dg ${params.disease_gp} \
    -t ${meta.threshold} \
    -l ${meta.linkage} \
    -c ${params.cdr3_sim_nproc} \
    -a ${params.auc_variable} \
    -v ${params.vdj_info} \
    -sc ${params.single_cell} \
    -r ${params.remove_dups} \
    ${args}

    """
}
