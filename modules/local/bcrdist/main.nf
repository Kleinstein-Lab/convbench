process BCRDIST{
    tag "${meta.id}"
    label 'process_high'

    container "docker.io/cfsullivan16/bcrdist:1.0.0dev"

    input:
    tuple val(meta), path(airr), path(embedding)

    output:
    path "tables/run_stats.tsv", emit: run_stats
    path "tables/auc_curve_vals_*.tsv", emit: auc_vals, optional: true
    path "tables/jaccard_plot_vals.tsv", emit: jaccard_vals, optional: true
    path "tables/bcrdist_clusters.tsv", emit: bcrdist_clusters
    path "tables/cluster_subject_freqs.tsv", emit: cluster_subject_freqs
    path "tables/cluster_summary.tsv", emit: cluster_summary, optional: true
    path "tables/*_seq_summary.tsv", emit: seq_summary
    path "tables/cluster_subj_summary.tsv", emit: cluster_subj_summary
    path "tables/fisher_table.tsv", emit: fisher_table
    path "tables/wilcox_res.tsv", emit: wilcox_res
    path "figures/*.png", emit: figs

    script:
    """
    # Assuming python in this container is located under /opt/conda/envs/bcrdist/bin/python
    # tcrdist3.py helper function file yet to be integrated

    bcrdist.R \
    -he "bin/tcrdist3.py" \
    -md $airr \
    -li ${params.lib_size_override} \
    -o . \
    -da ${params.da_variable} \
    -dg ${params.disease_gp} \
    -t 60 \
    -a ${params.auc_variable} \
    -py "/opt/conda/envs/bcrdist/bin/python" \
    -i FALSE \
    -c ${task.cpus} \
    -cs 1000 \
    -r 60 \
    -m "greedy"

    """
}

process BCRDIST_ASC{
    tag "${meta_id}_${asc_id}"
    label 'process_high'

    container "docker.io/cfsullivan16/bcrdist:1.0.0dev"

    input:
    tuple val(meta_id), val(asc_id), path(airr), path(embedding), path(library_sizes)

    output:
    tuple val(meta_id), path("tables/*_seq_summary.tsv"), emit: auc_input
    path "tables/run_stats.tsv", emit: run_stats
    path "tables/auc_curve_vals_*.tsv", emit: auc_vals, optional: true
    path "tables/jaccard_plot_vals.tsv", emit: jaccard_vals, optional: true
    path "tables/bcrdist_clusters.tsv", emit: bcrdist_clusters
    path "tables/cluster_subject_freqs.tsv", emit: cluster_subject_freqs
    path "tables/cluster_summary.tsv", emit: cluster_summary, optional: true
    path "tables/cluster_subj_summary.tsv", emit: cluster_subj_summary
    path "tables/fisher_table.tsv", emit: fisher_table
    path "tables/wilcox_res.tsv", emit: wilcox_res
    path "figures/*.png", emit: figs

    script:
    """
    # Assuming python in this container is located under /opt/conda/envs/bcrdist/bin/python
    # tcrdist3.py helper function file yet to be integrated

    bcrdist.R \
    -he "bin/tcrdist3.py" \
    -md $airr \
    -li $library_sizes \
    -o . \
    -da ${params.da_variable} \
    -dg ${params.disease_gp} \
    -t 60 \
    -a ${params.auc_variable} \
    -py "/opt/conda/envs/bcrdist/bin/python" \
    -i FALSE \
    -c ${task.cpus} \
    -cs 1000 \
    -r 60 \
    -m "greedy"

    """
}
