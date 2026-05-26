process BCRDIST{
    tag "${meta.id}"
    label 'process_medium'

    container "docker.io/cfsullivan16/bcrdist:1.0.0dev"

    input:
    tuple val(meta), path(airr), path(embedding)

    output:
    path "tables/run_stats.tsv", emit: run_stats
    path "tables/auc_curve_vals.tsv", emit: auc_vals, optional: true
    path "tables/jaccard_plot_vals.tsv", emit: jaccard_vals, optional: true
    path "tables/bcrdist_clusters.tsv", emit: bcrdist_clusters
    path "tables/cluster_summary.tsv", emit: cluster_summary, optional: true
    path "tables/seq_summary.tsv", emit: seq_summary
    path "tables/fisher_summary.tsv", emit: fisher_summary
    path "tables/fisher_table.tsv", emit: fisher_table
    path "figures/*.png", emit: figs

    script:
    """
    # Assuming python in this container is located under /opt/conda/envs/bcrdist/bin/python
    # tcrdist3.py helper function file yet to be integrated

    bcrdist.R \
    -he "bin/tcrdist3.py" \
    -md $airr \
    -o . \
    -da "status" \
    -dg "condition" \
    -t 60 \
    -si TRUE \
    -py "/opt/conda/envs/bcrdist/bin/python" \
    -i FALSE \
    -c 4 \
    -cs 1000 \
    -r 60 \
    -m "greedy"

    """
}
