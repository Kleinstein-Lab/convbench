process BCRDIST{
    tag "${meta.id}_${meta.group_id}"
    label 'process_high'

    container "docker.io/cfsullivan16/bcrdist:1.0.0"

    input:
    tuple val(meta), path(airr), path(embedding)

    output:
    tuple val(meta),
          path("tables/*_cluster_counts.tsv"),
          path("tables/*_memberships.tsv"),
          path("tables/*_tool_info.tsv"), emit: eval_input
    path "tables/bcrdist_clusters.tsv", emit: bcrdist_clusters
    path "tables/cluster_summary.tsv", emit: cluster_summary, optional: true
    path "figures/*.png", emit: figs, optional: true

    script:
    """
    # Assuming python in this container is located under /opt/conda/envs/bcrdist/bin/python
    # tcrdist3.py helper function file yet to be integrated
    bcrdist.R \
    -he "bin/tcrdist3.py" \
    -md $airr \
    -g ${meta.group_id} \
    -o . \
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
