process MILO{
    tag "${meta.id}_${meta.group_id}"
    label 'process_medium'

    container "docker.io/cfsullivan16/milo:1.0.0"

    input:
    tuple val(meta), path(airr), path(embedding)

    output:
    tuple val(meta),
          path("tables/*_cluster_counts.tsv"),
          path("tables/*_memberships.tsv"),
          path("tables/*_tool_info.tsv"),
          path("tables/*_milo.rds"), emit: eval_input
    path "tables/nhood_stats.tsv", emit: nhood_stats
    path "figures/*.png", emit: figs

    script:
    """
    milo.R \
    -d $embedding \
    -md $airr \
    -g ${meta.group_id} \
    -o . \
    -da ${params.da_variable} \
    -k ${params.milo_k_val} \
    -pr 0.1 \
    -a ${params.auc_variable} \
    -v ${params.vdj_info} \
    -sc ${params.single_cell} \
    -r ${params.remove_dups}

    """
}
