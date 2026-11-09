process CDR3_SIMILARITY{
    tag "${meta.id}_${meta.group_id}"
    label 'process_medium'

    container "docker.io/cfsullivan16/cdr3similarity:1.0.0"

    input:
    tuple val(meta), path(airr), path(embedding)

    output:
    tuple val(meta),
          path("tables/*_cluster_counts.tsv"),
          path("tables/*_memberships.tsv"),
          path("tables/*_tool_info.tsv"), emit: eval_input

    script:
    def args = task.ext.args ? task.ext.args : ""
    """
    cdr3_similarity.R \
    -md ${airr} \
    -g ${meta.group_id} \
    -o . \
    -t ${meta.threshold} \
    -l ${meta.linkage} \
    -c ${task.cpus} \
    -a ${params.auc_variable} \
    -v ${params.vdj_info} \
    -sc ${params.single_cell} \
    -r ${params.remove_dups} \
    ${args}

    """
}
