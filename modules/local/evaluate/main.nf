process EVALUATE{
    tag "${meta.tool_dir}_${meta.id}"
    label 'process_medium'

    container "docker.io/cfsullivan16/evaluate:1.1.0"

    input:
    tuple val(meta), path(counts), path(memberships), path(tool_info), path(milo_objects), path(library_sizes)

    output:
    path "tables/run_stats.tsv", emit: run_stats
    path "tables/cluster_test_results.tsv", emit: cluster_test_results
    path "tables/sequence_scores.tsv", emit: sequence_scores
    path "tables/cluster_subject_freqs.tsv", emit: cluster_subject_freqs
    path "tables/evaluation_curve_vals_*.tsv", emit: auc_vals, optional: true
    path "tables/purity_stats.tsv", emit: purity_stats, optional: true
    path "figures/*.png", emit: figs, optional: true

    script:
    // a single file is staged as a path rather than a list
    def counts_list = [counts].flatten()
    def memberships_list = [memberships].flatten()
    def tool_info_list = [tool_info].flatten()
    // milo_objects is empty for tools other than Milo
    def milo_list = [milo_objects].flatten()
    def milo_arg = milo_list ? "-mo ${milo_list.join(',')}" : ""
    """
    evaluate.R \
    -c ${counts_list.join(',')} \
    -m ${memberships_list.join(',')} \
    -i ${tool_info_list.join(',')} \
    ${milo_arg} \
    -li ${library_sizes} \
    -t ${meta.tool_type} \
    -n "${meta.tool_name}" \
    -o . \
    -da ${params.da_variable} \
    -dg ${params.disease_gp} \
    -cg ${params.control_gp} \
    -a ${params.auc_variable}

    """
}
