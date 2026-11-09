process EVALUATE{
    tag "${meta.id}"
    label 'process_medium'

    container "docker.io/cfsullivan16/evaluate:1.0.0"

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
    def as_list = { files -> [files].flatten() }
    // milo_objects is empty for tools other than Milo
    def milo_arg = as_list(milo_objects) ? "-mo ${as_list(milo_objects).join(',')}" : ""
    """
    evaluate.R \
    -c ${as_list(counts).join(',')} \
    -m ${as_list(memberships).join(',')} \
    -i ${as_list(tool_info).join(',')} \
    ${milo_arg} \
    -li ${library_sizes} \
    -t ${meta.tool_type} \
    -n "${meta.tool_name}" \
    -o . \
    -da ${params.da_variable} \
    -dg ${params.disease_gp} \
    -a ${params.auc_variable}

    """
}
