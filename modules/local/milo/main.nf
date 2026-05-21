process MILO{
    tag "${meta.id}"
    label 'process_medium'

    container "docker.io/cfsullivan16/milo:1.0.0dev"

    input:
    tuple val(meta), path(airr), path(embedding)

    output:
    path "tables/run_stats.tsv", emit: run_stats
    path "tables/auc_curve_vals.tsv", emit: auc_vals, optional: true
    path "tables/jaccard_plot_vals.tsv", emit: jaccard_vals, optional: true
    path "tables/seq_results.tsv", emit: seq_results
    path "tables/da_results.tsv", emit: da_results
    path "tables/nhood_stats.tsv", emit: nhood_stats
    path "tables/milo_nhoods.RDS", emit: milo_nhoods, optional: true
    path "figures/*.png", emit: figs

    script:
    """
    milo.R \
    -d $embedding \
    -md $airr \
    -o . \
    -da "status" \
    -k 20 \
    -pr 0.1 \
    -v TRUE \
    -sc FALSE \
    -si TRUE \
    -r FALSE \
    -w TRUE

    """
}
