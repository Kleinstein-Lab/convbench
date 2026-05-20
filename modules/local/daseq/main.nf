process DASEQ{
    tag "${meta.id}"
    label 'process_medium'

    container "docker.io/cfsullivan16/daseq:1.0.0dev"

    input:
    tuple val(meta), path(airr), path(embedding)

    output:
    path "tables/run_stats.tsv", emit: run_stats
    path "tables/da_score_auc_curve_vals.tsv", emit: da_score_auc_vals, optional: true
    path "tables/p_auc_curve_vals.tsv", emit: p_auc_vals, optional: true
    path "tables/jaccard_plot_vals.tsv", emit: jaccard_vals, optional: true
    path "tables/da_seqs.tsv", emit: da_seqs
    path "tables/region_stats.tsv", emit: region_stats
    path "tables/cluster_cts.tsv", emit: cluster_cts
    path "tables/UMAP_embeddings.rds", emit: daseq_umap
    path "tables/da_cells.rds", emit: da_cells, optional: true
    path "figures/*.png", emit: figs

    script:
    """
    daseq.R \
    -d $embedding \
    -md $airr \
    -o . \
    -da "status" \
    -m 10 \
    -t 25 \
    -x 250 \
    -v TRUE \
    -sc FALSE \
    -si TRUE \
    -r FALSE \
    -w TRUE

    """
}
