#!/usr/bin/env python

import pandas as pd
import numpy as np
from tcrdist.repertoire import TCRrep
import pwseqdist as pw
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

def run_tcrdist3(seqs,id_col,
                 cpus = 4,chunk_size=1000,radius=200,
                 #save_dist_mtx=False,save_path=None,
                 v2cdr = True):

    # Standardize the parameters as integer
    cpus = int(cpus)
    chunk_size = int(chunk_size)
    radius = int(radius)
    
    tr = TCRrep(
    cell_df = seqs,
    organism = 'human', 
    chains = ['beta'],
    db_file = 'combo_xcr_2024-03-05.tsv',
    compute_distances = False)

    if v2cdr == False:
        required_cols = ["pmhc_b_aa", "cdr1_b_aa", "cdr2_b_aa"]
        # Check if all required columns exist in the original dataframe
        if all(col in seqs.columns for col in required_cols):
            tr.clone_df = seqs
        else:
            # Identify exactly which ones are missing
            missing_col = [col for col in required_cols if col not in df.columns]
            print(f"Error: columns are missing: {missing_col} in the input dataframe. Inferred CDRs are used.")
            
    # Specifify the parameters
    
    # metrics_a = {
    #    "cdr3_a_aa" : pw.metrics.nb_vector_tcrdist,
    #    "pmhc_a_aa" : pw.metrics.nb_vector_tcrdist,
    #    "cdr2_a_aa" : pw.metrics.nb_vector_tcrdist,
    #    "cdr1_a_aa" : pw.metrics.nb_vector_tcrdist}
    
    metrics_b = {
        "cdr3_b_aa" : pw.metrics.nb_vector_tcrdist,
        "pmhc_b_aa" : pw.metrics.nb_vector_tcrdist,
        "cdr2_b_aa" : pw.metrics.nb_vector_tcrdist,
        "cdr1_b_aa" : pw.metrics.nb_vector_tcrdist }
    
    # weights_a= { 
    #    "cdr3_a_aa" : 3,
    #    "pmhc_a_aa" : 1,
    #    "cdr2_a_aa" : 1,
    #    "cdr1_a_aa" : 1}
    
    weights_b = { 
        "cdr3_b_aa" : 3,
        "pmhc_b_aa" : 0,
        "cdr2_b_aa" : 1,
        "cdr1_b_aa" : 1}
    
    # kargs_a = {  
    #    'cdr3_a_aa' : 
    #        {'use_numba': True, 
    #        'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 
    #        'dist_weight': 1, 
    #        'gap_penalty':4, 
    #        'ntrim':3, 
    #        'ctrim':2, 
    #        'fixed_gappos': False},
    #    'pmhc_a_aa' : {
    #        'use_numba': True,
    #        'distance_matrix': pw.matrices.tcr_nb_distance_matrix,
    #        'dist_weight':1,
    #        'gap_penalty':4,
    #        'ntrim':0,
    #        'ctrim':0,
    #        'fixed_gappos':True},
    #    'cdr2_a_aa' : {
    #        'use_numba': True,
    #        'distance_matrix': pw.matrices.tcr_nb_distance_matrix,
    #        'dist_weight': 1,
    #        'gap_penalty':4,
    #        'ntrim':0,
    #        'ctrim':0,
    #        'fixed_gappos':True},
    #    'cdr1_a_aa' : {
    #        'use_numba': True,
    #        'distance_matrix': pw.matrices.tcr_nb_distance_matrix,
    #        'dist_weight':1,
    #        'gap_penalty':4,
    #        'ntrim':0,
    #        'ctrim':0,
    #        'fixed_gappos':True}
    #    }
        
    kargs_b= {  
        'cdr3_b_aa' : 
            {'use_numba': True, 
            'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 
            'dist_weight': 1, 
            'gap_penalty':4, 
            'ntrim':3, 
            'ctrim':2, 
            'fixed_gappos': False},
        'pmhc_b_aa' : {
            'use_numba': True,
            'distance_matrix': pw.matrices.tcr_nb_distance_matrix,
            'dist_weight': 1,
            'gap_penalty':4,
            'ntrim':0,
            'ctrim':0,
            'fixed_gappos': True},
        'cdr2_b_aa' : {
            'use_numba': True,
            'distance_matrix': pw.matrices.tcr_nb_distance_matrix,
            'dist_weight':1,
            'gap_penalty':4,
            'ntrim':0,
            'ctrim':0,
            'fixed_gappos': True},
        'cdr1_b_aa' : {
            'use_numba': True,
            'distance_matrix': pw.matrices.tcr_nb_distance_matrix,
            'dist_weight':1,
            'gap_penalty':4,
            'ntrim':0,
            'ctrim':0,
            'fixed_gappos': True}
        }   

    # tr.metrics_a = metrics_a
    tr.metrics_b = metrics_b

    # tr.weights_a = weights_a
    tr.weights_b = weights_b

    # tr.kargs_a = kargs_a 
    tr.kargs_b = kargs_b

    # Assign computer resources
    tr.cpus = cpus
    tr.compute_sparse_rect_distances(radius = radius, chunk_size = chunk_size)

    # if save_dist_mtx:
      ### !!!!!! Make sure seq_id is always added !!!!!! ###
    #  seq_ids = tr.clone_df['seq_id'] 
    #  if save_path is not None:
    #    print('Saving the pairwise distance matrix as parquet...')
    #    dis_mtx = pd.DataFrame(tr.rw_beta.toarray(),index=seq_ids,columns=seq_ids)
    #    dis_mtx.to_parquet(save_path, engine='pyarrow', compression='snappy')
    #  else:
    #    print('Not saving because save_path was not provided')
    
    return tr



def greedy_clustering(mtx,threshold=120,recover0 = False):
    
    ### Calculate counts of neighbors within threshold (default 120) ###
    
    # Remove edges > threshold
    mtx.data[mtx.data > threshold] = 0
    mtx.eliminate_zeros()
    
    # Recover distance = -1 to 0 (real one)
    if recover0:
        mtx.data[mtx.data == -1] = 0
    
    # In a CSR matrix, the difference between consecutive indptr values suggests
    # how many non-zero elements are in that row (neighbors within threshold)
    degrees = np.diff(mtx.indptr)
    
    # Rank the seqs by neighbor counts
    nodes = np.argsort(degrees)[::-1]
    
    # Flag indicating if the node/seq has been clustered already
    # False means clustered and left out from the following search
    n_seqs = mtx.shape[0]
    flag = np.ones(n_seqs, dtype=bool)
    
    clusters = [] # Collect the clusters

    # Iterate the nodes/seqs as the clustering center
    for node in nodes:
        # Stop searching if it has no neighbors
        if degrees[node] <= 1:
            break
        
        # Skip it if it's clustered
        if not flag[node]:
            continue
            
        # 1. Find its neighbors
        start_idx = mtx.indptr[node]
        end_idx = mtx.indptr[node+1]
        neighbors = mtx.indices[start_idx:end_idx]
        
        # 2. Filter out neighbors that have already been clustered
        neighbors = neighbors[flag[neighbors]]
        
        # 3. Form the new cluster (using set to avoid self-loop)
        cluster = list(set([node] + neighbors.tolist()))
        if len(cluster) > 1:
            clusters.append(cluster)
        
        # 4. Marking the nodes/seqs clustered
        flag[cluster] = False

    n_clustered = np.sum(flag == False) # how many seqs clustered
        
    print(f"{len(clusters)} clusters found")
    print(f"{n_clustered} (out of {n_seqs}) sequences being clustered")
    return clusters



# Convert greedy clustering outputs (a list of clusters)
# to an array of cluster labels
def clusters_list2array(clusters, n_seqs):

    # 0 means not belonging to any clusters
    cluster_labels = np.zeros(n_seqs, dtype=int)
    
    # Assign the cluster labels from ID=1
    for cluster_id, cluster_nodes in enumerate(clusters, start=1):
        cluster_labels[cluster_nodes] = cluster_id
        
    # Find the nodes still labelled 0/unclustered
    no_clusters = (cluster_labels == 0)
    n_no_clusters = np.sum(no_clusters)
    
    # Assign unique IDs to the unclustered nodes
    if n_no_clusters > 0:
        start_new_id = len(clusters) + 1
        cluster_labels[no_clusters] = np.arange(
            start_new_id, 
            start_new_id + n_no_clusters
        )
        
    return cluster_labels



def single_linkage_clustering(mtx,threshold=120,recover0 = False):
    
    ### Calculate counts of neighbors within threshold (default 120) ###
    
    # Remove edges > threshold
    mtx.data[mtx.data > threshold] = 0
    mtx.eliminate_zeros()
    
    # Recover distance = -1 to 0 (real one)
    if recover0:
        mtx.data[mtx.data == -1] = 0

    # Nodes with any edges form clusters
    _, cluster_labels = connected_components(
        csgraph=mtx, 
        directed=False, 
        return_labels=True
    )

    # To match index in R
    cluster_labels = cluster_labels + 1

    # Clusters having at least 2 sequences
    _, label_counts = np.unique(cluster_labels, return_counts=True)
    n_clusters = np.sum(label_counts > 1)
    n_clustered = np.sum(label_counts[label_counts > 1])
    
    print(f"{n_clusters} clusters found")
    print(f"{n_clustered} (out of {mtx.shape[0]}) sequences being clustered")
    return cluster_labels


    

def tcrdist3_clusters(seqs,id_col,
                     cpus=4,chunk_size=1000,radius=150,
                     v2cdr=True,
                     mode='greedy',threshold=120,recover0=False,
                     save_dist_mtx=False,save_path=None):
    # Mode choices
    modes = ['greedy','single_linkage']
    if mode not in modes:
        raise ValueError(f"Invalid mode! Must be one of {modes}")
    print(f"Clustering in {mode} fashion!")

    # Standardize the parameters as integer
    cpus = int(cpus)
    chunk_size = int(chunk_size)
    radius = int(radius)
    threshold = int(threshold)
    
    print(f"Running tcrdist3 (radius = {radius})")
    tr = run_tcrdist3(seqs=seqs, id_col=id_col,
                      cpus=cpus, chunk_size=chunk_size, radius=radius,
                      #save_dist_mtx=save_dist_mtx, save_path=save_path,
                      v2cdr=v2cdr)
    
    dist_mtx = tr.rw_beta
    seqs_new = tr.clone_df
    
    if mode == 'greedy':
        print(f"Finding clusters (threshold = {threshold})")
        clusters_list = greedy_clustering(mtx = dist_mtx,threshold = threshold,recover0=recover0)
        clusters = clusters_list2array(clusters_list,n_seqs=dist_mtx.shape[0])
    elif mode == 'single_linkage':
        print(f"Finding clusters (threshold = {threshold})")
        clusters = single_linkage_clustering(dist_mtx,threshold = threshold,recover0=recover0)

    seqs_new['convergent_clone_id'] = clusters.astype(str)
    seqs_new['clone_id_full'] = seqs_new['convergent_clone_id'].astype(str) + '_' + seqs_new['subject_id'].astype(str)
    
    return(seqs_new)
