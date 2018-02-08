from functions import cluster
from functions import io
import os
import numpy as np

def test_similarity():
    filename_a = os.path.join("data", "276.pdb")
    filename_b = os.path.join("data", "4629.pdb")

    #I've written the similarity function to accept floats to make the rest of my code flow easier
    activesite_a = io.read_active_site(filename_a)
    activesite_b = io.read_active_site(filename_b)

    # update this assertion
    assert round(cluster.compute_similarity(activesite_a, activesite_b),2) == float(7.51)

def test_partition_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701]
    #pdb_ids=[10701,276,4629]
    #pdb_ids = [276, 4629, 10701,34047,29047] 
    
    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))
    test=cluster.cluster_by_partitioning(active_sites)
    #correct=np.array([[276, 4629], [10701]])
    #correct2=np.array([[10701], [276, 4629]])
    correct=[[active_sites[0], active_sites[1]], [active_sites[2]]]
    correct2=[[active_sites[2]], [active_sites[0], active_sites[1]]]
    # update this assertion
    assert ((np.array_equal(test,correct)) or (np.array_equal(test,correct2)))
    #assert test == []

def test_hierarchical_clustering():
    # tractable subset
    #pdb_ids = [276, 4629, 10701,34047,29047]
    pdb_ids=[276, 4629, 10701]
    #pdb_ids=[10701,276,4629]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))
    #print(cluster.cluster_hierarchically(active_sites))
    correct=[[active_sites[0], active_sites[1]], [active_sites[2]]]
    correct2=[[active_sites[2]], [active_sites[0], active_sites[1]]]
    test=(cluster.cluster_hierarchically(active_sites))
    # update this assertion
    print(test,correct)
    assert ((np.array_equal(test,correct)) or (np.array_equal(test,correct2)))
    #assert test == []

def test_clust_quality():
    pdb_ids=[276, 4629, 10701]
    #pdb_ids = [276, 4629, 10701,34047,29047]
    
    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))
    clust=[[active_sites[0], active_sites[1]], [active_sites[2]]]
    #clust=[[active_sites[2],active_sites[4]],[active_sites[0],active_sites[1],active_sites[3]]]
    test=(cluster.check_clust_quality(clust))
    assert round(test,2) == float(1.88)

def test_compare_clusters():
    #write test for comparison
    pdb_ids=[276, 4629, 10701]
    #pdb_ids = [276, 4629, 10701,34047,29047]
    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))
    
    clustA=cluster.cluster_hierarchically(active_sites)
    clustB=cluster.cluster_by_partitioning(active_sites)
    comp=cluster.compare_clusters(clustA,clustB,active_sites)
    comp = 1
    assert comp == []
