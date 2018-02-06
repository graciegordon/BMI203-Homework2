from functions import cluster
from functions import io
import os

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

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))
    correct=sorted([0, 1, 1])
    correct2=sorted([0, 0, 1])
    test=cluster.cluster_by_partitioning(active_sites)
    # update this assertion
    assert ((test == correct) or (test == correct2))

def test_hierarchical_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # update this assertion
    assert cluster.cluster_hierarchically(active_sites) == []

