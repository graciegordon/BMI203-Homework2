import sys
from functions import cluster
from functions import io
from functions import utils
import matplotlib.pyplot as plt

#from .io import read_active_sites, write_clustering, write_mult_clusterings
#from .cluster import cluster_by_partitioning, cluster_hierarchically
#def plotting_func(clusterA,clusterB):
def plot_clustser(clust,title):

    col=['red','blue','black']
    iter=0
    for i in clust:
        temp=[]
        #print(i)
        for x in i:
            hydro=cluster.compute_hydrophobicity_Index(x)
            temp.append(hydro)
        #print('convert active',x)
        plt.scatter(temp, temp, color=col[iter])
        iter+=1
    plt.title(title)
    plt.show()

# Some quick stuff to make sure the program is called correctly
if len(sys.argv) < 4:
    print("Usage: python -m hw2skeleton [-P| -H] <pdb directory> <output file>")
    sys.exit(0)

active_sites = io.read_active_sites(sys.argv[2])

# Choose clustering algorithm
if sys.argv[1][0:2] == '-P':
    print("Clustering using Partitioning method")
    clustering = cluster.cluster_by_partitioning(active_sites)
    io.write_clustering(sys.argv[3], clustering)

if sys.argv[1][0:2] == '-H':
    print("Clustering using hierarchical method")
    clusterings = cluster.cluster_hierarchically(active_sites)
    io.write_clustering(sys.argv[3], clusterings)
    #io.write_mult_clusterings(sys.argv[3], clusterings)

if sys.argv[1][0:2] == '-A':
    clusteringPart = cluster.cluster_by_partitioning(active_sites)
    clusteringHi = cluster.cluster_hierarchically(active_sites)
    
    QualPart=cluster.check_clust_quality(clusteringPart)
    QualHi=cluster.check_clust_quality(clusteringHi)
    
    RandIndex=cluster.compare_clusters(clusteringPart,clusteringHi,active_sites)

    print("Quatity Partition, Agg :", QualPart, QualHi)
    print("RandIndex",RandIndex)
   

    plot_clustser(clusteringPart, 'partition clustering')
    plot_clustser(clusteringHi, 'hierarchical clustering')
