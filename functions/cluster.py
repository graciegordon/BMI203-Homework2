#This script will implement K-means
#from .utils import Atom, Residue, ActiveSite
from functions import utils
import math
import numpy as np
import random

def compute_hydrophobicity_Index(site):
    #calculate Hydrophobicity Index for aa. Source: sigma-aldrich
    hydrophobicityIndexNeutral={'PHE':100,'ILE':99,'TRP':97,'LEU':97,'VAL':76, 'MET':74,'TYR':63,'CYS':49,'ALA':41,'THR':13,'HIS':8,'GLY':0,'SER':-5,'GLN':-10,'ARG':-14, 'LYS':-23,'ASN':-28,'GLU':-31,'PRO':-46,'ASP':-55}
    
    H_index=0
    numRes=0
    for res in site.residues:
        sl = slice(0,3)
        tempaa=str(res)
        aa=tempaa[sl]
        H_index+=float(hydrophobicityIndexNeutral[aa])
        numRes+=1
    H_index_final=H_index/numRes
    return H_index_final

def compute_similarity(site_a, site_b):
    """
    Compute the similarity between two given ActiveSite instances.
    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    """

    similarity = 0.0

    #similarity distance computed with Euclidean Distance
    siteA_idx=compute_hydrophobicity_Index(site_a)
    siteB_idx=compute_hydrophobicity_Index(site_b)
    similarity=math.sqrt((siteA_idx-siteB_idx)**2)

    return similarity


def cluster_by_partitioning(active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.
    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    
    #Use K-means Partitioning Algorithm
    #set number of clusters must be less than number of clusters
    numClusters=2
    numSites=0
    for x in active_sites:
        numSites+=1
    
    #create matrix with feature values
    features=[]
    for z in active_sites:
        features.append(compute_hydrophobicity_Index(z))

    #initialize centroids
    centroid=[]
    
    for i in range(numClusters):
        #tempCentroid=random.choice(active_sites)
        #print('pos before',pos)
        if i ==0:
            pos=np.random.choice(features, 1)
        while pos in centroid:
            pos=np.random.choice(features, 1)
            print('pos',pos)
        pos=float(str(pos[0]))
        centroid.append(pos)
    Cidx=[]
    
    #define number of iterations
    for iteration in range(0,1):
        #assign clusters
        #list of indicies of closest cluster for each active site
        idxtemp=[0]*numSites
        
        #for each active site, assign it a cluster based on distance
        for j in range(0,len(idxtemp)):
            minimum=float("inf")
            #site=features[j]
            site=active_sites[j]
            for i in range(0,len(centroid)):
                #this is a stupid work around do better to get back to active sites
                for item in active_sites:
                    if centroid[i]==(compute_hydrophobicity_Index(item)):
                            tempActive=item
                            test=compute_similarity(site, tempActive)
                            if test <= minimum:
                                minimum=test
                                idxtemp[j]=i
        '''   
        #test if clustering is the same, if it is return 
        if idxtemp == Cidx:
            C=[]
            for item in idxtemp:
                C.append(centroid[item])
            print('c', C)
            print('idx',idxtemp)
            return sorted(idxtemp)
            break
        elif idxtemp != Cidx:
            Cidx=idxtemp
        '''

        #recalculate centers for each cluster
        for j in range(0,len(centroid)):
            nums=idxtemp.count(j)
            vect=[0]*nums
            counter=0
            #collect all clusters of the same and get the average point
            for i in range(0,len(Cidx)):
                print(Cidx[i],j)
                if (Cidx[i]==j):
                    vect[counter]=float(features[i])
                    counter+=1
            if vect==[]:
                vect.append(0)
            avg=np.mean(vect)
            centroid[j]=float(avg)
        
        #test if clustering is the same, if it is return 
        if idxtemp == Cidx:
            C=[]
            for item in idxtemp:
                C.append(centroid[item])
            print('c', C)
            print('idx',idxtemp)
            return sorted(idxtemp)
            break
        elif idxtemp != Cidx:
            Cidx=idxtemp

    print('idx',idxtemp)
    print('idx',sorted(idxtemp))
    return sorted(idxtemp)
    
    
def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #
    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """
    print('aggs')
    #implement Agglomerative Algorithm
    listoflists=[]
    #determine how many clusters to stop at
    numberclust=2
    print('active site', active_sites) 
    #put each object in its own cluster
    for site in active_sites:
        print(site)
        print(site.residues)
        a_list=[]
        a_list.append(site)
        print(a_list)
        listoflists.append(a_list)
        #a_list=[]
    print('listsolists',listoflists)
    #iterate through list compare distance and merge cells with closest distance
    #do until number of clusters desired is reached
    #for i in range(4):
    while len(listoflists)!=numberclust:
        #find minimum distance between the clusters
        #get first to compare
        bestclusterstomerge=float("inf")
        for clust in listoflists:
            #get second cluster to compare
            for compare in listoflists:
                #lowestdistforKNN=float("inf")
                if compare != clust:
                    #compare each element in list to find the minimum distance between the cluster (KNN)
                    lowestdistforKNN=float("inf")
                    for i in clust:
                        for j in compare:
                            disttemp=compute_similarity(i,j)
                            #save lowest distance found between two clusters
                            if disttemp<lowestdistforKNN:
                                lowestdistforKNN=disttemp
                    #save lowest distance for these clusters
                    lowestdistfortwoclusters=lowestdistforKNN
                    print('low cluster',lowestdistfortwoclusters)
                    print('last best',bestclusterstomerge)
                    #if this is the lowest distance found so far, save this info
                    if lowestdistfortwoclusters<bestclusterstomerge:
                        bestclusterstomerge=lowestdistfortwoclusters
                        clusterA=clust
                        clusterB=compare
                    print('best this clust so far',lowestdistfortwoclusters)
                    print('best of all', bestclusterstomerge,clusterA,clusterB)
        #merge two closest clusters
        print('to merge',clusterA,clusterB,bestclusterstomerge)
        idxlist=[]
        idxlist.append(listoflists.index(clusterA))
        idxlist.append(listoflists.index(clusterB))
        idxlist=sorted(idxlist)
    
        mergeitems=listoflists[idxlist[0]]+listoflists[idxlist[1]]
        print(mergeitems)
        listoflists.remove(listoflists[idxlist[0]])
        idxlist[1]=idxlist[1]-1
        listoflists.remove(listoflists[idxlist[1]])
        listoflists=listoflists+[mergeitems]
        print(listoflists) 
    print('final',listoflists)
    return listoflists
