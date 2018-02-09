#This script will implement K-means
#from .utils import Atom, Residue, ActiveSite
from functions import utils
import math
import numpy as np
import random
def compute_hydrophobicity_Index(site):
    #calculate Hydrophobicity Index for aa. Source: sigma-aldrich
    hydrophobicityIndexNeutral={'PHE':100,'ILE':99,'TRP':97,'LEU':97,'VAL':76, 'MET':74,'TYR':63,'CYS':49,'ALA':41,'THR':13,'HIS':8,'GLY':0,'SER':-5,'GLN':-10,'ARG':-14, 'LYS':-23,'ASN':-28,'GLU':-31,'PRO':-46,'ASP':-55}
    
    #for each residue in the active site, sum hydrophobicity score and average
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
    ###################################
    #Use K-means Partitioning Algorithm 
    ###################################

    #set number of clusters must be less than number of clusters
    numClusters=2
    numSites=0
    for x in active_sites:
        numSites+=1    
    #create matrix with feature values
    features=[]
    for z in active_sites:
        features.append(compute_hydrophobicity_Index(z))
   
    Cidx=[]
    #initialize centroids
    centroid=[]
    #make sure clusters are not the same
    centroid=np.random.choice(features,numClusters,replace=False)
    while len(set(centroid)) != len(centroid):
        centroid=np.random.choice(features,numClusters,replace=False)
    
    #define number of iterations
    for iteration in range(0,50):
        
        ################
        #assign clusters
        ################
        
        #list of centroids of closest cluster for each active site (matched index for feature,active and idx)
        idxtemp=[0]*numSites
        for num in range(len(idxtemp)):
            mindist=float('inf')
            #calculate distance to each centroid and assign to lowest value
            for cent in centroid:
                dist=math.sqrt((cent-features[num])**2)
                if dist<mindist:
                    mindist=dist
                    idxtemp[num]=cent
        
        ###################################################################
        #average elements assigned to each vector and recalculate centroids
        ###################################################################

        tempclust=[]
        centid=0
        for item in centroid:
            #get number of items assigned to each cluster and make vector
            nums=idxtemp.count(item)
            avgvect=[0]*nums
            clustvect=[0]*nums
            count=0
            i=0
            #collect all values for the cluster
            for j in (idxtemp):
                if j == item:
                    #append distance to recalculate centroid
                    #append active site numbers to show these clusters
                    avgvect[count]=features[i]
                    clustvect[count]=active_sites[i]
                    count+=1
                i+=1
            #calculate new centroid location
            newcent=np.mean(avgvect)
            centroid[centid]=newcent
            centid+=1
            #create list of list for current clusters
            tempclust.append(clustvect)

        #####################################################
        #determine if clusters stayed the same this iteration
        #if so return current clsuters, else keep iterating
        #####################################################

        if tempclust == Cidx:
            clustfinal=tempclust
            #print('final found',clustfinal)
            return clustfinal
        else:
            Cidx=tempclust
    
    #max iter reached
    clustfinal=tempclust
    print('final reached',clustfinal)
    return clustfinal
    
    
def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #
    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """
    
    ##################################
    #implement Agglomerative Algorithm
    ##################################

    listoflists=[]
    #determine how many clusters to stop at
    numberclust=2
    
    #put each object in its own cluster
    for site in active_sites:
        a_list=[]
        a_list.append(site)
        listoflists.append(a_list)
   
    ############################################################################
    #iterate through list compare distance by nearest neighbor metric
    #and merge cells with closest distance
    #do until number of clusters desired is reached
    ###########################################################################

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
                    #if this is the lowest distance found so far, save this info
                    if lowestdistfortwoclusters<bestclusterstomerge:
                        bestclusterstomerge=lowestdistfortwoclusters
                        clusterA=clust
                        clusterB=compare
        #merge two closest clusters
        #determine the indicies to merge
        idxlist=[]
        idxlist.append(listoflists.index(clusterA))
        idxlist.append(listoflists.index(clusterB))
        idxlist=sorted(idxlist)
    
        #merge the indicies given
        #for ind in idxlist:
        mergeitems=listoflists[idxlist[0]]+listoflists[idxlist[1]]
        listoflists.remove(listoflists[idxlist[0]])
        idxlist[1]=idxlist[1]-1
        listoflists.remove(listoflists[idxlist[1]])
        listoflists=listoflists+[mergeitems]
    #print('final',listoflists)
    return listoflists

def check_clust_quality(clusterlist):
    #Use Sum of Distances Method to check cluster
    totclust=0
    for clust in clusterlist:
        clustdist=0
        #sum and normalize the distances between elements in the cluster
        for i in range(len(clust)):
            for j in range(i+1,len(clust)):            
                disttemp=compute_similarity(clust[i],clust[j])
                clustdist+=disttemp
        clustdist=clustdist/len(clust)
        #sum cluster total 
        totclust+=clustdist 
    #average the total sum of distances for all clusters by number of clusters
    final=totclust/len(clusterlist)
    #return the total cluster average
    return final

def convert_sites_to_clust_idx(clust,active_sites):
    #convert sites to 1/2 to assign clusters
    idx=[]
    for num in active_sites:
        for i in range(len(clust)):
            if num in clust[i]:
                idx.append(i+1)
            
    return idx

def guarentee_site_order(cluster):
    pair={}
    originalidx=0
    #sort clusters by lowest hydrophobicity index to guarentee they are ordert the same between groups
    for item in cluster:
        tot=0
        for element in item:
            num=compute_hydrophobicity_Index(element)
            tot+=num
        tot=tot/len(item)
        pair[tot]=originalidx
        originalidx+=1

    clustnew=[]
    for key in sorted(pair):
        temp=pair[key]
        clustnew.append(cluster[temp])
    
    return clustnew

def compare_clusters(clusterA,clusterB,active_sites):
    #clusterA will be lower hydrophobicity and B will be higher 
   
    #sort sites so clusters are ordered by average hydorphobicity
    clusterA=guarentee_site_order(clusterA)
    clusterB=guarentee_site_order(clusterB)

    #input two lists of lists with clusters and return cluster assignments
    CidxA=convert_sites_to_clust_idx(clusterA,active_sites)
    CidxB=convert_sites_to_clust_idx(clusterB,active_sites)
    #print('clusta',CidxA)
    #print('clustb',CidxB)

    #Compare each cluster to determine False/True Negatives/Positives
    ##Implement Rand Index to compare two clusters
    TP=0
    TN=0
    FP=0
    FN=0
    for i in range(len(CidxA)):
        for j in range(i+1,len(CidxA)):
            ##True postive
            if (CidxA[i]==CidxA[j]) and (CidxB[i]==CidxB[j]):
                TP+=1
            ##False negative
            if (CidxA[i]!=CidxA[j]) and (CidxB[i]==CidxB[j]):
                FN+=1
            ##False positive
            if (CidxA[i]==CidxA[j]) and (CidxB[i]!=CidxB[j]):
                FP+=1
            ##True negative
            if (CidxA[i]!=CidxA[j]) and (CidxB[i]!=CidxB[j]):
                TN+=1

    RandIdx=(TP+TN)/(TP+FP+TN+FN)

    #print("RI",RandIdx)
    return RandIdx
