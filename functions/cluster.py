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
    print(features)
   
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
        
        print('centroid',centroid)

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
        print('idx',idxtemp)
        
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
                print('comp',j,item)
                if j == item:
                    #append distance to recalculate centroid
                    #append active site numbers to show these clusters
                    avgvect[count]=features[i]
                    clustvect[count]=active_sites[i]
                    count+=1
                i+=1
            print('vect',avgvect)
            newcent=np.mean(avgvect)
            print("new",newcent)
            centroid[centid]=newcent
            print("new cent vect",centroid)
            centid+=1
            #create list of list for current clusters
            tempclust.append(clustvect)
        print(tempclust)

        #####################################################
        #determine if clusters stayed the same this iteration
        #if so return current clsuters, else keep iterating
        #####################################################

        if tempclust == Cidx:
            clustfinal=tempclust
            print('final found',clustfinal)
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
    print('aggs')
    
    ##################################
    #implement Agglomerative Algorithm
    ##################################

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
        #determine the indicies to merge
        idxlist=[]
        idxlist.append(listoflists.index(clusterA))
        idxlist.append(listoflists.index(clusterB))
        idxlist=sorted(idxlist)
    
        #merge the indicies given
        #for ind in idxlist:
        mergeitems=listoflists[idxlist[0]]+listoflists[idxlist[1]]
        print(mergeitems)
        listoflists.remove(listoflists[idxlist[0]])
        idxlist[1]=idxlist[1]-1
        listoflists.remove(listoflists[idxlist[1]])
        listoflists=listoflists+[mergeitems]
        print(listoflists) 
    print('final',listoflists)
    return listoflists

def check_clust_quality(clusterlist):
    print(clusterlist)
    totclust=0
    for clust in clusterlist:
        clustdist=0
        for i in range(len(clust)):
            for j in range(i+1,len(clust)):            
                disttemp=compute_similarity(clust[i],clust[j])
                print('sim',disttemp)
                clustdist+=disttemp
                print('clust dist',clustdist)
        clustdist=clustdist/len(clust)
        print(clustdist)
        totclust+=clustdist 
    
    print('all clust',totclust)
    final=totclust/len(clusterlist)
    print(final)
    return final


    print(clusterlist)
    return clusterlist



