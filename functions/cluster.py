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
    #set number of clusters
    numClusters=2
    numSites=0
    for x in active_sites:
        numSites+=1
    
    #create matrix with feature values
    features=[]
    for z in active_sites:
        features.append(compute_hydrophobicity_Index(z))

    print('make feat',features)

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
        
        #pos=pos.tolist()
        
        pos=float(str(pos[0]))
        print('pos after',pos)
        #print('assign cluster')
        #print(pos)
        centroid.append(pos)
    print('centroid',centroid)
    Cidx=[]
    print('sites', active_sites)
    
    #define number of iterations
    for iteration in range(0,1):

        #assign clusters
        #list of indicies of closest cluster for each active site
        idxtemp=[0]*numSites
        
        #for each active site, assign it a cluster based on distance
        for j in range(0,len(idxtemp)):
            print('test access',active_sites[j].residues)
            minimum=float("inf")
            #site=features[j]
            site=active_sites[j]
            print('len',len(centroid))
            for i in range(0,len(centroid)):
                #this is a stupid work around do better
                for item in active_sites:
                    if centroid[i]==(compute_hydrophobicity_Index(item)):
                            tempActive=item
                test=compute_similarity(site, tempActive)
                print('test',test)
                if test <= minimum:
                    minimum=test
                    idxtemp[j]=i
        
        print('temp ind',idxtemp)
        print("ctemp",Cidx)
        #if iteration==0:
         #   Cidx=idxtemp
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
        
        #print('c',C)
        #recalculate centers for each cluster
        for j in range(0,len(centroid)):
            nums=idxtemp.count(j)
            print('num', nums)
            vect=[0]*nums
            print('vect',vect)
            counter=0
            print('len',len(Cidx))
            for i in range(0,len(Cidx)):
                print(Cidx[i],j)
                if (Cidx[i]==j):
                    vect[counter]=float(features[i])
                    print('feat',features[i])
                    print(vect)
                    counter+=1
            print('vect2',vect)
            if vect==[]:
                vect.append(0)
            avg=np.mean(vect)
            centroid[j]=float(avg)
            print(avg)
        print('recalc',centroid)
        #print('c', C)
        print('idx',idxtemp)
        print('idx',sorted(idxtemp))
        return sorted(idxtemp)
        """
        iter=0
        for m in centroid:
            #count number of samples in each cluster
            counts=0
            clust=[]
            for item in active_sites:
                if item == m:
                    counts+=1
                    clust.append(item)
            centroid[iter]=(1/counts)*sum(clust)
        """ 
    

def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #
    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    #implement Agglomerative Algorithm
    listoflists=[]
    '''
    #put each object in its own cluster
    for site in active_sites:
        a_list.append(site)
        listoflists.append(a_list)
        a_list=[]
    
    #iterate through list compare distance and merge cells with closest distance
    for item in listoflists:
        #compare first item to all others merge shortest distance

    #inplement Nearest Neighbor approach
    '''


    return []
