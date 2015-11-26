# -*- coding: utf-8 -*-
"""
Created on ........ Tue Nov 17 20:26:48 2015
Last modified on .. Thu Nov 19 23:18:39 2015

Hierarchical and kmekclusters clustering implementation for RNA/DNA reads.

@author: Caleb Andrade
"""

from Cluster import Cluster
import time
import random
import math

# DIM = 1 manhattan metric; DIM = 2 euclidean metric
DIM = 1


def readFastq(filename):
    """
    Parse read and quality strings from a FASTQ file with sequencing reads.
    @author: Ben Langmead & Jacob Pritt.
    
    Input: file path
    Output: A list of reads, the list of qualities
    """
    sequences = []
    qualities = []
    
    with open(filename) as fh:
        while True:
            fh.readline() # skip name line
            seq = fh.readline().rstrip() #read base sequence
            fh.readline() # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
            
    return sequences, qualities
    

def closestPair(cluster_list):
    """
    Compute the closest pair from a list of clusters
    Note: Brute force method
    
    Input: List of clusters
    Output: Pair of clusters whose distance is minimal
    """
    minimum = float('inf')
    # loop through all possible pairs
    for idx in range(len(cluster_list) - 1):
        for jdx in range(idx + 1, len(cluster_list)):
            # calculate distance between pair of clusters
            dist = cluster_list[idx].distance(cluster_list[jdx])
            # pick the best so far
            if dist < minimum:
                minimum = dist
                best = (cluster_list[idx], cluster_list[jdx])
    #print "Minimum distance: ", minimum
    return best


def hierClustering(cluster_list, num_clusters, closestPair):
    """
    Compute a hierarchical clustering of a set of clusters.
    Note: the function may mutate cluster_list.
    
    Input: List of clusters, integer number of clusters, closest pair function.
    Output: List of clusters whose length is num_clusters
    """
    while len(cluster_list) > num_clusters:
        # find closest pair
        temp = closestPair(cluster_list)
        cluster1, cluster2 = temp[0], temp[1]
        # pop closest pair
        cluster_list.remove(cluster1)
        cluster_list.remove(cluster2)
        # merge closest pair and append to cluster_list
        cluster1.mergeClusters(cluster2)
        cluster_list.append(cluster1)
    return cluster_list
    
    
def autHierClustering(cluster_list, closestPair, threshold):
    """
    Compute a hierarchical clustering of a set of clusters.
    Note: number of clusters depend on a specified error threshold.
    
    Input: List of clusters, closest pair function.
    Output: List of clusters.
    """
    new_cluster_list = []
    while len(cluster_list) > 1:
        # find closest pair
        temp = closestPair(cluster_list)
        cluster1, cluster2 = temp[0], temp[1]
        # copy and merge to test distortion
        temp1, temp2 = cluster1.copy(), cluster2.copy()
        temp1.mergeClusters(temp2)
        # if within threshold procede
        if temp1.clusterError() < threshold:
            # pop closest pair
            cluster_list.remove(cluster1)
            cluster_list.remove(cluster2)
            # merge closest pair and append to cluster_list
            cluster1.mergeClusters(cluster2)
            cluster_list.append(cluster1)
        elif cluster1.clusterError() > cluster2.clusterError():
            new_cluster_list.append(cluster1)
            cluster_list.remove(cluster1)
        else:
            new_cluster_list.append(cluster2)
            cluster_list.remove(cluster2)
    if len(cluster_list) == 1:
        new_cluster_list.append(cluster_list.pop())
    return new_cluster_list
    

def kmeansClustering(cluster_list, num_clusters, num_iterations):
    """
    Compute the k-means clustering of a set of clusters (RNA/DNA reads)
    Note: the function may not mutate cluster_list
    
    Input: List of clusters, integer num_clusters = k and number of iterations
    Output: List of clusters whose length is num_clusters
    """
    kclusters = [] # this list to store k clusters to compare with (non-mutable)
    centroids = [] # this list to store the initial k centroids (average stats vectors)
    temp = list(cluster_list)
    random.shuffle(temp)
    
    # sample k initial random clusters to define initial centroids
    for cluster in temp[:num_clusters]:
        kclusters.append(cluster.copy())
        centroids.append(cluster.getAvgStats())
        
    for dummy_i in range(num_iterations):
        clusters = []
        # initialize new empty cluster objects at the centroids
        for idx in range(num_clusters):
            cluster = Cluster([])
            cluster.avg_stats_vectors = list(centroids[idx])
            clusters.append(cluster)
        
        # for every cluster in cluster_list
        for num in range(len(cluster_list)):
            best = (float('inf'), -1)
            # compare distance to every centroid at kclusters
            for idx in range(num_clusters):
                temp = cluster_list[num].distance(kclusters[idx])
                if temp < best[0]:
                    best = (temp, idx)
            # merge cluster to best centroid in list of mutable clusters
            clusters[best[1]].mergeClusters(cluster_list[num])
        
        # make a copy of re-calculated centroids: kclusters and centroids.
        for idx in range(num_clusters):
            kclusters[idx] = clusters[idx].copy()
            centroids[idx] = (clusters[idx].getAvgStats())
    return kclusters
    
#******************************************************************************
#                          T E S T I N G    Z O N E     
#******************************************************************************
def speciesPercent(results, num_reads, num_species):
    """
    Computes the percentage distribution of each species in cluster.
    
    Input: list of clusters, number of reads per species, number of species.
    Output: list of percentages for each cluster.
    """
    percentages = []
    #3 loop all clusters in results
    for cluster in results:
        ## initialize a percentages list per cluster
        count = num_species*[0]
        for label in cluster.getIDs():
            ## check to which species this label belongs to and increase count
            index = label / num_reads
            count[index] += 1
        ## convert counts to relative percentages
        if cluster.getSize() != 0:
            temp = [round(float(100*num) / cluster.getSize(), 1) for num in count]
            percentages.append(temp)
        else: 
            percentages.append(count)
    
    return percentages
    
def printResults(results, num_reads, num_species, tic, toc):
    """
    Display results.
    
    Input: List of clusters, number of reads per species, number of species,
    initial time, finishing time.
    """
    count = 1
    percent = speciesPercent(results, num_reads, num_species)
    for cluster in results:
        print '\n\n                       C L U S T E R ', count, '\n'
        info = "Species: "
        info += str(percent[count-1][0]) + '%  '
        info += str(percent[count-1][1]) + '%  '
        info += str(percent[count-1][2]) + '%  '
        info += ",  Reads: " + str(cluster.getSize())
        info += ",  Distortion: " + str(round(cluster.clusterError(), DIM)) + '\n'
        print info
        IDs = sorted(list(cluster.getIDs()))
        print IDs
        count += 1
    print "\nRunning time: ", toc - tic
    

def main(num_reads, num_clusters):
    """
    Testing.
    
    Input: number of reads per species to be clustered, number of clusters.
    Output: clusters of readings IDs
    """
    read0 = readFastq('ads1_week4_reads.fq')[0]
    read1 = readFastq('ERR037900_1.first1000.fastq')[0]
    read2 = readFastq('ERR266411_1.for_asm.fastq')[0]
    
    ## Mix reads
    all_reads = read0[:num_reads] + read1[:num_reads] + read2[:num_reads]
    
    ## label reads
    reads = [(idx, all_reads[idx]) for idx in range(len(all_reads))]
    
    ## create initial list of clusters
    cluster_list1 = [Cluster([read]) for read in reads] 
    cluster_list2 = [Cluster([read]) for read in reads]        
    cluster_list3 = [Cluster([read]) for read in reads]  
    
    print "\n********************** KMEANS CLUSTERING ************************"
    tic = time.clock()
    iterations = 1 + int(math.log(len(cluster_list1), 2))
    results = kmeansClustering(cluster_list1, num_clusters, iterations)
                               
    toc = time.clock()
    printResults(results, num_reads, 3, tic, toc)
        
    print "\n******************** HIERARCHICAL CLUSTERING ********************"
    tic = time.clock()
    results = hierClustering(cluster_list2, num_clusters, closestPair)
    toc = time.clock()
    printResults(results, num_reads, 3, tic, toc)
    
    print "\n****************** AUT HIERARCHICAL CLUSTERING ******************"
    tic = time.clock()
    results = autHierClustering(cluster_list3, closestPair, 20)
    toc = time.clock()
    printResults(results, num_reads, 3, tic, toc)

#main(10, 3)    
#main(100, 3)
#main(100, 10)
#main(500, 30)
    


            
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    