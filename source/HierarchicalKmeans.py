# -*- coding: utf-8 -*-
"""
Created on ......... Tue Nov 17 20:26:48 2015
Last modified on ... Mon Nov 30 23:52:44 2015
Hierarchical and kmeans clustering implementation for RNA/DNA reads.
@author: Caleb Andrade
"""

from Cluster import Cluster
import random

# DIM = 1 manhattan metric; DIM = 2 euclidean metric
DIM = 1


def readFastq(filename, limit = float('inf')):
    """
    Parse read and quality strings from a FASTQ file with sequencing reads.
    @author: Ben Langmead & Jacob Pritt.
    
    Input: file path, limit of number of reads to extract from file
    Output: The list of (id, reads), the list of qualities.
    """
    sequences = []
    qualities = []
    count = 0 # counts lines
    
    with open(filename) as fh:
        while count < limit:
            first_line = fh.readline()
            name = first_line[1:].rstrip() # name line
            seq = fh.readline().rstrip() # read base sequence
            fh.readline() # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append((name, seq))
            qualities.append(qual)
            count += 1
            
    return sequences, qualities


def pairedFastq(filename1, filename2, limit = float('inf')):
    """
    Parse paired-end reads from a pair of FASTQ filehandles
    For each pair, we return a name, the nucleotide string
    for the first end, the quality string for the first end,
    the nucleotide string for the second end, and the
    quality string for the second end.
    
    @author: Ben Langmead
    @adaptation: Caleb Andrade
    """
    reads = []
    count = 0 # counts lines
    
    fh1 = open(filename1) 
    fh2 = open(filename2)
    
    while count < limit:
        first_line_1 = fh1.readline()
        first_line_2 = fh2.readline()
        if len(first_line_1) == 0:
            break  # end of file
        name_1, name_2 = first_line_1[1:].rstrip(), first_line_2[1:].rstrip()
        seq_1, seq_2 = fh1.readline().rstrip(), fh2.readline().rstrip()
        fh1.readline()  # ignore line starting with +
        fh2.readline()  # ignore line starting with +
        fh1.readline()  # ignore qualities line
        fh2.readline()  # ignore qualities line
        reads.append([(name_1, seq_1), (name_2, seq_2)])
        count += 1
    
    return reads


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


def hierClustering(cluster_list, num_clusters):
    """
    Compute a hierarchical clustering of a set of clusters.
    Note: the function may mutate cluster_list.
    
    Input: List of clusters, integer number of clusters.
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
    
    
def autoClustering(cluster_list, threshold):
    """
    Compute a hierarchical clustering of a set of clusters.
    Note: number of clusters depend on a specified error threshold.
    
    Input: List of clusters, error threshold.
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
    

def kmeansClustering(cluster_list, k, iterations, shuffle = True):
    """
    Compute the k-means clustering of a set of clusters (RNA/DNA reads)
    Note: the function may not mutate cluster_list
    
    Input: List of clusters, k number of clusters, iterations, select randomly initial clusters?
    Output: List of clusters whose length is num_clusters
    """
    kclusters = [] # this list to store k clusters to compare with (non-mutable)
    centroids = [] # this list to store the initial k centroids (average stats vectors)
    if shuffle:
        # shuffle cluster list
        random.shuffle(cluster_list) 
    else:
        # sort it by size
        cluster_list.sort(key = lambda cluster: cluster.getSize(), reverse = True)

    # k initial clusters to define initial centroids
    for cluster in cluster_list[:k]:
#        print "Cluster size: ", cluster.getSize()
        kclusters.append(cluster.copy())
        centroids.append(cluster.getAvgStats())
        
    for iteration in range(iterations):
        clusters = []
        # initialize new empty cluster objects at the centroids
        for idx in range(k):
            cluster = Cluster([])
            cluster.avg_stats_vectors = list(centroids[idx])
            clusters.append(cluster)
        
        # for every cluster in cluster_list
        for num in range(len(cluster_list)):
            best = (float('inf'), -1)
            # compare distance to every centroid at kclusters
            for idx in range(k):
                temp = cluster_list[num].distance(kclusters[idx])
                if temp < best[0]:
                    best = (temp, idx)
            # merge cluster to best centroid in list of mutable clusters
            clusters[best[1]].mergeClusters(cluster_list[num])
        
        # make a copy of re-computed centroids: kclusters and centroids.
        for idx in range(k):
            kclusters[idx] = clusters[idx].copy()
            centroids[idx] = (clusters[idx].getAvgStats())
    
    return kclusters

def printResults(clusters, n):
    """
    Display results.
    
    Input: List of clusters.
    """
    mean_error = 0
    for cluster in clusters:
        temp = cluster.clusterError()
        mean_error += temp*cluster.getSize()
    
    print "\nWeighted Average Error: ", round(mean_error / n, 2) 
    
