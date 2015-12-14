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
    Output: The list of (id, reads).
    """
    sequences = []
    count = 0 # counts lines
    
    with open(filename) as fh:
        while count < limit:
            first_line = fh.readline()
            name = first_line[1:].rstrip() # name line
            seq = fh.readline().rstrip() # read base sequence
            fh.readline() # skip placeholder line
            fh.readline().rstrip() # base quality line, ignore it
            if len(seq) == 0:
                break
#            name += str(count)
            sequences.append((name, seq))
            count += 1
            
    return sequences


def closestPair(cluster_list):
    """
    Compute the closest pair from a list of clusters. Brute force method
    
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

    return best


def hierClustering(cluster_list, num_clusters):
    """
    Hierarchical clustering for a set of clusters. Mutate cluster_list.
    
    Input: List of initialized clusters, desired number of clusters.
    Output: List of clusters of size num_clusters.
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
    Hierarchical clustering for a set of clusters.
    Note: number of clusters depend on a specified error threshold.
    
    Input: List of clusters, error threshold.
    Output: List of clusters.
    """
    new_cluster_list = []
    while len(cluster_list) > 1:
        # find closest pair
        cluster1, cluster2 = closestPair(cluster_list)
        # copy and merge to test distortion
        temp1, temp2 = cluster1.copy(), cluster2.copy()
        temp1.mergeClusters(temp2)
        # if within threshold procede
        if temp1.clusterError() < threshold:
            # remove closest pair from cluster_list
            cluster_list.remove(cluster1)
            cluster_list.remove(cluster2)
            # merge closest pair and append it to cluster_list
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
    Compute the k-means clustering of a set of clusters (reads/kmers)
    Note: the function may not mutate cluster_list
    
    Input: List of clusters, k number of clusters, iterations, 
    select initial clusters: randomly or by size?
    Output: List of clusters.
    """
    kclusters = [] # this list to store k clusters to compare with (non-mutable)
    centroids = [] # this list to store the initial k centroids (average stats vectors)
    
    if shuffle:
        # shuffle cluster list
        random.shuffle(cluster_list) 
    else:
        # sort by size
        cluster_list.sort(key = lambda cluster: cluster.getSize(), reverse = True)

    # k initial clusters to define initial centroids
    for cluster in cluster_list[:k]:
        kclusters.append(cluster.copy())
        centroids.append(cluster.getAvgAbundance())
        
    for iteration in range(iterations):
        clusters = []
        # initialize new empty cluster objects at the centroids
        for idx in range(k):
            cluster = Cluster([])
            cluster.avg_abundance_vectors = list(centroids[idx])
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
            centroids[idx] = (clusters[idx].getAvgAbundance())
    
    return kclusters


def avgError(clusters, n):
    """
    Compute weighted average error.
    """
    mean_error = 0
    for cluster in clusters:
        temp = cluster.clusterError()
        mean_error += temp*cluster.getSize()
    
    return mean_error / n
    
