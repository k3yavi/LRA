# -*- coding: utf-8 -*-
"""
Created on ......... Mon Nov 30 22:42:01 2015
Last modified on ... Thu Dec 03 10:59:05 2015

Three-step-clustering implementation.

@author: Caleb Andrade
"""

import time
from HierarchicalKmeans import *
from LocalitySensitive import *
from Cluster import Cluster
from matplotlib import pyplot as plt


def main():
    args = parse_args()
    reads = readFastq(args.infile, 10000)
    
    # parameters
    N = len(reads)
    min_bucket_size = 3
    kmeans_loops = 10
    error_t = 25 # error threshold (standard deviation)
        
    # first stage
    initial_clusters, k = firstStage(reads, min_bucket_size, N)
    del reads
    
    # second stage
    clusters = secondStage(initial_clusters, k, N, kmeans_loops)
    del initial_clusters

    # third stage
    new_clusters = thirdStage(clusters, N, error_t)

    # write clusters to file
    fileClusters(new_clusters)


def parse_args():
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('infile', help='the FASTQ read file')
        return parser.parse_args()


def statMeasures(values):
    """
    Compute the mean and variance of a list of values.
    
    Input: list of values.
    Output: average and standard deviation.
    """
    average = float(sum(values)) / len(values)
    variance = [(average - values[i])**2 for i in range(len(values))]
    
    return average, (float(sum(variance)) / (len(values) - 1))**0.5
    

def firstStage(reads, min_bucket_size, n):
    """
    Compute LSH for reads/kmers. Hash = concatenation of random primitive hashes.
        
    Input: list of reads/kmers, minimum size of buckets, number of reads.
    Output: list of clustered reads/kmers (as Cluster objects)
    """
    print "\n******************************************************************"
    print "\n        F I R S T  S T A G E: LOCALITY-SENSITIVE HASHING "
    
    clusters = []
    length = len(reads[0][1]) # length of read/kmer's as a bit-vector
    histogram = []
    k = 0
            
    tic = time.clock()
    # sample primitive hashes indices, and apply the respective LSH hash
    result = hashLSH(reads, random.sample(range(4*length), length)).values() 
    # initialize Cluster objects with hash buckets
    for bucket in result:
        clusters.append(Cluster(bucket))
        # append cluster sizes to histogram that are above a given threshold
        if min_bucket_size < len(bucket):
            # discard cluster sizes for histogram that are > 1% of total reads
            if len(bucket) < n*0.01: 
                histogram.append(len(bucket))
            else:
                k += 1 # increase the number of k clusters
    toc = time.clock()
            
    # how many initial clusters based on average and std. dev. of histogram
    avg, std_dev = statMeasures(histogram)
    for value in sorted(histogram, reverse = True):
        if value < avg + 2*std_dev:
            break
        k += 1 # increase the number of k clusters
        
    print "\nTotal number of reads : ", n
    print "Number of clusters:     ", len(clusters)
    print "Running time:           ", round(toc-tic, 2), " s"
    print "k: ", k
    print "\nCluster size average and std dev: ", round(avg, 2),", ", round(std_dev, 2)
    plt.plot(histogram)
    plt.show()
            
    return clusters, k


def secondStage(cluster_list, k, n, loops):
    """
    Second stage: apply kmeans clustering using k largest
    clusters from first stage.
    
    Input: list of clusters, number of k clusters, total number of reads,
    number of kmeans loops.
    Output: list of Cluster objects.
    """
    print "\n        S E C O N D  S T A G E: KMEANS-CLUSTERING "
    tic = time.clock()
    clusters = kmeansClustering(cluster_list, k, loops, False)
    toc = time.clock()
    printResults(clusters, n)
    print "Number of clusters:     ", len(clusters)
    print "\nRunning time:           ", round(toc-tic, 2), " s"
    
    return clusters
    
    
def thirdStage(clusters, n, error_t):
    """
    After kmeans we obtain k clusters, in this stage, we attempt
    to further merge clusters, but without surpassing an error threshold. 
    
    Input: list of clusters, total number of reads, error threshold.
    Output: list of cluster objects.
    """        
    print "\n        T H I R D  S T A G E: HIERARCHICAL-CLUSTERING"
    tic = time.clock()
    # Apply automatic hierarchical clustering for clusters from stage 2. 
    new_clusters = autoClustering(clusters, error_t)
    printResults(new_clusters, n)
    toc = time.clock()
    
    print "Number of clusters:     ", len(new_clusters)
    print "\nRunning time:           ", round(toc-tic, 2), " s"
    
    return new_clusters   
    

def fileClusters(cluster_list):
    """
    Writes to a file clusters' reads IDs.
    """
    text_file = open("clusters.txt", "w")
    for cluster in cluster_list:
        for ID in cluster.getIDs():
            text_file.write(ID + '\n')
        text_file.write('\n')
    text_file.close()


if __name__ == '__main__':
    main()
