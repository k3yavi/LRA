# -*- coding: utf-8 -*-
"""
Created on ......... Mon Nov 30 22:42:01 2015

3LayerClustering implementation.

@author: Caleb Andrade
"""

import time
from HierarchicalKmeans import *
from LocalitySensitive import *
from Cluster import Cluster
from matplotlib import pyplot as plt


def main():
    args = parse_args()
    reads = readFastq(args.infile)
    
    # parameters
    N = len(reads)
    min_bucket_size = 2
    max_bucket_size = 0.01*N
    kmeans_loops = 5
    num_hashes = 50 # number of primitive hashes
        
    # first stage
    initial_clusters, k = firstStage(reads, min_bucket_size, N, num_hashes, max_bucket_size)
    del reads
    
    # second stage
    clusters, error_t = secondStage(initial_clusters, k, N, kmeans_loops)
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
    Compute mean and standard deviation of a list of values.
    """
    average = float(sum(values)) / len(values)
    variance = [(average - values[i])**2 for i in range(len(values))]
    
    return average, (float(sum(variance)) / (len(values) - 1))**0.5
    

def firstStage(reads, min_bucket_size, n, num_hashes, max_bucket_size):
    """
    Compute LSH for reads. Hash = concatenation of random primitive hashes.
        
    Input: list of reads, minimum size of buckets, number of reads, 
    number of primitive hashes, maximum size of buckets.
    Output: list of clustered reads (as Cluster objects)
    """
    print "\n******************************************************************"
    print "\n        F I R S T  S T A G E: LOCALITY-SENSITIVE HASHING "
    
    clusters = []
    histogram = []
    length = len(reads[0][1]) 
            
    tic = time.clock()
    # sample primitive hashes indices, apply the respective LSH hash
    buckets = hashLSH(reads, random.sample(range(4*length), num_hashes)) 
    # initialize Cluster objects with hash buckets
    for bucket in buckets.values():
        cluster = Cluster(bucket)
        clusters.append(cluster)
        # append to histogram cluster sizes that are above a given threshold
        if min_bucket_size < len(bucket):
            # discard cluster sizes for histogram that are too big
            if len(bucket) < max_bucket_size:
                histogram.append(len(bucket))
#                print
#                for cluster in bucket:
#                    print cluster
    
    toc = time.clock()
    avg_error = avgError(clusters, n)       
    
    # how many initial clusters based on average of  size histogram
    k = 0
    avg, std_dev = statMeasures(histogram)
    for value in sorted(histogram, reverse = True):
        k += 1 # increase the number of k clusters
        if value < avg:
            break 
    
    print "\nNumber of reads:        ", n
    print "\nNumber of LSH buckets:  ", len(clusters)
    print "Bucket 'avg size':      ", round(avg, 2)
    print "Weighted Average Error: ", round(avg_error, 2) 
    print "Running time:           ", round(toc-tic, 2), " s"
    print "kmeans clusters (k):    ", k
    plt.plot(histogram)
    plt.show()
            
    return clusters, k


def secondStage(cluster_list, k, n, loops):
    """
    Kmeans clustering using k largest clusters from first stage.
    
    Input: list of clusters, number of k clusters, total number of reads,
    number of kmeans loops.
    Output: list of Cluster objects.
    """
    print "\n\n        S E C O N D  S T A G E: KMEANS-CLUSTERING "
    
    tic = time.clock()
    clusters = kmeansClustering(cluster_list, k, loops, False)
    toc = time.clock()
    avg_error = avgError(clusters, n)
    print "\nNumber of clusters:     ", len(clusters)
    print "Weighted Average Error: ", round(avg_error, 2) 
    print "Running time:           ", round(toc-tic, 2), " s"
    print "Kmeans loops:           ", loops
    
    return clusters, avg_error
    
    
def thirdStage(clusters, n, error_t):
    """
    Further clustering, without surpassing an error threshold. 
    
    Input: list of clusters, total number of reads, error threshold.
    Output: list of cluster objects.
    """        
    print "\n\n        T H I R D  S T A G E: HIERARCHICAL-CLUSTERING"
    
    tic = time.clock()
    # Apply automatic hierarchical clustering for clusters from stage 2. 
    new_clusters = autoClustering(clusters, error_t)
    avg_error = avgError(new_clusters, n)
    toc = time.clock()
    
    print "\nNumber of clusters:     ", len(new_clusters)
    print "Weighted Average Error: ", round(avg_error, 2) 
    print "Running time:           ", round(toc-tic, 2), " s"
    
    return new_clusters   
    

def fileClusters(cluster_list):
    """
    Writes to a file.
    """
    text_file = open("clusters.txt", "w")
    for cluster in cluster_list:
        for ID in cluster.getIDs():
            text_file.write(ID + '\n')
        text_file.write('\n')
    text_file.close()


if __name__ == '__main__':
    main()
