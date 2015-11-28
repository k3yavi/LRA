# -*- coding: utf-8 -*-
"""
Created on ......... Wed Nov 25 23:35:58 2015
Last modified on ... Fri Nov 27 21:55:55 2015

Locality-Sensitive Hashing applied to cluster RNA/DNA kmers.

@author: Caleb Andrade.
"""

from HierarchicalKmeansCluster import *

def kmerHashMap(reads, k):
    """
    Create a hash map between kmers and readings. 
    Note: Keys are kmers. Values are sets of those reads containing them.
    
    Input: list of reads, length of kmers.
    Output: hash map.
    """
    kmers_dict = {}
    # loop through all reads
    for i in range(len(reads)):
        # loop read's bases, except for the last k, to obtain its kmers
        for j in range(1+len(reads[i])-k):
            kmer = reads[i][j:k+j]
            if kmers_dict.has_key(kmer):
                kmers_dict[kmer].add(i)
            else:
                kmers_dict[kmer] = set([i])
    
    return kmers_dict
    

def unaryMap(kmer):
    """
    Map kmer/read to a unary representation.
    
    Input: kmer/read.
    Output: concatenation of unary maps of kmer/read's bases.
    """
    bit_vector = ''
    for base in kmer:
        bit_vector += UNARY_MAP[base]
    return bit_vector
    

def simHash(bit_vector, ran_numbers):
    """
    Similarity hash.
    
    Input: bit vector, sorted list of k random numbers.
    Output: simHash index.
    """
    ans = ''
    for num in ran_numbers:
        ans += bit_vector[num]
    return ans
    

def hashLSH(kmers, ran_numbers):
    """
    Build the similarity hash table according to the given indices.
    
    Input: set of kmers/reads (ID, 'ACCTGCA...'), list of random integers.
    Output: kmer/read similarity hash map.
    """
    hashMap = {}
    for kmer in kmers:
        temp = simHash(unaryMap(kmer[1]), ran_numbers)
        if hashMap.has_key(temp):
            hashMap[temp].append(kmer)
        else:
            hashMap[temp] = [kmer]
    return hashMap

def compress(bucket):
    """
    Take the average of reads in bucket.
    
    Input: Bucket of raw reads as strings.
    Output: Average read as string.
    """
    acgt = {0:'A', 1:'C', 2:'G', 3:'T', 4:'T'}
    sumVector = 4*len(bucket[0])*[0]
    average =''    
    
    for read in bucket:
        temp = unaryMap(read)
        for idx in range(len(temp)):
            sumVector[idx] += int(temp[idx])
    
    best = 0
    index = 0
    count = 0
    
    for idx in range(len(sumVector)):
        if sumVector[idx] > best:
            index = idx
            best = sumVector[idx]
        count += 1
        if count == 4:
            average += acgt[index%4]
            best = 0
            index = 0
            count = 0
            
    return average
#******************************************************************************
# TESTING ZONE
#******************************************************************************
UNARY_MAP = {'A':'1000', 'C':'0100', 'G':'0010', 'T':'0001', 'N':'0000'}

read0 = readFastq('ads1_week4_reads.fq')[0]
read1 = readFastq('ERR037900_1.first1000.fastq')[0]
read2 = readFastq('ERR266411_1.for_asm.fastq')[0]


def main1(list_reads, bucket_t, error_t, kmers = True):
    """
    Test LSH for reads/kmers. Hash = concatenation of random primitive hashes.
    For kmers the number of primitive hashes is set to 32, 4 cycles of hashing.
    For reads the number of primitive hashes is set to 50, 8 cycles of hashing.
    
    Input: list of reads/kmers, kmers or reads?
    Output: list of clusters, list of reads not clustered.
    """
    print "\n******************************************************************"
    print "\n        F I R S T  S T A G E: LOCALITY-SENSITIVE HASHING "
    if kmers:
        k = 16 # are we analyzing kmers?
    else:
        k = 25 # are we analyzing reads?
    # initialization
    count = 0
    reads = [(idx, list_reads[idx]) for idx in range(len(list_reads))]
    clusters = []
    reads_clustered = 0
    bit_length = 4*len(list_reads[0]) # length of the read/kmer as a bit-vector
    indices = set([i for i in range(bit_length)]) # indices to sample for primitive hashes
    
    tic = time.clock()
    # hash, filter and rehash those isolated reads 
    while count < bit_length / k: # loop to re-apply LSH 4 or 8 cycles of hashing.
        ran_numbers = set(random.sample(indices, k)) # sample indices
        indices = indices.difference(ran_numbers) # substract sampled indices
        result = hashLSH(reads, ran_numbers) 
        temp = []
        for bucket in result.values(): # filter unclustered reads
            if len(bucket) > bucket_t:
                cluster = Cluster(bucket)
                error = cluster.clusterError()
                if error < error_t:
                    clusters.append(Cluster(bucket))
                    reads_clustered += len(bucket)
#                    print "\nCluster error: ", error
#                    for item in bucket:
#                        print item[1]
                else:
                    temp += bucket
            else:
                temp += bucket                
        reads = temp
        count += 1
    toc = time.clock()
    
    print "\nSize of batch of reads: ", len(list_reads)
    print "Number of clusters:     ", len(clusters)
    print "Total reads clustered:  ", reads_clustered
    print "\nRunning time:           ", round(toc-tic, 2), " s"
    
    return clusters, reads


def main2(list_reads, bucket_t, error_t, stage3 = False):
    """
    Test a two stage method for clustering reads: Apply first LSH to 
    generate K initial clusters, and the rest of reads not clustered. Then 
    apply kmeans clustering, where iterations = log(n).
    
    Input: list of reads.
    Output: list of cluster objects.
    """
    # first stage of clustering
    results = main1(list_reads, bucket_t, error_t, False)

    print "\n        S E C O N D  S T A G E: KMEANS-CLUSTERING "
    tic = time.clock()
    not_clustered = [Cluster([read]) for read in results[1]]
    # initializing Cluster objects
    cluster_list = results[0] + not_clustered
    # seting the number of iterations as log(n)
    iterations = 1 + int(math.log(len(list_reads), 2))
    results2 = kmeansClustering(cluster_list, len(results[0]), iterations, False)
    toc = time.clock()
    printResults(results2, list_reads)
    print "Number of clusters:     ", len(results2)
    print "\nRunning time:           ", round(toc-tic, 2), " s"
    
    if not stage3:
        return results2
        
    print "\n        T H I R D  S T A G E: HIERARCHICAL-CLUSTERING"
    tic = time.clock()
    results3 = autHierClustering(results2, 20)
    printResults(results3, list_reads)
    toc = time.clock()
    print "Number of clusters:     ", len(results3)
    print "\nRunning time:           ", round(toc-tic, 2), " s"
    
    return results3    
    

def main3(k, list_reads, bucket_t, error_t):
    """
    Test a two stage method for clustering reads: Apply first kmeans
    with k initial clusters, then apply automatic hierarchical
    clustering with a threshold = 20
    
    Input: list of reads.
    Output: list of cluster objects.
    """
    results = kmeans(1000, 50)
    
    print "\n        S E C O N D  S T A G E: HIERARCHICAL-CLUSTERING"
    tic = time.clock()
    results2 = autHierClustering(results, 20)
    printResults(results2, list_reads)
    toc = time.clock()
    print "Number of clusters:     ", len(results2)
    print "\nRunning time:           ", round(toc-tic, 2), " s"
    
    return results2
    
    
def main4(kmers, bucket_t):
    """
    Compute an LSH hash for a list of kmers, using 4 cycles of 32 primitive 
    hash functions, sampling every time 32 indices without replacement for
    128-bit vectors.
    
    Input: List of kmers.
    Output: List of clusters of kmers. Those clusters with less than a bucket 
    threshold of kmers are ignored.
    """
    kmerHash1 = kmerHashMap(kmers, 32) # kmers-reads hash map
    kmers = kmerHash1.keys()
    results = main1(kmers, bucket_t, 20)
    print "\nKMER APPROXIMATE CLUSTERING VIA LSH (PRIMITIVE HASHES PARADIGM)"
    
    return results[0], kmers


x = main1(read1, 5, 20, False)
y = main2(read1, 5, 20)
z = main2(read1, 5, 20, True)
u = main3(50, read1, 5, 20)
w = main4(read1, 5)

def printReads(cluster_list, reads):
    """
    Print clusters of reads with cluster error.
    Input: list of clusters
    """
    for cluster in cluster_list:
        print "\nCluster error: ", cluster.clusterError()
        for idx in cluster.getIDs():
            print reads[idx]

#printReads(x[0], read1)
#printReads(y, read1)
#printReads(z, read1)
#printReads(u, read1)
#printReads(w[0], w[1])

def main5():
    """
    Testing compress.
    """
    bucket = ['CCTAACANTAACCCTAACCCCTAACCCTAACC',
              'CCTAACCCTAAACCCTAAACCCTAAACCCTAA',
              'CCTAACCCTAAACCTAACCTCTCACCCTAACC',
              'CCTAACCCCTAACCCTAACCCCTAACCCCTNA',
              'CCTAACCCTAACCCTAAACCCTAAACCCTAAA',
              'CCTAACCCTAACCCTTAACCCTTAANCCTTAA',
              'CCTAACCCTAACCCTAACANCCCTAACCCTAA',
              'CCTAACCCTAACCCTAACCCTAACCACCCTAA',
              'CTTAACCCTTAANCCTTAACCCTTAACCCTAA',
              'CCTAACCCTAACCCTAACCCTTAACCCTTAAN',
              'CCTAACCCTAACCCTAACCCTTAACCCTAACN',
              'CCTAACCCTNACCCTAACCCTTAACCCTAACC',
              'CCTAACCCTTACCCTTAACCCTTAACCCTNAA']
     
    print compress(bucket)
        
main5()
