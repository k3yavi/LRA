# -*- coding: utf-8 -*-
"""
Locality-Sensitive Hashing applied to cluster RNA/DNA kmers.
Created on Wed Nov 25 23:35:58 2015
@author: Caleb Andrade.
"""

import random

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
    Map kmer to a unary representation.
    
    Input: kmer.
    Output: concatenation of unary maps of kmer's bases.
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
    Builds the similarity hash table according to the given indices.
    
    Input: set of kmers, list of random integers.
    Output: kmer similarity hash map.
    """
    hashMap = {}
    for kmer in kmers:
        temp = simHash(unaryMap(kmer), ran_numbers)
        if hashMap.has_key(temp):
            hashMap[temp].append(kmer)
        else:
            hashMap[temp] = [kmer]
    return hashMap


#******************************************************************************
# TESTING ZONE
#******************************************************************************
UNARY_MAP = {'A':'1000', 'C':'0100', 'G':'0010', 'T':'0001', 'N':'0000'}

read0 = readFastq('ads1_week4_reads.fq')[0]
read1 = readFastq('ERR037900_1.first1000.fastq')[0]
read2 = readFastq('ERR266411_1.for_asm.fastq')[0]


def main1(k):
    """
    Testing LSH for 32-mers. Hash is a concatenation of random primitive hashes.
    
    Input: number of primitive hashes (best results ~ k =32).
    """
    kmerHash1 = kmerHashMap(read1, 32) # kmers-reads hash map
    kmers = kmerHash1.keys()
    # each kmer maps to a 128-bit vector, so we need a 128-long random numbers list
    ran_numbers = sorted(random.sample([i for i in range(128)], k))
    result = hashLSH(kmers, ran_numbers) # apply LSH hash
    
    print "\nThose buckets with more than 4 clustered kmers: \n"
    for bucket in result.values():
        if len(bucket) > 4:
            for kmer in bucket:
                print kmer
            print
    
    print "\nSize of reads data set: ", len(read1)
    print "\nInitial count of different kmers: ", len(kmers)
    print "\nNumber of kmer clusters: ", len(result)
    print "\nIndices of primitive hashes: \n\n", ran_numbers


def main2(k, iterations):
    """
    Testing LSH for reads. Hash is a concatenation of random primitive hashes.
    
    Input: number of primitive hashes (best results ~ k = 25).
    """
    # initialization
    count = 0
    temp = []
    reads = read1
    clusters = []
    reads_clustered = 0
    # hash, filter and rehash those isolated reads, loop # iterations.
    while count < iterations:
        ran_numbers = sorted(random.sample([i for i in range(400)], k))
        result = hashLSH(reads, ran_numbers)
        
        for bucket in result.values():
            if len(bucket) > 4:
                print
                clusters.append(bucket)
                reads_clustered += len(bucket)
                for kmer in bucket:
                    print kmer
            else:
                temp += bucket                
        reads = list(temp)
        temp = []
        count += 1
    
    print "\nSize of reads data set: ", len(read1)
    print "\nNumber of read clusters: ", len(clusters)
    print "\nTotal reads clustered: ", reads_clustered
    

#main1(32)

#main2(25, 25)

    
    



        
    




































