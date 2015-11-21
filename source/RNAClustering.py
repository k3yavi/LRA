# -*- coding: utf-8 -*-
"""
Created on ........ Mon Nov 09 00:03:48 2015
Last modified on .. Wed Nov 18 19:43:10 2015

RNA Clustering project.
Some methods useful for DNA and RNA reads sequencing & clustering algorithms.

@author: Caleb Andrade (Unless stated otherwise)
"""

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
    

def editDistance(x, y):
    """
    Compute edit distance between two strings using Dynamic Programming.
    @author: Ben Langmead & Jacob Pritt.
    
    Input: string x, string y.
    Output: edit distance.
    """
    D = [] # define array to store values
    
    for i in range(len(x)+1): # create rows
        D.append([0]*(len(y)+1))
    for i in range(len(x)+1): # initialize first column
        D[i][0] = i
    for i in range(len(y)+1): # initialize first row
        D[0][i] = i
    
    # Dynamic programming stage
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min (distHor, distVer, distDiag)

    return D[-1][-1]
    
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
    
def readStats(reads):
    """
    Create a list of stats vectors for reads, each vector has the form:
    (#A's, #C's, #G's, #T's, d(A's), d(C's), d(G's), d(T's), #switches, errors)
    
        * # A's means the total number of A's in the read.
        * d(A's) means the total distance in between A's.
        * #switches means the total number of base alternations in a read.
        * errors means the number of undefined 'N' bases in a read.
        
    Note: With this information we can encode a read as a 10-tuple, so that 
    the manhattan/euclidean metrics can be used to measure closeness between
    the properties encoded in the stats vector of any two reads.
    
    Input: list of reads.
    Output: list of stats vectors corresponding to reads in the list of reads.
    """
    read_stats = []
    for read in reads:
        # initialize counters for each base (A, C, G, T)
        count = {'A':0, 'C':0, 'G':0, 'T':0}
        # 'A':[a, b]: a = distance in between A's, b = index of last 'A'
        dist = {'A':[0, None], 'C':[0, None], 'G':[0, None], 'T':[0, None]}
        num_switch = 0 # counter to keep track of alternation of bases
        prev_base = ''
        index = 0 # to indicate current's base index
        errors = 0 # number of N's in the read
        # loop through read's bases 
        for base in read:
            # update number of switches
            if base != prev_base:
                num_switch += 1
                prev_base = base
            # if an 'N' appears, we skip it
            if base == 'N':
                index += 1
                errors += 1
                continue
            # update number of bases
            count[base] += 1
            # update distance in between same bases
            if dist[base][1] != None:
                dist[base][0] += index - dist[base][1] - 1
            dist[base][1] = index
            index += 1
        read_stats.append((count['A'], count['C'], count['G'], count['T'],
                        dist['A'][0], dist['C'][0], dist['G'][0], dist['T'][0],
                        num_switch - 1, errors))
    return read_stats
    
def manhattan(stats_a, stats_b):
    """
    Compute manhattan distance of two read, represented as stats vectors.
    
    Input: two stats vectors (10-tuple).
    Output: manhattan distance of stats vectors.
    """
    manhattan = 0
    for idx in range(len(stats_a)):
        manhattan += abs(stats_a[idx] - stats_b[idx])
    return manhattan
    
def euclidean(stats_a, stats_b):
    """
    Compute euclidean distance of two reads, represented as stats vectors.
    
    Input: two stats vectors (10-tuple).
    Output: euclidean distance of stats vectors.
    """
    euclidean = 0
    for idx in range(len(stats_a)):
        euclidean += (stats_a[idx] - stats_b[idx])**2
    return euclidean**0.5
    
