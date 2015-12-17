# -*- coding: utf-8 -*-
"""
Created on ......... Wed Nov 25 23:35:58 2015
Last modified on ... Mon Dec 14 02:32:44 2015

Locality-Sensitive Hashing applied to cluster RNA/DNA reads.

@author: Caleb Andrade.
"""

from HierarchicalKmeans import *

UNARY_MAP = {'A':'1000',
             'a':'1000',
             'C':'0100',
             'c':'0100',
             'G':'0010',
             'g':'0010',
             'T':'0001',
             't':'0001',
             'N':'0000',
             'n':'0000'}
             

def unaryMap(read):
    """
    Map read to a binary representation.
    
    Input: read.
    Output: concatenation of unary maps of read's bases.
    """
    bit_vector = ''
    for base in read:
        bit_vector += UNARY_MAP[base]
    
    return bit_vector
    

def simHash(bit_vector, ran_numbers):
    """
    Input: bit vector, sorted list of k random numbers.
    Output: similarity hash bin index.
    """
    ans = ''
    for num in ran_numbers:
        # prevent an error if num is out of range for bit_vector
        if num < len(bit_vector):
            ans += bit_vector[num]
    
    return ans
    

def hashLSH(reads, ran_numbers):
    """
    Build the similarity hash according to the given indices.
    
    Input: list of reads (ID, 'ACCTGCA...'), list of random integers.
    Output: reads' similarity hash map.
    """
    hashMap = {}
    for read in reads:
        index = simHash(unaryMap(read[1]), ran_numbers)
        if hashMap.has_key(index):
            hashMap[index].append(read)
        else:
            hashMap[index] = [read]
    
    return hashMap


def consensus(bucket):
    """
    Compute the consensus of a bucket of reads.
    
    Input: Bucket of raw reads as strings.
    Output: Majority vote read as string.
    """
    acgt = {0:'A', 1:'C', 2:'G', 3:'T'}
    sumVector = 4*len(bucket[0])*[0]
    average =''    
    
    for read in bucket:
        temp = unaryMap(read) # convert read to bit representation
        for idx in range(len(temp)):
            # compute base's frequency across reads, same position
            sumVector[idx] += int(temp[idx]) 

    best = 0
    index = 0
    count = 0
    
    # loop through the vector that recorded the frequencies of bases
    for idx in range(len(sumVector)):
        # find largest value in a set of four contiguous positions
        if sumVector[idx] > best:
            index = idx
            best = sumVector[idx]
        count += 1
        if count == 4:
            if best == 0:
                average += 'N' # if best value was zero
            else:
                average += acgt[index%4] # add the base with highest frequency
            
            best = 0
            index = 0
            count = 0
            
    return average
