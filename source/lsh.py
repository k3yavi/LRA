#!/usr/bin/env python

import numpy
from collections import namedtuple
from lsh_functions import *

"""
Locality Sensitive Hash (LSH).
"""

# All defaults are stored here
Defaults = namedtuple('Defaults', ['kmer_size', 'hash_size'])
defaults = Defaults(kmer_size=35, hash_size=20)


def main():
    args = parse_args()
    print "HERE"
    # Hash input files.
    matrix = loc_hash(args.kmer_size, args.hash_size, args.file1, args.file2)
    # Pass matrix to gensim.models.LsiModel


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Locality Sensitive Hash (LSH)')
    parser.add_argument('-k', '--kmer-size', action='store', dest='kmer_size',
                        type=int, default=defaults.kmer_size,
                        help='kmer length, default: {}'.format(defaults.kmer_size))
    parser.add_argument('-s', '--hash-size', action='store', dest='hash_size',
                        type=int, default=defaults.hash_size,
                        help='Hash size, default: {}'.format(defaults.hash_size))
    parser.add_argument('file1', action='store',
                        type=argparse.FileType('r'),
                        help='FASTQ file 1 of paired end read')
    parser.add_argument('file2', action='store',
                        type=argparse.FileType('r'),
                        help='FASTQ file 2 of paired end read')
    return parser.parse_args()


def loc_hash(kmer_size, hash_size, file1, file2):
    """
    Performs Locality Sensitive Hashing.
    :param kmer_size: length of kmers.
    :param hash_size: length of hash?
    :param file1: FASTQ file 1 of paired end reads.
    :param file2: FASTQ file 2 of paired end reads.
    :return: Numpy matrix/gensim corpus to pass to gensim?
    """
    print file1
    print file2
    reads = readFastq(file1, file2)
    buildKmerListForReads(reads, kmer_size)
    translateKmerList(reads)
    abundanceMatrix = produceAbundanceMatrix(reads, kmer_size, hash_size)
    clusters = makeClusters(abundanceMatrix)
    print len(clusters)
    print clusters
    #print reads
    #pass


if __name__ == '__main__':
    main()
