#__author__ = 'soumadipmukherjee'
from __future__ import print_function
import random
import numpy as np
from collections import defaultdict
import sys

translate_table = {
        'A': (1, 0),
        'T': (-1, 0),
        'C': (0, 1),
        'G': (0, -1)
}


def reverseComplement(seq):
    seq_dict = {'A':'T','T':'A','G':'C','C':'G', 'N':'N'}
    return ''.join([seq_dict[base] for base in reversed(seq)])

#@profile
def readFastq(filename, filename2, k_length):
    #sequence will be an array of read objects
    #the object will have the ID, the original sequence, and the reverse complement
    sequences = []
    fh = filename
    fh2 = filename2
    while True:
        read_obj = {
            'id': None,
        }
        #Do all of the reading here. READ FOUR lines for each READ
        id = fh.readline() #read the id
        id_2= fh2.readline().rsplit()
        seq = fh.readline().rstrip() #read base sequence
        seq_2 = fh2.readline().rstrip()
        fh.readline() # skip placeholder line
        fh.readline() # base quality line
        fh2.readline()
        fh2.readline()
        #End of Reading here

        if len(seq) == 0 or len(seq_2)== 0:
            break
        #Make the mates upper case
        seq = seq.upper()
        seq_2 = seq_2.upper()

        read_obj['id'] = id[:-1]

        mate1 = seq + reverseComplement(seq_2)
        mate2 = reverseComplement(seq) + seq_2
        skipRange = xrange(len(seq)-k_length+1, len(seq))
        kmer_array = [None] * (len(mate1)-k_length+1)
        for index in xrange(len(mate1)-k_length+1):
            k1 = mate1[index:index+k_length]
            k2 = mate2[index:index+k_length]

            if index in skipRange:
                continue

            if 'N' in k1 and 'N' in k2:
                continue
            elif 'N' in k1:
                kmer_array[index] = k2
            elif 'N' in k2:
                kmer_array[index] = k1
            elif k1 <= k2:
                kmer_array[index] = k1
            else:
                kmer_array[index] = k2
        read_obj['kmer_vectors'] = filter(lambda a: a != None, kmer_array)
        sequences.append(read_obj)
    return sequences


def generateNRandomVectors(n, k_length):
    """
    Generate n random kmers of length k_length.
    The matrices will be n x k_length.
    :return: tuple of 2 arrays.
    """
    listOfOptions = [(1, 0), (-1, 0), (0, 1),  (0, -1)]
    random_kmers_x = np.zeros((n, k_length), dtype=np.int8)
    random_kmers_y = np.zeros((n, k_length), dtype=np.int8)

    for x in range(0, n):
        for y in range(0, k_length):
            encoded_base = random.choice(listOfOptions)
            # Hold the random base's x into the first matrix
            random_kmers_x[x][y] = encoded_base[0]
            # And the random base's y into the second matrix
            random_kmers_y[x][y] = encoded_base[1]

    return random_kmers_x, random_kmers_y


def produceKmerMatrix(kmers, k_length):
    """
    Generate kmer matrices in base-encoded form.
    The matrices will be k_length x len(kmers).
    :return: tuple of 2 arrays.
    """
    kmer_matrix_x = np.zeros((k_length, len(kmers)), dtype=np.int8)
    kmer_matrix_y = np.zeros((k_length, len(kmers)), dtype=np.int8)
    for kindex, kmer in enumerate(kmers):
        for ntindex, base in enumerate(kmer):
            encoded_base = translate_table[base]
            kmer_matrix_x[ntindex][kindex] = encoded_base[0]
            kmer_matrix_y[ntindex][kindex] = encoded_base[1]

    return kmer_matrix_x, kmer_matrix_y


def multiplyVec(vec1, vec2):
    dotproduct = 0
    for i in xrange(len(vec1)):
        dotproduct += (vec2[i][0] * vec1[i][0]) + (vec2[i][1] * vec1[i][1])
    if dotproduct >= 0:
        return '1'
    else:
        return '0'


#@profile
def produceAbundanceMatrix(reads, k_length, h_size):
    rand_kmers_x, rand_kmers_y = generateNRandomVectors(h_size, k_length)
    abundanceMatrix = []
    read_index = 0
    for r_obj in reads:
        print ('\r\rReading another one... {}'.format(read_index+1), end='', file=sys.stderr)
        kmers = r_obj['kmer_vectors']
        kmers_x, kmers_y = produceKmerMatrix(kmers, k_length)

        # X coordinates of encoded bases
        XMatrix = np.dot(rand_kmers_x, kmers_x)
        # Y coordinates of encoded bases
        YMatrix = np.dot(rand_kmers_y, kmers_y)
        # Sum of X,Y
        RMatrix = np.add(XMatrix, YMatrix)

        # Convert each column into an integer and count
        kmer_count = defaultdict(int)
        for x in range(0, RMatrix.shape[1]):
            col = RMatrix[:, x]
            bit_val = 0
            for index, product in enumerate(col):
                if product >= 0:
                    bit_val |= 1 << index

            kmer_count[bit_val] += 1

        abundanceMatrix.append(kmer_count)
        read_index += 1

    print ('\nDone hashing files.')
    return abundanceMatrix


def makeClusters(aMatrix, robj):
    clusters = {}
    for row_index, row in enumerate(aMatrix):
        binIndices = tuple(sorted(row.keys()))
        if binIndices in clusters:
            clusters[binIndices].append(row_index)
        else:
            clusters[binIndices] = [row_index]
    print (len(clusters))
    with open('cluster.txt', 'w') as filen:
        for _, v in clusters.iteritems():
            for readIndices in v:
                filen.write(robj[readIndices]['id'] + '\n')
            filen.write('\n')
