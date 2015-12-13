#__author__ = 'soumadipmukherjee'
from __future__ import print_function
import random
import numpy as np
import math
import sys

translate_table = {
        'A': (1, 0),
        'T': (-1, 0),
        'C': (0, 1),
        'G': (0, -1)
    }

def reverseComplement(seq):
    seq_dict = {'A':'T','T':'A','G':'C','C':'G', 'N':'N'}
    return "".join([seq_dict[base] for base in reversed(seq)])

#@profile
def readFastq(filename, filename2, k_length):
    #sequence will be an array of read objects
    #the object will have the ID, the original sequence, and the reverse complement
    sequences = []
    fh = filename
    fh2 = filename2
    while True:
        read_obj = {
            "id": None,
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

        read_obj["id"] = id[:-1]

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
        read_obj["kmer_vectors"] = filter(lambda a: a != None, kmer_array)
        sequences.append(read_obj)
    return sequences


def generateNRandomVectors(n, k_length):
    """
    Generate n random kmers of length k_length. The matrix will actually
    be n x (2 * k_length).
    :return: the random elements in an array
    """
    #k_length *= 2
    listOfOptions = [(1, 0), (-1, 0), (0, 1),  (0, -1)]
    random_vectors = np.zeros((n, k_length), dtype=object)

    for x in range(0, n):
        for y in range(0, k_length): #, 2):
            random_vectors[x][y] = random.choice(listOfOptions)
            # random_vectors[x][y] = val[0]
            # random_vectors[x][y + 1] = val[1]
    return random_vectors


def produceKmerMatrix(kmers, k_length):
    kmerMatrix = np.zeros((len(kmers), k_length), dtype=object)
    for kindex, kmer in enumerate(kmers):
        for ntindex, val in enumerate(kmer):
            kmerMatrix[kindex][ntindex] = translate_table[val]
    return kmerMatrix


def multiplyVec(vec1, vec2):
    dotproduct = 0
    for i in xrange(len(vec1)):
        dotproduct += (vec2[i][0] * vec1[i][0]) + (vec2[i][1] * vec1[i][1])
    if dotproduct >= 0:
        return '1'
    else:
        return '0'

@profile
def produceAbundanceMatrix(reads, k_length, h_size):
    randomVectorMatrix = generateNRandomVectors(h_size, k_length)
#    print (randomVectorMatrix)
    abundanceMatrix = []
    #abundanceMatrix = np.zeros((len(reads),2**h_size), dtype=np.int32)
    #print randomVectorMatrix
    #print randomVectorMatrix.shape
    bitVector = []
    l_i = []
    g_j = []
    MAX_CONSTANT_32 = 4294967295.0
    read_index = 0
    #print ("reached here")
    for r_obj in reads:
        print ("\r\rReading another one... {}".format(read_index+1), end="", file=sys.stderr)
        #print (r_obj["kmer_vectors"])
        kmers = r_obj["kmer_vectors"]
        kmerMatrix = produceKmerMatrix(kmers, k_length)

        #print kmerMatrix.shape
        #RMatrix = np.dot(randomVectorMatrix, kmerMatrix)
        RMatrix = np.zeros((randomVectorMatrix.shape[0], kmerMatrix.shape[0]), dtype=np.str)

        for i in xrange(RMatrix.shape[0]):
            for j in xrange(RMatrix.shape[1]):
                RMatrix[i, j] = multiplyVec(randomVectorMatrix[i, :], kmerMatrix[j, :])

        for x in range(0, RMatrix.shape[1]):
            bitVector.append(''.join(RMatrix[:, x]))

        kmer_count = {}
        for index in bitVector:
            if(index in kmer_count):
                kmer_count[index] += 1
            else:
                kmer_count[index] = 1
        abundanceMatrix.append(kmer_count)

        read_index+=1
    print ("Done Reading DAMM FILES")
    return abundanceMatrix

    # print abundanceMatrix[0][0]
    # for i in range (0, len(reads)):
    #     for j in range(0, 2**h_size):
    #         print str(i), str(j)
    #         #print abundanceMatrix[i][j]
    #     print
    #     print
    # print "DONE"
    # for i in range(0,len(reads)):
    #     r_hits = []
    #     for j in range(0,2**h_size):
    #         if abundanceMatrix[i][j] != 0:
    #             r_hits.append(abundanceMatrix[i][j])
    #     #print sum(c_hits)
    #   l_i.append(math.sqrt(sum(r_hits))/(2**h_size))
    #
    #print l_i
    # print "DONE1"
    # for j in range(0,2**h_size):
    #     col_hit_count = 0
    #     for i in range(0,len(reads)):
    #         if abundanceMatrix[i][j] != 0:
    #             col_hit_count+=1
    #     if col_hit_count == 0:
    #         g_j.append(MAX_CONSTANT_32)
    #     else:
    #         g_j.append(math.log((len(reads)/col_hit_count),2))
    #         #print
    # #print g_j
    # print "DONE2"
    # #print l_i
    # for i in range(0, len(reads)):
    #     for j in range(0, 2**h_size):
    #         #print abundanceMatrix[i][j], g_j[j], l_i[i]
    #         if(abundanceMatrix[i][j] > 100):
    #             print "HERE"
    #             print g_j[j], l_i[i]
    #         abundanceMatrix[i][j]*= (g_j[j]/l_i[i])
    #
    # #print abundanceMatrix[0][readIndices[len(readIndices) -1]]
    # print "DONE3"
    #for i in range (0, len(reads)):
    #     for j in range(0, 2**h_size):
    #         print abundanceMatrix[i][j],
    #     print
    # return



def makeClusters(aMatrix, robj):
    clusters = {}
    for row_index, row in enumerate(aMatrix):
        binIndices = ''.join(sorted(row.keys()))
        if(binIndices in clusters.keys()):
            clusters[binIndices].append(row_index)
        else:
            clusters[binIndices] = [row_index]
    print (len(clusters))
    with open('cluster.txt', 'w') as filen:
        for _,v in clusters.iteritems():
            for readIndices in v:
                filen.write(robj[readIndices]["id"] + '\n')
            filen.write('\n')
