#__author__ = 'soumadipmukherjee'
from __future__ import print_function
import random
import numpy as np
import math
import sys
import gc

translate_table = {
        'A': (1, 0),
        'T': (-1, 0),
        'C': (0, 1),
        'G': (0, -1)
    }

#function returns the vector representation of the kmer, according to our mapping
#@profile
#def translateKmer(kmer):
#    translate_table = {
#        'A': 1,
#        'T': -1,
#        'C': 2,
#        'G': -2
#    }
#    vector = [None] * len(kmer)
#    index = 0
#    for k in kmer:
#        vector[index] = translate_table[k]
#        index += 1
#    return vector


def reverseComplement(seq):
    seq_dict = {'A':'T','T':'A','G':'C','C':'G', 'N':'N'}
    return "".join([seq_dict[base] for base in reversed(seq)])

#@profile
def readFastq(filename, filename2, k_length):
    """
    Parse read and quality strings from a FASTQ file with sequencing reads.
    @author: Ben Langmead & Jacob Pritt.
    Input: file path
    Output: A list of reads, the list of qualities
    """
    #sequence will be an array of read objects
    #the object will have the ID, the original sequence, and the reverse complement
    sequences = []
#    compMapping = {
#        'A' : 'T',
#        'a' : 'T',
#        'T' : 'A',
#        't' : 'A',
#        'C' : 'G',
#        'c' : 'G',
#        'G' : 'C',
#        'g' : 'C',
#        'N' :'N'
#        'n': 'N'
#    }
    fh = filename
    fh2 = filename2
    #with open(filename) as fh:
    #    with open(filename2) as fh2:
    while True:
        read_obj = {
            "id": None,
            #"seq": None,
            #"rev_comp_seq": None
            #"mate1":None,
            #"mate2":None
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

        #Store the mate ids in the structure to retun
        #read_obj["mate1_id"] = id[0] # skip name line
        #read_obj["mate2_id"] = id_2[0]

        #produce the id of the Reads by removing the las part
        #end =read_obj["mate1_id"].rindex(':')#get the last index of the colon
        read_obj["id"] = id[:-1]

        #store the length of mate1 so we do kmer has over the bridge of mate1 and mate2
        #read_obj["mate1_len"] = len(seq)


        #construct Read1
        #Reverse complement mate2 and concatenat it to mate1
#        reverse_2 = seq_2[::-1]
#        reverse_comp_r2=""
#        for x in range(0, len(reverse_2)):
#            try:
#                reverse_comp_r2 += compMapping[reverse_2[x]]
#            except KeyError as e:
#                #print e
#                continue

        #read_obj["read1"] = mate1


        #construct Read2
        #Reverse complement mate1 and concatenate it with mate2
#        reverse_1 = seq[::-1]
#        reverse_comp_r1= ""
#        for x in range(0, len(reverse_1)):
#            try:
#                reverse_comp_r1 += compMapping[reverse_1[x]]
#            except KeyError as e:
#                #print e
#                continue

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
        kmer_array = filter(lambda a: a != None, kmer_array)
        #read_obj["read2"] = mate2

        #bridgeIndex = len(mate1) - 1
        #kmer_array =  []
        ##go from 0 to bridgeIndex
        ##Dont consider Kmers with 'N' base in them
        #for x in range(0,((bridgeIndex+1)-k_length)):
        #    k1 = mate1[x:x+k_length]
        #    k2 = mate2[x:x+k_length]
        #    b1 = 'N' in k1
        #    b2 = 'N' in k2
        #    if b1 == True and b2 == True:
        #        continue
        #    elif b1 == True and b2 == False:
        #        kmer_array.append(k2)
        #    elif b1 == False and b2 == True:
        #        kmer_array.append(k1)
        #    elif b1 <= b2:
        #        kmer_array.append(k1)
        #    else:
        #        kmer_array.append(k2)
        ##from bridge index to the end of the read
        #for x in range(bridgeIndex,((len(mate1))-k_length)):
        #    k1 = mate1[x:x+k_length]
        #    k2 = mate2[x:x+k_length]
        #    b1 = 'N' in k1
        #    b2 = 'N' in k2
        #    if b1 == True and b2 == True:
        #        continue
        #    elif b1 == True and b2 == False:
        #        kmer_array.append(k2)
        #    elif b1 == False and b2 == True:
        #        kmer_array.append(k1)
        #    elif b1 <= b2:
        #        kmer_array.append(k1)
        #    else:
        #        kmer_array.append(k2)
        #r_obj["kmer_list"] = kmer_array
        #kmer_vectors = []
        #for k in kmer_array:
        #    kmer_vectors.append(translateKmer(k))

        read_obj["kmer_vectors"] = kmer_array
        sequences.append(read_obj)
    return sequences
#def buildKmerListForReads(reads, k_length):
#    for r_obj in reads:
#        r1 = r_obj["read1"]
#        r2 = r_obj["read2"]
#        bridgeIndex = r_obj["mate1_len"] - 1
#        kmer_array =  []
#        #go from 0 to bridgeIndex
#        #Dont consider Kmers with 'N' base in them
#        for x in range(0,((bridgeIndex+1)-k_length)):
#            k1 = r1[x:x+k_length]
#            k2 = r2[x:x+k_length]
#            b1 = 'N' in k1
#            b2 = 'N' in k2
#            if b1 == True and b2 == True:
#                continue
#            elif b1 == True and b2 == False:
#                kmer_array.append(k2)
#            elif b1 == False and b2 == True:
#                kmer_array.append(k1)
#            elif b1 <= b2:
#                kmer_array.append(k1)
#            else:
#                kmer_array.append(k2)
#        #from bridge index to the end of the read
#        for x in range(bridgeIndex,((len(r1))-k_length)):
#            k1 = r1[x:x+k_length]
#            k2 = r2[x:x+k_length]
#            b1 = 'N' in k1
#            b2 = 'N' in k2
#            if b1 == True and b2 == True:
#                continue
#            elif b1 == True and b2 == False:
#                kmer_array.append(k2)
#            elif b1 == False and b2 == True:
#                kmer_array.append(k1)
#            elif b1 <= b2:
#                kmer_array.append(k1)
#            else:
#                kmer_array.append(k2)
#        r_obj["kmer_list"] = kmer_array
#        r_obj["read1"] = ""
#        r_obj["read2"] = ""
#    return


#transform the list of kmers to vectors representations
#def translateKmerList(reads):
#    for r_obj in reads:
#        kmer_list = r_obj["kmer_list"]
#        kmer_vectors = []
#        for k in kmer_list:
#                kmer_vectors.append(translateKmer(k))
#
#        r_obj["kmer_vectors"] = kmer_vectors
#    return

#function will generate N random kmers of length k_length and return the random elements
 #in and array
def generateNRandomVectors(n, k_length):
    listOfOptions = [(1, 0), (-1, 0), (0, 1),  (0, -1)]
    random_vectors = np.zeros((n, k_length), dtype=object)

    for x in range(0, n):
        for y in range(0,k_length):
            random_vectors[x][y] = random.choice(listOfOptions)

    return random_vectors

#def produceRandomMatrix(randomVectors, k_length):
#    randomVectorsMatrix = np.zeros((len(randomVectors), k_length), dtype= np.int8)
#    i = 0;
#    j = 0;
#    for r in randomVectors:
#        for val in r:
#            randomVectorsMatrix[i][j] = val
#            j+=1
#
#        j = 0
#        i +=1
#    return randomVectorsMatrix

def produceKmerMatrix(kmers, k_length):
    kmerMatrix = np.zeros((k_length, len(kmers)), dtype=object)
    for kindex, kmer in enumerate(kmers):
        for ntindex, val in enumerate(kmer):
            kmerMatrix[ntindex][kindex] = translate_table[val]

    return kmerMatrix


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
        RMatrix = np.zeros((randomVectorMatrix.shape[0], kmerMatrix.shape[1]), dtype=np.str)

        for i in xrange(RMatrix.shape[0]):
            for j in xrange(RMatrix.shape[1]):
                RMatrix[i, j] = multiplyVec(randomVectorMatrix[i, :], kmerMatrix[:, j].transpose())
        #RMatrix = np.dot(randomVectorMatrix , kmerMatrix)
        #RMatrix = [[sum(a*b for a,b in zip(X_row,Y_col)) for Y_col in zip(*kmerMatrix)] for X_row in randomVectorMatrix]
        # shape  = RMatrix.shape
        #print shape
        for x in range(0, RMatrix.shape[1]):
            #col = RMatrix[:,x]
            ##normalized = math.sqrt(sum([y ** 2  for y in col]))
            ###print normalized
            #normalized = np.sqrt(np.sum(np.square(col), axis=0))
            #cosine = col / normalized

            #print "Normalized", str(normalized)
            #print cosine
            #newCol = ''
            #for c in cosine:
            #    if c < -1 or c > 1:
            #        print ("error")
            #    theta = np.arccos(c)
            #    # print theta
            #    if theta < (math.pi/2):
            #        newCol += '1'
            #    else:
            #        newCol += '0'
            bitVector.append(''.join(RMatrix[:, x]))

        kmer_count = {}
        for index in bitVector:
            if(index in kmer_count):
                kmer_count[index] += 1
            else:
                kmer_count[index] = 1
        #kmer_list = []
        #for k,v in kmer_count.iteritems():
        #    kmer_list.append((k,v))
        abundanceMatrix.append(kmer_count)

        #abundanceMatrix[read_index][index] += 1

        read_index+=1
    print ("Done Reading DAMM FILES")
    #print (abundanceMatrix)
    return abundanceMatrix
    # file.write(str(abundanceMatrix))
    # file.close()

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



def makeClusters(aMatrix):
    clusters = {}
    for row_index, row in enumerate(aMatrix):
        binIndices = ''.join(sorted(row.keys()))
        if(binIndices in clusters.keys()):
            clusters[binIndices].append(row_index)
        else:
            clusters[binIndices] = [row_index]
    print (len(clusters))
#    for _,v in clusters.iteritems():
#        print (v, end= "\n")
#    print (clusters)
#    print (len(clusters))
#    return clusters








#def prettyPrintObject(sequences):
#    for x in sequences:
#        print "#" * 100
#        print "ID: ", x["id"]
#        print "ID1: ", x["mate1_id"]
#        print "ID2: ", x["mate2_id"]
#        print "Read1: ", x["read1"]
#        print "Read2: ", x["read2"]
#        print "Mate1 Length: ", x["mate1_len"]
#        # print "S : ", x["seq"]
#        # print "RC: ", x["rev_comp_seq"]
#        print "#" * 100
#        print ""
#
#prettyPrintObject(readFastq("dummyTest.fastq"))

#rettyPrintObject(readFastq("r1_short.fq", "r2_short.fq"))

# r = readFastq("r1_short.fq", "r2_short.fq")
# print "Finished reading files"
# #print r
# # prettyPrintObject(r)
# buildKmerListForReads(r, 33)
# print "Finished building kmwe List for each read"
#
# # #print r
# #
# translateKmerList(r)
# print "Finished Translating kmer list"
# # print r
#
# #randomVecs = generateNRandomVectors(10, 4)
#
# #produceRandomMatrix(randomVecs, 4)
# produceAbundanceMatrix(r, 33, 29)


#print translateKmer("ATCG")
