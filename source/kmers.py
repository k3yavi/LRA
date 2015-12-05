#!/user/bin/python
import numpy
import itertools
import math
import pyfasta
__author__ = 'soumadipmukherjee'
def performKmerHash(klength,reads):
    k_hash = {}
    #r[0] is the ID of the read
    #r[1] is the actual full read
    #for each read in the reads array we build the kmer hash
    for r in reads:
        full_read = r[1] #assign the full read to the variable
        u_limit = len(full_read) - klength + 1 #the upper limit below loop
        r.append(set())
        for x in range(0, u_limit):
            #print full_read[x:(x+klength)]
            r[2].add(full_read[x:(x+klength)])
            try:
                #if the key is present append the read id to the SET
                k_hash[full_read[x:(x+klength)]].add(r[0])
            except KeyError:
                #if key isnt present, create a new entry
                k_hash[full_read[x:(x+klength)]] = set([r[0]])
        #inner for
    #outer for
    #print k_hash
    #print reads
    return k_hash, reads

def cosineDistance(reads):
    matrix_const = len(reads)
    cos_matrix = numpy.zeros((matrix_const, matrix_const))#numpy.empty((matrix_const, matrix_const)) #constructy 2D array
    #cos_matrix[:] = numpy.Infinity    #initialize array with infinity
    #we will have a read by read length array. The bottom left and the diagnol are all Infinity
    #the closer the number is to ONE, the closer the points are in K-Dimensional space
    #the closer the number is to ZERO, the more distance there is K-Dimensional space
    for x in range(0,matrix_const):
        for y in range((x+1),matrix_const):
            kmers_x = reads[x][2] #get set of kmers in read_x
            kmers_y =  reads[y][2] #get set of kmers in read_y
            intersection =  len(kmers_x.intersection(kmers_y))#find the number of intersections
            numerator = intersection
            denom_x = len(kmers_x) #find length of kmers_x, this will be used for magnitude
            denom_y = len(kmers_y) #find length of kmers_y, this will be used for magnitude
            cosTheta = numerator / (math.sqrt(denom_x) * math.sqrt(denom_y))
            cos_matrix[x][y] = cosTheta
    print cos_matrix
    return cos_matrix


def hierarchicalClustering(k_hash, r_hash, timesCluster):
    matrix_const = len(r_hash)
    for t in range(0, timesCluster):
        matrix = numpy.zeros((matrix_const, matrix_const))
        #print matrix
        for i in k_hash:
            #print i, k_hash[i]
            nCr =  list(itertools.combinations(k_hash[i], 2))
            #print nCr
            for p in nCr:
                #print p
                x = p[0]
                y = p[1]
                try:
                    matrix[x][y] += 1
                except IndexError:
                    print "INDEX ERROR"
                    print p
        if matrix.argmax() == 0:
            print "No more clustering can be performed"
            return # if the there is nothing left to combine then return
        max_index =  numpy.unravel_index(matrix.argmax(), matrix.shape)

        print "Hierarchically combining : ", max_index
        if max_index[0] < max_index[1]:
            get_read = r_hash[max_index[1]]
            #del r_hash[max_index[1]]
            r_hash[max_index[1]] = "DNE"
            for e in get_read:
                r_hash[max_index[0]].add(e)
            for i in k_hash:
                if max_index[1] in k_hash[i]:
                    k_hash[i].remove(max_index[1])
        else:
            get_read = r_hash[max_index[1]]
            #del r_hash[max_index[0]]
            r_hash[max_index[0]] = "DNE"
            for e in get_read:
                r_hash[max_index[1]].add(e)
            for i in k_hash:
                if max_index[0] in k_hash[i]:
                    k_hash[i].remove(max_index[0])
    #print nCr
    #print k_hash


def main():
    #reads have a unique ID and the read itself as coming in a tuple
    reads = [[0, "ATATCGCGAT"], [1, "CGCGATCG"], [2, "CTCGATCGATCAGACAGTACAGTACA"] ]
    #initially stores the read
    r_hash = []#{}
    for r in reads:
        r_hash.append(set([r[0]]))

    print r_hash
    k_hash, reads = performKmerHash(2, reads)
    print k_hash
    print reads
    cosineDistance(reads)
    hierarchicalClustering(k_hash, r_hash, 1)
    clusters = []
    for r in r_hash:
        if r != "DNE":
            clusters.append(r)
    for i in range(0, len(clusters)):
        print "Cluster ", str(i+1), "-> Read IDs: ", ", ".join(str(e) for e in clusters[i])#clusters[i]
    return

main()