# -*- coding: utf-8 -*-
"""
Created on ........ Thu Nov 12 12:39:23 2015
Last modified on .. Sun Dec 13 17:57:46 2015

Cluster class and its helper functions for RNA/DNA reads clustering.

@author: Caleb Andrade
"""

# DIM = 1 manhattan metric; DIM = 2 euclidean metric
DIM = 1

#******************************************************************************
# HELPER FUNCTIONS
#******************************************************************************
def codonAbundance(reads):
    """
    Create a list of codon abundance vectors of reads.
    """
    # codon abundance vectors list    
    abundance_vectors = []
    
    # disambiguation hash for lower/upper case bases
    base_map = {'a':'A', 'c':'C', 'g':'G', 't':'T', 'n':'N',
                'A':'A', 'C':'C', 'G':'G', 'T':'T', 'N':'N'}
    
    # codon hash, map codons to abundance codon vector indices
    codon_map = {'AAA':0, 'AAC':1, 'AAG':2, 'AAT':3,
             'ACA':4, 'ACC':5, 'ACG':6, 'ACT':7,
             'AGA':8, 'AGC':9, 'AGG':10, 'AGT':11,
             'ATA':12, 'ATC':13, 'ATG':14, 'ATT':15,
             'CAA':16, 'CAC':17, 'CAG':18, 'CAT':19,
             'CCA':20, 'CCC':21, 'CCG':22, 'CCT':23,
             'CGA':24, 'CGC':25, 'CGG':26, 'CGT':27,
             'CTA':28, 'CTC':29, 'CTG':30, 'CTT':31,
             'GAA':32, 'GAC':33, 'GAG':34, 'GAT':35,
             'GCA':36, 'GCC':37, 'GCG':38, 'GCT':39,
             'GGA':40, 'GGC':41, 'GGG':42, 'GGT':43,
             'GTA':44, 'GTC':45, 'GTG':46, 'GTT':47,
             'TAA':48, 'TAC':49, 'TAG':50, 'TAT':51,
             'TCA':52, 'TCC':53, 'TCG':54, 'TCT':55,
             'TGA':56, 'TGC':57, 'TGG':58, 'TGT':59,
             'TTA':60, 'TTC':61, 'TTG':62, 'TTT':63}
             
    for read in reads:
        codon_vector = 64*[0]
        for i in range(1 + len(read) - 3):
            cod0 = base_map[read[i]]
            cod1 = base_map[read[i+1]]
            cod2 = base_map[read[i+2]]
            # ignore a codon with 'N'
            if cod0 != 'N' and cod1 != 'N' and cod2 != 'N':
                cod = cod0 + cod1 + cod2
                codon_vector[codon_map[cod]] += 1
                
        abundance_vectors.append(codon_vector)
        
    return abundance_vectors
    

def distance(vector_a, vector_b, dim = DIM):
        """
        Compute Minkowski's distance between two vectors.
        Note: dim = 1, then manhattan distance; dim = 2 then euclidean.
        """
        distance = 0
        
        for idx in range(len(vector_a)):
            distance += abs(vector_a[idx] - vector_b[idx])**dim
        
        if dim == 1:
            return distance
        
        return distance**(1 / float(dim))

#******************************************************************************
# CLUSTER CLASS
#******************************************************************************
class Cluster(object):
    """
    Class to model RNA/DNA reads clusters. The object itself does not 
    store reads, but rather, an euclidean embedding, namely, the average
    of their codon abundance vectors. 
    """
    
    def __init__(self, reads = []):
        """
        Input: a list of reads, empty by default. Each read is assumed to 
        be a tuple as follows: (ID, 'ACGTCAG...'). Each ID should be unique.
        """
        self.sum_abundance_vectors = 64*[0] # coordinate-wise sum of abundance vectors
        self.avg_abundance_vectors = 64*[0] # average of abundance vectors
        self.reads_IDs = set([]) # IDs of reads in cluster
        self.abundance_vectors = []
                
        if reads != []:
            # create abundance vectors and loop them
            self.abundance_vectors = codonAbundance([read[1] for read in reads])
            for i in range(len(self.abundance_vectors)):
                item = self.abundance_vectors[i]
                # loop through every entry of an abundance vector
                for idx in range(64):
                    # update cluster's sum abundance vector
                    self.sum_abundance_vectors[idx] += item[idx]
            # extract ID's
            self.reads_IDs = set([read[0] for read in reads]) # save ID's
            # update average abundance vector
            self.update()
      
      
    def __str__(self):
        """
        As string.
        """
        string = "\nAverage abundance vector: \n"
        temp = []
        for idx in range(64):
            temp.append(round(self.getAvgAbundance()[idx], 1))
        string += str(tuple(temp))
        string += "\nReads' IDs in cluster: \n"
        string += str(self.getIDs())
        string += "\nCluster error: \n" + str(round(self.clusterError(), DIM))
        
        return string
        
            
    def getIDs(self):
        """
        Get set of reads' IDs.
        """
        return set(self.reads_IDs)
        
        
    def getSize(self):
        """
        Get size of cluster.
        """
        return len(self.reads_IDs)
        
    
    def getCodonAbundance(self):
        """
        Get list of abundance vectors.
        """
        return list(self.abundance_vectors)
        
        
    def getAvgAbundance(self):
        """
        Get average abundance vector.
        """
        return list(self.avg_abundance_vectors)
        
        
    def getSumAbundance(self):
        """
        Get the sum of all abundance vectors.
        """
        return list(self.sum_abundance_vectors)
        
        
    def distance(self, other_cluster, dim = DIM):
        """
        Compute Minkowski's distance between self's and other cluster's 
        avg_abundance_vectors.
        """
        return distance(self.getAvgAbundance(), other_cluster.getAvgAbundance(), dim)
        
    
    def update(self):
        """
        Update average abundance vector.
        """
        for idx in range(64):
            avg = float(self.sum_abundance_vectors[idx]) 
            if self.getSize() != 0:
                avg = avg / self.getSize()
            self.avg_abundance_vectors[idx] = avg         
    

    def mergeClusters(self, other_cluster):
        """
        Mutate self into a new cluster by merging it with other_cluster.
        """
        if self.isEqualTo(other_cluster):
            return
        # The sets of both clusters' ID's are merged.
        if self.getIDs().intersection(other_cluster.getIDs()) != set([]):
            raise Warning('Some reads had the same ID in clusters to merge')
        self.reads_IDs = self.reads_IDs.union(other_cluster.getIDs())
        # The set of abundance vector of self merges with other's
        self.abundance_vectors += other_cluster.getCodonAbundance()
        # The sum_abundance_vectors vectors' of both clusters are added.
        for idx in range(64):
            self.sum_abundance_vectors[idx] += other_cluster.sum_abundance_vectors[idx]
        # The avg_abundance_vectors vector is updated.       
        self.update()
                        
        
    def copy(self):
        """
        Make an exact copy of self.
        """
        copy = Cluster()
        copy.sum_abundance_vectors = self.getSumAbundance()
        copy.avg_abundance_vectors = self.getAvgAbundance()
        copy.reads_IDs = self.getIDs()
        copy.abundance_vectors = list(self.getCodonAbundance())
        
        return copy       
        
        
    def isEqualTo(self, other_cluster):
        """
        Is self equal to other_cluster?
        """
        if type(other_cluster) != type(self):
            raise Exception('Different object types!')
        if self.getSize() != other_cluster.getSize():
            return False
        if self.getSumAbundance() != other_cluster.getSumAbundance():
            return False
        if self.getIDs() != other_cluster.getIDs():
            return False
        if self.getCodonAbundance() != other_cluster.getCodonAbundance():
            return False
        
        return True
        
   
    def clusterError(self, dim = DIM):
        """
        Standard deviation of abundance vectors from average abundance vector.
        """
        std_dev = 0
        if self.getSize() == 0:
            return 0
        for vector in self.getCodonAbundance():
            std_dev += distance(vector, self.getAvgAbundance(), dim)**2
        
        return (float(std_dev) / self.getSize())**0.5
        