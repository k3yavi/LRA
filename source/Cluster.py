# -*- coding: utf-8 -*-
"""
Created on ........ Thu Nov 12 12:39:23 2015
Last modified on .. Thu Nov 19 23:18:39 2015

Cluster class and its helper functions for RNA/DNA reads clustering.

@author: Caleb Andrade
"""

# DIM = 1 manhattan metric; DIM = 2 euclidean metric
DIM = 1

#******************************************************************************
# HELPER FUNCTIONS
#******************************************************************************
def statsVectors(reads):
    """
    This function creates a list of stats vectors, one per read in 'reads'.
    Each stats vector is defined as:
    
    (#A's, #C's, #G's, #T's, d(A's), d(C's), d(G's), d(T's), #switches, errors)
    
        * By # A's we mean the total number of A's in the read.
        * By d(A's) we mean the total number of positions in between A's.
        * By #switches we mean the total number of base alternations in a read.
        * By errors we mean the number of undefined 'N' bases in a read.
        
    With this information we can encode a read as a 10-tuple, so that the
    manhattan/euclidean metrics can be used in a 10-dimensional space to
    measure closeness of the properties encoded in the stats vector of any 
    two reads.
    """
    stats_vectors = []
    for read in reads:
        # initialize read's counters of each base
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
        stats_vectors.append((count['A'], count['C'], count['G'], count['T'],
                        dist['A'][0], dist['C'][0], dist['G'][0], dist['T'][0],
                        num_switch - 1, errors))
    return stats_vectors
    

def distance(vector_a, vector_b, dim = DIM):
        """
        Compute Minkowski's distance between two vectors.
        Note: dim = 1, then manhattan distance; dim = 2 then euclidean.
        
        Input: two tuples with real values as entries.
        Output: Minkowski's distance between vector_a and vector_b.
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
    Class to model RNA/DNA reads clustering. The class itself does not 
    store reads, but rather tuple representations of them, namely, the 
    average of their stats vector. The stats vector description is found
    in the 'statsVectors' function's comment.
    """
    
    def __init__(self, reads = []):
        """
        Initializes a Cluster object for DNA/RNA reads.
        
        Input: a list of 'reads', empty by default. Each read in 'reads' 
        is assumed to be a tuple as follows: (ID, 'ACGTCAG...').
        Note: each ID should be unique.
        """
        self.sum_stats_vectors = 10*[0] # total sum of stats vectors
        self.avg_stats_vectors = 10*[0] # average of stats vectors
        self.reads_IDs = set([]) # IDs of reads in cluster
        self.stats_vectors = []
        self.weight = 0
                
        if reads != []:
            # create stats vectors and loop through them
            self.stats_vectors = statsVectors([read[1] for read in reads])
            for item in self.stats_vectors:
                # loop through every entry of a stats vector
                for idx in range(10):
                    # update cluster's sum stats vector
                    self.sum_stats_vectors[idx] += item[idx]
                self.weight += sum(item[0:4])
            # extract ID's
            self.reads_IDs = set([read[0] for read in reads])
            # update average stats vector
            self.update()
      
      
    def __str__(self):
        """
        As string.
        """
        string = "\nAverage stats vector: "
        temp = []
        for idx in range(10):
            temp.append(round(self.getAvgStats()[idx], 1))
        string += str(tuple(temp))
        string += "\nReads' IDs in cluster: "
        string += str(self.getIDs())
        string += "\nCluster error: " + str(round(self.clusterError(), DIM))
        
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
        
    
    def getWeight(self):
        """
        Get weight of cluster (sum of all bases of its readings).
        """
        return self.weight
        
        
    def getStatsVectors(self):
        """
        Get list of stats vectors.
        """
        return list(self.stats_vectors)
        
        
    def getAvgStats(self):
        """
        Get average stats vector.
        """
        return list(self.avg_stats_vectors)
        
        
    def getSumStats(self):
        """
        Get the sum of all stats vectors.
        """
        return list(self.sum_stats_vectors)
        
        
    def distance(self, other_cluster, dim = DIM):
        """
        Compute Minkowski's distance between self's and other cluster's 
        avg_stats_vectors.
        """
        return distance(self.getAvgStats(), other_cluster.getAvgStats(), dim)
        
    
    def update(self):
        """
        Update average stats vector and weight.
        """
        for idx in range(10):
            avg = float(self.sum_stats_vectors[idx]) 
            if self.getSize() != 0:
                avg = avg / self.getSize()
            self.avg_stats_vectors[idx] = avg         
        
        self.weight = sum(self.sum_stats_vectors[:4])
          

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
        # The set of stats vector of self merges with other's
        self.stats_vectors += other_cluster.getStatsVectors()
        # The sum_stats_vectors vectors' of both clusters are added.
        for idx in range(10):
            self.sum_stats_vectors[idx] += other_cluster.sum_stats_vectors[idx]
        # The avg_stats_vectors vector and weight are updated.       
        self.update()
                        
        
    def copy(self):
        """
        Make an exact copy of self.
        """
        copy = Cluster()
        copy.sum_stats_vectors = self.getSumStats()
        copy.avg_stats_vectors = self.getAvgStats()
        copy.reads_IDs = self.getIDs()
        copy.stats_vectors = list(self.getStatsVectors())
        
        return copy       
        
        
    def isEqualTo(self, other_cluster):
        """
        Is self equal to other_cluster?
        """
        if type(other_cluster) != type(self):
            raise Exception('Impossible comparison: different object types!')
        if self.getSize() != other_cluster.getSize():
            return False
        if self.getSumStats() != other_cluster.getSumStats():
            return False
        if self.getIDs() != other_cluster.getIDs():
            return False
        if self.getStatsVectors() != other_cluster.getStatsVectors():
            return False
        
        return True
        
    def clusterError(self, dim = DIM):
        """
        Standard deviation of stats vectors from average stats vector.
        """
        std_dev = 0
        if self.getSize() == 0:
            return 0
        for vector in self.getStatsVectors():
            std_dev += distance(vector, self.getAvgStats(), dim)**2
        
        return (float(std_dev) / self.getSize())**0.5
        
        
