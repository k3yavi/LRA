__author__ = 'soumadipmukherjee'
def readFastq(filename):
    """
    Parse read and quality strings from a FASTQ file with sequencing reads.
    @author: Ben Langmead & Jacob Pritt.

    Input: file path
    Output: A list of reads, the list of qualities
    """

    #sequence will be an array of read objects
    #the object will have the ID, the original sequence, and the reverse complement
    sequences = []
    #qualities = []
    compMapping = {
        'A' : 'T',
        'a' : 'T',
        'T' : 'A',
        't' : 'A',
        'C' : 'G',
        'c' : 'G',
        'G' : 'C',
        'g' : 'C',
        'n' : 'N',
        'N' : 'N'
    }

    with open(filename) as fh:
        while True:
            read_obj = {
                "id": None,
                "seq": None,
                "rev_comp_seq": None
            }
            id = fh.readline().rsplit() #read the id
            seq = fh.readline().rstrip() #read base sequence
            fh.readline() # skip placeholder line
            fh.readline() # base quality line
            if len(seq) == 0:
                break
            read_obj["id"] = id[0] # skip name line
            read_obj["seq"] = seq.upper()
            reversed = seq[::-1]
            #print reversed
            rev_comp = ""

            for x in range(0, len(reversed)):
                try:
                    #print reversed[x], compMapping[reversed[x]]
                    rev_comp += compMapping[reversed[x]]
                except KeyError as e:
                    print e
                #rev_comp += compMapping[reversed[x]]

            read_obj["rev_comp_seq"] = rev_comp




            sequences.append(read_obj)

    return sequences#, qualities

def prettyPrintObject(sequences):
    for x in sequences:
        print "#" * 100
        print "ID: ", x["id"]
        print "S : ", x["seq"]
        print "RC: ", x["rev_comp_seq"]
        print "#" * 100
        print ""
prettyPrintObject(readFastq("dummyTest.fastq"))