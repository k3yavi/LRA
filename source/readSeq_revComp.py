__author__ = 'soumadipmukherjee'
def readFastq(filename, filename2):
    """
    Parse read and quality strings from a FASTQ file with sequencing reads.
    @author: Ben Langmead & Jacob Pritt.
    Input: file path
    Output: A list of reads, the list of qualities
    """
    #sequence will be an array of read objects
    #the object will have the ID, the original sequence, and the reverse complement
    sequences = []
    compMapping = {
        'A' : 'T',
        'a' : 'T',
        'T' : 'A',
        't' : 'A',
        'C' : 'G',
        'c' : 'G',
        'G' : 'C',
        'g' : 'C'
    }
    with open(filename) as fh:
        with open(filename2) as fh2:
            while True:
                read_obj = {
                    "id": None,
                    #"seq": None,
                    #"rev_comp_seq": None
                    "mate1":None,
                    "mate2":None
                }
                id = fh.readline().rsplit() #read the id
                id_2= fh2.readline().rsplit()
                seq = fh.readline().rstrip() #read base sequence
                seq_2 = fh2.readline().rstrip()
                fh.readline() # skip placeholder line
                fh.readline() # base quality line
                fh2.readline()
                fh2.readline()

                if len(seq) == 0 or len(seq_2)== 0:
                    break

                seq = seq.upper()
                seq_2 = seq_2.upper()
                read_obj["id1"] = id[0] # skip name line
                read_obj["id2"] = id_2[0]
                end =read_obj["id1"].rindex(':')
                print "END: ", end
                read_obj["id"] = read_obj["id1"][:end]
                read_obj["r1_len"] = len(seq)
                #read_obj["seq"] = seq.upper()

                #construct Mate1
                #Reverse complement read2 and concatenat it to read1
                reverse_2 = seq_2[::-1]
                reverse_comp_r2=""
                for x in range(0, len(reverse_2)):
                    try:
                        reverse_comp_r2 += compMapping[reverse_2[x]]
                    except KeyError as e:
                        #print e
                        continue
                mate1 = seq + reverse_comp_r2

                read_obj["mate1"] = mate1


                #construct Mate2
                #Reverse complement read1 and concatenate it with read2
                reverse_1 = seq[::-1]
                reverse_comp_r1= ""
                for x in range(0, len(reverse_1)):
                    try:
                        reverse_comp_r1 += compMapping[reverse_1[x]]
                    except KeyError as e:
                        #print e
                        continue

                mate2 = reverse_comp_r1 + seq_2

                read_obj["mate2"] = mate2

                #reversed = seq[::-1]
                #print reversed
                # rev_comp = ""
                # for x in range(0, len(reversed)):
                #     try:
                #         #print reversed[x], compMapping[reversed[x]]
                #         rev_comp += compMapping[reversed[x]]
                #     except KeyError as e:
                #         #if the letter is not in the MAP disreguard it as a read all together.
                #         print e
                #     #rev_comp += compMapping[reversed[x]]
                # read_obj["rev_comp_seq"] = rev_comp
                sequences.append(read_obj)
    return sequences#, qualities

def prettyPrintObject(sequences):
    for x in sequences:
        print "#" * 100
        print "ID: ", x["id"]
        print "ID1: ", x["id1"]
        print "ID2: ", x["id2"]
        print "Mate1: ", x["mate1"]
        print "Mate2: ", x["mate2"]
        print "R1 Length: ", x["r1_len"]
        # print "S : ", x["seq"]
        # print "RC: ", x["rev_comp_seq"]
        print "#" * 100
        print ""
#prettyPrintObject(readFastq("dummyTest.fastq"))

prettyPrintObject(readFastq("r1_short.fq", "r2_short.fq"))