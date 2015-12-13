__author__ = 'soumadipmukherjee'
#!/usr/bin/env python
from collections import namedtuple
Defaults = namedtuple('Defaults', ['offset'])
defaults = Defaults(offset=1000)
def main():
    print "HELLO WORLD"
    args = parse_args()
    generateNewFastq(args.offset, args.file1, args.fout)
    #print args
    return

def generateNewFastq(offset, file1, fout):
    #print offset
    #print file1
    #print fout
    count = 0
    #fout = open()
    while True:
        id = file1.readline()
        seq = file1.readline()
        plus = file1.readline()
        quality = file1.readline()
        if len(id) == 0 or len(seq) == 0 or len(plus) == 0 or len(quality) == 0:
            break

        if count % offset == 0:
            print count
            fout.write(id)
            fout.write(seq)
            fout.write(plus)
            fout.write(quality)
        count += 1

    fout.close()
    file1.close()

    return
def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Generate fastq file with offset reads from input files')
    parser.add_argument('-o', '--file-offset-size', action='store',
                       dest='offset', type=int, default=defaults.offset,
                       help="offset, default: {}".format(defaults.offset))
    parser.add_argument('file1', action='store', type=argparse.FileType('r'),
                       help="FASTQ input file")
    parser.add_argument('fout', action='store', type=argparse.FileType('w'),
                        help="FASTQ output file")
    return parser.parse_args()

if __name__=='__main__':
    main()