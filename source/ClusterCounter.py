#!/usr/bin/env python
from collections import namedtuple

Defaults = namedtuple('Defaults', ['count'])
defaults = Defaults(count=1)


def main():
    args = parse_args()
    #print args
    one_c, t_c = countClusters(args.file1, args.count)
    print "Clusters with less then "+str(args.count+ 1)+" reads: " + str(one_c)
    print "Total clusters: " + str(t_c)


def countClusters(clusters, count_min ):
    print clusters
    count = 0
    runcount = 0
    cluster_with_one = 0
    while True:
        r = clusters.readline()
        if r == '\n':
            count+=1
            if runcount <= count_min:
                cluster_with_one += 1
            runcount=0
        else:
            runcount +=1
        if len(r) == 0:
            return cluster_with_one, count


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Clusters counter with small reads')
    parser.add_argument('-c', '--cluster count', action='store', dest='count',
                        type=int, default=defaults.count,
                        help='cluster count, default:{}'.format(defaults.count))
    parser.add_argument('file1', action='store',
                       type=argparse.FileType('r'),
                       help="Cluster file")
    return parser.parse_args()


if __name__ =="__main__":
    main()
