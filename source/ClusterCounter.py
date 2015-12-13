__author__ = 'soumadipmukherjee'
#!/usr/bin/env python
from collections import namedtuple
Defaults = namedtuple('Defaults', ['count'])
defaults = Defaults(count=1)

def main():
    args = parse_args()
    print args
    one_c, t_c = countClusters(args.file1)
    print one_c
    print t_c
    return
def countClusters(clusters):
    print clusters
    count = 0
    runcount = 0
    cluster_with_one = 0;
    while True:
        r = clusters.readline()
        if r == '\n' or len(r) == 1:
            count+=1

            if runcount == 1:
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