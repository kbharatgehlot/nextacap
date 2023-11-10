#! /usr/bin/env python3
import numpy as np
import sys
from argparse import ArgumentParser
parser = ArgumentParser("Convert sage solutions to numpy, create file with eff _nr");


parser.add_argument('-c','--cluster_file', help='sagecal cluster file',dest='cluster_file')
parser.add_argument('--eff_nr', help='output effective solutions file',default = "eff_nr.npy")

def main(argv):
    args=parser.parse_args(argv)
    nrs = []
    for data in open(args.cluster_file):
        if data.strip() and not data.strip()[0]=="#":
            nr = int(data.strip().split()[1])
            nrs.append(nr)
    np.save(args.eff_nr,nrs)

if __name__ == '__main__':
    main(sys.argv[1:])    
