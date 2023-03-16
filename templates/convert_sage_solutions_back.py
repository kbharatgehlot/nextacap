#! /usr/bin/env python3
import glob
import numpy as np
import os
import sys
from argparse import ArgumentParser

parser = ArgumentParser("Convert sage solutions to numpy")

parser.add_argument("-o", "--obsid", help="obs Id", dest="obsid", required=True)
parser.add_argument(
    "-p",
    "--indir",
    help="directory for the input solfiles",
    dest="indir",
    default="/data/users/lofareor/sarod/pipeline/",
)
parser.add_argument(
    "-d",
    "--data",
    help="input numpy files",
    dest="data",
    default="/data/users/lofareor/NCP/L90490/L90490_sol_250k/L90490.npy",
)
parser.add_argument(
    "-l",
    "--label",
    nargs="?",
    help="extra label of the data (e.g. _NP)",
    dest="label",
    default="",
    const="",
)
parser.add_argument(
    "-m",
    "--metadata",
    help="input numpy z files",
    dest="metadata",
    default="/data/users/lofareor/NCP/L90490/L90490_sol_250k/L90490.npz",
)
parser.add_argument(
    "--pid", type=int, help="process step id (default =2)", dest="pid", default=2
)
parser.add_argument(
    "-n",
    "--nodelist",
    nargs="+",
    help="nodes. On EOR-cluster: specify consecutive multiple nodes as e.g. 116..131 (will use 'node116' to \node'131)",
    dest="nodelist",
    default=[os.getenv("HOSTNAME")],
)


def main(argv):
    args = parser.parse_args(argv)
    # alldatas = []
    freqs = []
    new_data = np.load(args.data).transpose((2, 4, 0, 1, 3, 5))
    metadata = np.load(args.metadata)
    freqs = metadata["freqs"] * 1.0e-6
    nodes = [
        "node%03d" % i
        for j in args.nodelist
        if ".." in j
        for i in range(int(j.split("..")[0]), int(j.split("..")[1]) + 1)
    ]
    for inode in nodes:
        print(
            "getting data",
            "/net/%s/%s/%s*%03d%s.MS.solutions"
            % (inode, args.indir, args.obsid, args.pid, args.label),
        )
        myl = glob.glob(
            "/net/%s/%s/%s*%03d%s.MS.solutions"
            % (inode, args.indir, args.obsid, args.pid, args.label)
        )
        for fname in myl:
            print(inode, fname)
            # if True:
            try:
                myf = open(fname)
                newfname = fname.replace(".MS.", "_filtered.MS.")
                header = ""
                for i in range(3):
                    a = myf.readline()
                    header += a
                freq, bw, timestep, nStations, nClust, nClustEff = tuple(
                    [float(i) for i in a.split()]
                )
                print(header.strip())
                myf.close()
                data = np.loadtxt(
                    fname,
                    skiprows=3,
                    usecols=tuple(range(0, int(nClustEff) + 1)),
                    unpack=True,
                )
                freqidx = np.argmin(np.abs(freqs - freq))
                print("freqidx", freqidx, freqs[freqidx], freq)
                print(
                    "filling data",
                    data.shape,
                    new_data.shape,
                    new_data[freqidx].reshape((int(nClustEff), -1)).shape,
                )
                data[1:] = new_data[freqidx].reshape((int(nClustEff), -1))
                np.savetxt(
                    newfname,
                    data.T,
                    fmt="%d " + " %.6e" * int(nClustEff),
                    header=header.strip(),
                    comments="",
                )
            except:
                print("filename", fname, "not correct!!")
                continue


if __name__ == "__main__":
    main(sys.argv[1:])
