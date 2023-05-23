#! /usr/bin/env python3
from argparse import ArgumentParser
import casacore.tables as pt
import numpy as np
import sys
import os

parser = ArgumentParser("Create DATA like columns to multiple MS")
parser.add_argument(
    "-i",
    "--MSlist",
    nargs="+",
    help="Measurement sets or file with list of MS",
    dest="MSlist",
    default="",
)
parser.add_argument(
    "-c",
    "--colnames",
    nargs="+",
    help="lsit of columns to create",
    dest="colnames",
    default="",
)


def main(argv):
    args = parser.parse_args(argv)
    colnames = args.colnames
    MSlist = [i for i in args.MSlist if pt.tableexists(i)]
    if len(MSlist) == 0:
        MSlists = []
        for ifile in args.MSlist:
            if os.path.isfile(ifile):
                myf = open(ifile)
                MSlists.append([i.strip() for i in myf if pt.tableexists(i.strip())])
    else:
        MSlists = [MSlist]
    MSlist = MSlists[0]
    for ms in MSlist:
        myt = pt.table(ms, readonly=False)
        data = myt.getcol("DATA")
        desc = myt.getcoldesc("DATA")
        for col in colnames:
            if col in myt.colnames():
                # myt.removecols(col)
                print(f"MS already has {col} column. Not removing/readdding")
            else:
                desc["name"] = col
                myt.addcols(desc)
                myt.putcol(col, np.zeros_like(data))
                print(f"{col} column added to {ms}")

        myt.close()


if __name__ == "__main__":
    main(sys.argv[1:])
