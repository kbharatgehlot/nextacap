#! /usr/bin/env python3
from argparse import ArgumentParser
import pyrap.tables as tab
import sys
import os
import numpy as np
import subprocess

parser = ArgumentParser("get clip value from XY,YX rms of data")
parser.add_argument(
    "-i",
    "--MSlist",
    nargs="+",
    help="Measurement sets or file with list of MS",
    dest="MSlist",
    default="",
)
parser.add_argument(
    "-t",
    "--table_name",
    help="table name,default=table.f4_TSM0",
    dest="table_name",
    default="table.f4_TSM0",
)
parser.add_argument(
    "-o",
    "--outname",
    help="output name,default=table.f4_TSM0_copy",
    dest="outname",
    default="table.f4_TSM0_copy",
)


def main(argv):
    args = parser.parse_args(argv)
    MSlist = [i for i in args.MSlist if tab.tableexists(i)]
    if len(MSlist) == 0:
        MSlists = []
        for ifile in args.MSlist:
            if os.path.isfile(ifile):
                myf = open(ifile)
                MSlists.append([i.strip() for i in myf if tab.tableexists(i.strip())])
    else:
        MSlists = [MSlist]
    for MSlist in MSlists:
        for i in MSlist:
            print("copying")
            myinput = i + "/" + args.table_name
            myoutput = i + "/" + args.outname
            print(myinput, "to", myoutput)
            subprocess.call(["cp", myinput, myoutput])


print("done")

if __name__ == "__main__":
    main(sys.argv[1:])
