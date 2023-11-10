#! /usr/bin/env python3
from argparse import ArgumentParser
import casacore.tables as pt
import numpy as np
import sys
import os

parser = ArgumentParser(
    "add the the direction dependent model data column to measurement sets by subtracting the calibrated data residuals column from the uncalibrated data column"
)
parser.add_argument(
    "-i",
    "--MSlist",
    nargs="+",
    help="Measurement sets or file with list of MS",
    dest="MSlist",
    default="",
)
parser.add_argument(
    "-d",
    "--data_colname",
    help="list of columns to create",
    dest="data_colname",
    default="",
)
parser.add_argument(
    "-c",
    "--corrected_colname",
    help="list of columns to create",
    dest="corrected_colname",
    default="",
)
parser.add_argument(
    "-m",
    "--model_colname",
    help="list of columns to create",
    dest="model_colname",
    default="",
)


def main(argv):
    args = parser.parse_args(argv)

    data_colname = args.data_colname
    corrected_colname = args.corrected_colname
    model_colname = args.model_colname

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

        if model_colname in myt.colnames():
            myt.removecols(model_colname)
            print(f"MS already has {model_colname} column. Removing it first")

        data = myt.getcol(data_colname)
        corrected_data = myt.getcol(corrected_colname)
        model = data - corrected_data

        desc = myt.getcoldesc(data_colname)
        desc["name"] = model_colname
        myt.addcols(desc)
        myt.putcol(model_colname, model)

        print(f"{model_colname} column added to {ms}")
        myt.close()


if __name__ == "__main__":
    main(sys.argv[1:])
