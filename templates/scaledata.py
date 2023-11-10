#! /usr/bin/env python3
import casacore.tables as pt
import sys
import numpy as np
from argparse import ArgumentParser
import os

parser = ArgumentParser("scale data ")
parser.add_argument(
    "-i",
    "--MSlist",
    nargs="+",
    help="Measurement sets or file with list of MS",
    dest="MSlist",
    default="",
)
parser.add_argument(
    "-f",
    "--scalevalue",
    help="scalevalue, if <-1 estimate from data. default=-1",
    dest="scalevalue",
    type=float,
    default=-1,
)


def main(argv):
    args = parser.parse_args(argv)
    MSlist = [i for i in args.MSlist if pt.tableexists(i)]
    if len(MSlist) == 0:
        MSlists = []
        for ifile in args.MSlist:
            if os.path.isfile(ifile):
                myf = open(ifile)
                MSlists.append([i.strip() for i in myf if pt.tableexists(i.strip())])
    else:
        MSlists = [MSlist]
    for MSlist in MSlists:
        for i in MSlist:
            myt = pt.table(i, readonly=False)
            data = myt.getcol("DATA")
            if args.scalevalue <= 0:
                flags = myt.getcol("FLAG")
                mdata = np.ma.array(data, mask=flags)
                if (
                    np.sum(flags[:10000, :, (0, 3)])
                    / (np.prod(np.array(flags[:10000, :, (0, 3)].shape)))
                    > 0.3
                ):
                    scalef = np.ma.average(np.ma.absolute(mdata[:10000, :, (0, 3)]))
                else:  # too many datapoints flagged, use all data
                    scalef = np.ma.average(np.ma.absolute(mdata))
                print("this data has scale", i, scalef)
                # this block might need some optimising butit is not an issue..
                if scalef > 100:
                    scalef = 1e-5
                elif scalef < 1e-3:
                    scalef = 1e4
                elif scalef < 5e-2:
                    scalef = 10
                else:
                    scalef = 1
                    myt.close()
                    continue
            else:
                scalef = args.scalevalue

            print("scaling", i, "with", scalef)
            data *= scalef
            myt.putcol("DATA", data)
            myt.close()


if __name__ == "__main__":
    main(sys.argv[1:])

# TODO:
# if the gains are not between 0.5 adnd 10, this script scales the gains to there so that sagecal can work. This is due to the digital and station gains which are unknown.
# the correlator changed and that is why the digital gains chained... we don't know the dates when the correlator changed.

# include the test_scale script. To see if the data was already scaled or not.

# In the end the ultimate check is how stable the gains look
