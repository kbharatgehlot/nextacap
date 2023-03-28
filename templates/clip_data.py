#! /usr/bin/env python3
import pyrap.tables as tab
from argparse import ArgumentParser
import sys
import os
import numpy as np

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
    "-f",
    "--clipvalue",
    help="clipvalue, default=70",
    dest="clipvalue",
    type=float,
    default=70,
)
parser.add_argument(
    "-c", "--column", help="column default=DATA", dest="column", default="DATA"
)
parser.add_argument(
    "-o",
    "--outputcolumn",
    help="out put column default=MODEL_DATA",
    dest="outputcolumn",
    default="MODEL_DATA",
)
parser.add_argument(
    "-r", "--reset_flags", help="set all flags to false", action="store_true"
)
parser.add_argument(
    "-t", "--flag_intrastations", help="flag all intrastations", action="store_true"
)
parser.add_argument(
    "-l",
    "--flag_longbaselines",
    help="flag stations outside 10km of core",
    action="store_true",
)
parser.add_argument(
    "-b", "--flag_badbaselines", help="flag bad CS CS stations", action="store_true"
)
parser.add_argument(
    "-m", "--flag_remote", help="flag all remote stations", action="store_true"
)
parser.add_argument(
    "-n",
    "--do_not_clip",
    help="do not clip nor create data columns",
    action="store_true",
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
            if args.reset_flags:
                print("resetting flags")
                tab.taql("UPDATE $i SET FLAG=False")
            if args.flag_intrastations:
                print("flagging intrastrations")
                tab.taql(
                    r'UPDATE $i SET FLAG=True WHERE mscal.baseline("/(.*)HBA0&\1HBA1/")'
                )  # mind the r
            if args.flag_longbaselines:
                print("flagging long baselines")
                freq = tab.table(i + "/SPECTRAL_WINDOW").getcol("REF_FREQUENCY")[0]
                # tab.taql(r'UPDATE $i SET FLAG=True WHERE mscal.baseline("RS208HBA,RS210HBA,RS310HBA,RS407HBA,RS409HBA,RS508HBA,RS509HBA")')#mind the r
                if freq < 134e6:
                    tab.taql(
                        r'UPDATE $i SET FLAG=True WHERE mscal.baseline("RS406HBA")'
                    )  # mind the r
                elif freq < 146e6:
                    tab.taql(
                        r'UPDATE $i SET FLAG=True WHERE mscal.baseline("RS306HBA,RS406HBA")'
                    )  # mind the r
                else:
                    tab.taql(
                        r'UPDATE $i SET FLAG=True WHERE mscal.baseline("RS406HBA,RS306HBA,RS106HBA")'
                    )  # mind the r

                # tab.taql(r'UPDATE $i SET FLAG=True WHERE mscal.baseline(">10000")')#mind the r
            if args.flag_badbaselines:
                print(
                    "flagging bad CS-CS baselines, only flagging CS013 when needed now"
                )
                obsid = int(i.split("/")[-1].split("_")[0][1:])
                if obsid < 100000:
                    print("obsid %d:  flagging bad CS013" % obsid)
                    tab.taql(
                        r'UPDATE $i SET FLAG=True WHERE mscal.baseline("CS013HBA0,CS013HBA1")'
                    )
                # tab.taql(r'UPDATE $i SET FLAG=True WHERE mscal.baseline("/CS002HBA0&CS005HBA1|CS002HBA1&CS005HBA1|CS002HBA0&CS006HBA1|CS002HBA1&CS006HBA1|CS006HBA0&CS007HBA0|CS003HBA0&CS004HBA1/")')
            if args.flag_remote:
                print("flagging remote stations")
                tab.taql(
                    r'UPDATE $i SET FLAG=True WHERE mscal.baseline("R*&&*")'
                )  # mind the r
            # make sure all non-dutch stations are flaggged
            tab.taql(r'UPDATE $i SET FLAG=True WHERE mscal.baseline("![CR]*&&")')
            if not args.do_not_clip:
                data = tab.table(i).getcol(args.column)
                nanflags = np.isnan(data)
                flags = tab.table(i).getcol("FLAG")
                flags = np.logical_or(nanflags, flags)
                # overwrite old flags
                nflags = np.abs(data[:, :, :]) > args.clipvalue
                # data[nflags]=args.clipvalue
                flags = np.logical_or(nflags, flags)
                print("flagging", np.sum(nflags) / 4, "data points")
                myt = tab.table(i, readonly=False)
                myt.putcol("FLAG", flags)
                if not "MODEL_DATA" in myt.colnames():
                    ddesc = myt.getcoldesc("DATA")
                    ddesc["name"] = "MODEL_DATA"
                    myt.addcols(ddesc)
                    myt.putcol("MODEL_DATA", data)
                if not "CORRECTED_DATA" in myt.colnames():
                    ddesc = myt.getcoldesc("DATA")
                    ddesc["name"] = "CORRECTED_DATA"
                    myt.addcols(ddesc)
                    myt.putcol("CORRECTED_DATA", data)
                # myt.putcol(args.outputcolumn,data)
                myt.flush()


if __name__ == "__main__":
    main(sys.argv[1:])
