#! /usr/bin/env python3
import glob
import numpy as np
import os
import sys
from argparse import ArgumentParser
import PosTools

parser = ArgumentParser("Convert sage solutions to numpy")

parser.add_argument("-o", "--obsid", help="obs Id", dest="obsid", required=True)
parser.add_argument(
    "-m", "--MS", nargs="+", help="MS for metadata", dest="MS", required=True
)
parser.add_argument(
    "-c",
    "--cluster",
    nargs=2,
    help="cluster indices (start,end)",
    dest="cluster",
    default=[0, 1],
    type=int,
)
parser.add_argument(
    "-p",
    "--indir",
    help="directory for the input",
    dest="indir",
    default="/data/users/lofareor/sarod/pipeline/",
)
parser.add_argument(
    "--pid", type=int, help="process step id (default =2)", dest="pid", default=2
)
parser.add_argument(
    "-d",
    "--outdir",
    help="directory for the output",
    dest="outdir",
    default="/net/node131/data/users/lofareor/NCP/numpy_data/",
)
parser.add_argument(
    "-n",
    "--nodelist",
    nargs="+",
    help="nodes. On EOR-cluster: specify consecutive multiple nodes as e.g. 116..131 (will use 'node116' to \node'131)",
    dest="nodelist",
    default=[os.getenv("HOSTNAME")],
)
parser.add_argument("--eff_nr", help="eff solution file", default="eff_nr.npy")


def main(argv):
    args = parser.parse_args(argv)
    alldatas = []
    freqs = []
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)
    nodes = [
        "node%03d" % i
        for j in args.nodelist
        if ".." in j
        for i in range(int(j.split("..")[0]), int(j.split("..")[1]) + 1)
    ]
    nodes += [j for j in args.nodelist if not ".." in j]
    timerange, timestep, pointing, stations, station_pos = PosTools.getMSinfo(
        args.MS[0]
    )
    for inode in nodes:
        fnames = "/net/%s/%s/%s*%03d*.MS.solutions" % (
            inode,
            args.indir,
            args.obsid,
            args.pid,
        )
        if args.pid == 0:
            fnames = "/net/%s/%s/%s*MS.solutions" % (inode, args.indir, args.obsid)
        if args.pid < 0:
            fnames = "/net/%s/%s/%s*%02d*.MS.solutions" % (
                inode,
                args.indir,
                args.obsid,
                -1 * args.pid,
            )
        print("getting data", fnames)
        myl = glob.glob(fnames)
        myl = [i for i in myl if not "filtered" in i]
        for fname in myl:
            print(inode, fname)
            try:
                myf = open(fname)
                for i in range(3):
                    a = myf.readline()
                freq, bw, timestep, nStations, nClust, nClustEff = tuple(
                    [float(i) for i in a.split()]
                )
                myf.close()
                data = np.loadtxt(
                    fname,
                    skiprows=3,
                    usecols=tuple(range(1, int(nClustEff) + 1)),
                    unpack=True,
                )
                datas = []
                for i in range(int(nClustEff)):
                    a = data[i].reshape((-1, int(nStations), 4, 2))
                    cdata = a[:, :, :, 0] + 1.0j * a[:, :, :, 1]
                    datas.append(cdata)
            except:
                print("filename", fname, "not correct!!")
                continue
            freqs.append(freq)
            alldatas.append(datas)
    eff_nr = np.load(args.eff_nr)
    has_eff_nr = False
    if np.sum(eff_nr) == nClustEff:
        has_eff_nr = True
    mysorted0, mysorted1 = zip(*sorted(zip(freqs, alldatas)))
    if not os.path.isdir("%s/%s" % (args.outdir, args.obsid)):
        os.makedirs("%s/%s" % (args.outdir, args.obsid))

    cdata = (np.array(mysorted1)[:, args.cluster[0] : args.cluster[1]]).transpose(
        (2, 3, 0, 4, 1)
    )
    data = np.zeros(cdata.shape + (2,), dtype=np.float64)
    data[:, :, :, :, :, 0] = np.real(cdata)
    data[:, :, :, :, :, 1] = np.imag(cdata)
    meantimestep = (timerange[1] - timerange[0]) / (data.shape[0] - 1)
    if has_eff_nr:
        np.savez(
            "%s/%s" % (args.outdir, args.obsid),
            freqs=np.array(mysorted0) * 1.0e6,
            timerange=timerange,
            timestep=timestep,
            meantimestep=meantimestep,
            stations=stations,
            stat_pos=station_pos,
            pointing=pointing,
            eff_nr=eff_nr,
        )
    else:
        np.savez(
            "%s/%s" % (args.outdir, args.obsid),
            freqs=np.array(mysorted0) * 1.0e6,
            timerange=timerange,
            timestep=timestep,
            meantimestep=meantimestep,
            stations=stations,
            stat_pos=station_pos,
            pointing=pointing,
        )

    np.save("%s/%s" % (args.outdir, args.obsid), data)


if __name__ == "__main__":
    main(sys.argv[1:])
