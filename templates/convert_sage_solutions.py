"""Conversion of sagecal solution to numpy arrays"""

#! /usr/bin/env python3

import os
import sys
import glob
from argparse import ArgumentParser

import numpy as np

import PosTools

parser = ArgumentParser("Convert sagecal solutions to numpy")

parser.add_argument(
    "-o",
    "--obsid",
    help="Observation ID, for example L254871",
    dest="obsid",
    required=True,
)
parser.add_argument(
    "-m",
    "--MS",
    nargs="+",
    help="The full path to the measurement set files that contain the metadata information, e.g. /path/to/file/*.MS",
    dest="MS",
    required=True,
)
parser.add_argument(
    "-c",
    "--cluster",
    nargs=2,
    help="The indices of the directions (clusters) to read from sagecal solution file (start_index, end_index). Default: 0 1, which is the first direction/cluster for direction-independent calibration.",
    dest="cluster",
    default=[0, 1],
    type=int,
)
parser.add_argument(
    "-p",
    "--indir",
    help="The directory of the sagecal solutions files.",
    dest="indir",
    default="/data/users/lofareor/sarod/pipeline/",
)
parser.add_argument(
    "--pid",
    type=int,
    help="Process step ID. Default: 2, which is after direction-independent calibration. For unprocessed data use 0.)",
    dest="pid",
    default=2,
)
parser.add_argument(
    "-d",
    "--outdir",
    help="The directory to save the output numpy file.",
    dest="outdir",
    default="/net/node131/data/users/lofareor/NCP/numpy_data/",
)
parser.add_argument(
    "-n",
    "--nodelist",
    nargs="+",
    help="The list of nodes where sagecal solutions files are saved. On EOR-cluster: specify consecutive multiple nodes as e.g. 116..131 (will use 'node116' to \node'131)",
    dest="nodelist",
    default=[os.getenv("HOSTNAME")],
)
parser.add_argument(
    "--eff_nr",
    help="Effective solution file. This file contains information about the number of solutions per time interval per cluster.",
    default="eff_nr.npy",
)


def main(argv):
    args = parser.parse_args(argv)
    alldatas = []
    freqs = []

    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    # Get the list of nodes
    nodes = [
        "node%03d" % i
        for j in args.nodelist
        if ".." in j
        for i in range(int(j.split("..")[0]), int(j.split("..")[1]) + 1)
    ]

    nodes += [j for j in args.nodelist if not ".." in j]

    # Get metadata information from MS files
    timerange, timestep, pointing, stations, station_pos = PosTools.getMSinfo(
        args.MS[0]
    )

    # Get full path to sagecal solution files, looping over nodes means looping over frequency bands
    for inode in nodes:
        fnames = "/net/%s/%s/%s*%03d*.MS.solutions" % (
            inode,
            args.indir,
            args.obsid,
            args.pid,
        )

        # If process ID is zero, skip process ID in path to file
        if args.pid == 0:
            fnames = "/net/%s/%s/%s*MS.solutions" % (inode, args.indir, args.obsid)
        if args.pid < 0:
            fnames = "/net/%s/%s/%s*%02d*.MS.solutions" % (
                inode,
                args.indir,
                args.obsid,
                -1 * args.pid,
            )

        # Get a list of files from each node
        myl = glob.glob(fnames)
        myl = [i for i in myl if not "filtered" in i]

        for fname in myl:
            print(inode, fname)
            try:
                print("Successfully opened file")
                myf = open(fname)

                # The first lines contain some metadata such as frequencies, bandwidth, time resolution, number of stations, number of clusters, number
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

                # If we only calibrated into 1 dirction / 1 cluster
                if nClustEff == 1.0:
                    a = data.reshape((-1, int(nStations), 4, 2))
                    cdata = a[:, :, :, 0] + 1.0j * a[:, :, :, 1]
                    datas.append(cdata)

                # Else, loop over the number of diretions / clusters
                else:
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

    # Sort data in ascending frequency order
    mysorted0, mysorted1 = zip(*sorted(zip(freqs, alldatas)))

    if not os.path.isdir("%s/%s" % (args.outdir, args.obsid)):
        os.makedirs("%s/%s" % (args.outdir, args.obsid))

    # The calibration data are in shape (nfreq, ncluster, ntime, nstations, npol), after transpose (ntime, nstations, nfreq, npol, ncluster)
    cdata = (np.array(mysorted1)[:, args.cluster[0] : args.cluster[1]]).transpose(
        (2, 3, 0, 4, 1)
    )

    # For some reason the real part and imaginary part are saved to last axis with dimension 2
    data = np.zeros(cdata.shape + (2,), dtype=np.float64)
    data[:, :, :, :, :, 0] = np.real(cdata)
    data[:, :, :, :, :, 1] = np.imag(cdata)

    meantimestep = (timerange[1] - timerange[0]) / (data.shape[0] - 1)

    # Check if the number of solutions per time interval per cluster is different for each cluster
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
