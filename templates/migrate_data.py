#!/usr/bin/env python3

import glob
import os
import sys
from argparse import ArgumentParser

parser = ArgumentParser(description="Migrate data")
parser.add_argument(
    "-a",
    "--from_nodes",
    help="List of nodes where the data is now",
    nargs="+",
    required=True,
)
parser.add_argument(
    "-b",
    "--to_nodes",
    help="List of nodes where to migrate the data to",
    nargs="+",
    required=True,
)
parser.add_argument(
    "-d", "--data", help="glob parttern of the data to be moved", required=False
)
parser.add_argument("-x", "--directory_A", help="from this directory", required=False)
parser.add_argument("-y", "--directory_B", help="to this directory", required=False)
parser.add_argument(
    "-r",
    "--remove_dir_first",
    help="Will remove `--directory_B` first before copying data. use cautiously!!!!",
    action="store_true",
)
parser.add_argument(
    "-n",
    "--nfiles_per_worker",
    type=int,
    help="number of files o distribute to each worker node. the master node gets the remainder",
)

parser.add_argument(
    "-t",
    "--dry_run",
    help="test first before actually moving the data",
    action="store_true",
)


def chunks(lst, n):
    """Yield successive n-sized chunks of msfiles from a list of all msfiles."""
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def remove_dir(nodeB, datadirB):
    assert os.path.isdir(f"/net/node{nodeB}{datadirB}")
    rm_cmd = f"ssh node{nodeB} 'rm -r {datadirB}'"
    print(f"running {rm_cmd}")
    os.system(rm_cmd)


def main(args):
    if len(args.from_nodes) == len(args.to_nodes):
        print("scenario 1")
        for n, (nodeA, nodeB) in enumerate(zip(args.from_nodes, args.to_nodes)):
            assert len(args.from_nodes) == len(args.to_nodes), "Mismatch between nodes"
            datadirA = f"{args.directory_A}"
            datadirB = f"{args.directory_B}"

            assert os.path.isdir(
                f"/net/node{nodeA}{datadirA}"
            ), f"/net/node{nodeA}{datadirA} not a valid path"
            files = glob.glob(f"/net/node{nodeA}{datadirA}/{args.data}")

            assert len(
                files
            ), f"no files found at /net/node{nodeA}{datadirA}/{args.data}"

            print(len(files), f"files found in node{nodeA}")

            if args.remove_dir_first:
                print("removing DirectoryB firts")
                if not args.dry_run:
                    remove_dir(nodeB, datadirB)

            mkdir_cmd = f"ssh node{nodeB} 'mkdir -p {datadirB}'"
            mv_cmd = f"ssh node{nodeB} 'scp -r node{nodeA}:{datadirA}/{args.data} {datadirB}'"

            print(f"running {mkdir_cmd}")
            if not args.dry_run:
                os.system(mkdir_cmd)

            print(f"running {mv_cmd}")
            if not args.dry_run:
                os.system(mv_cmd)

    elif len(args.from_nodes) == 1:
        print("Scenario 2: Moving data from single node to multiple nodes")
        assert args.nfiles_per_worker, "Provide a number of files for each worker"

        datadirA = args.directory_A
        nodeA = args.from_nodes[0]

        files = sorted(glob.glob(f"/net/node{nodeA}{datadirA}/{args.data}"))

        assert len(files), f"no files found at /net/node{nodeA}{datadirA}/{args.data}"

        print(f"Found {len(files)} files in node{nodeA}")

        files_per_node = list(chunks(files, args.nfiles_per_worker))

        assert len(files_per_node) == len(
            args.to_nodes
        ), f"chuncks ({len(files_per_node)}) less than nodes {len(args.to_nodes)}, reduce files per node accordingly"

        print(args.to_nodes)
        for n, nodeB in enumerate(args.to_nodes):
            nodeB = int(nodeB.strip("[").strip("]").strip(","))
            datadirB = f"{args.directory_B}"
            files_chunk = files_per_node[n]

            if args.remove_dir_first:
                print("removing DirectoryB firts")
                if not args.dry_run:
                    remove_dir(nodeB, datadirB)

            mkdir_cmd = f"ssh node{nodeB} 'mkdir -p {datadirB}'"
            print(f"running {mkdir_cmd}")
            if not args.dry_run:
                os.system(mkdir_cmd)

            # assert os.path.isdir(datadirB)

            for fyl in files_chunk:
                # fyl2 = fyl.split("/")[-1].replace("001", "002")
                mv_cmd = f"ssh node{nodeB} 'scp -r node{nodeA}:{fyl} {datadirB}'"
                print(f"running {mv_cmd}")
                if not args.dry_run:
                    os.system(mv_cmd)

    else:
        print("scenario 3")
        assert args.nfiles_per_worker, "Provide a number of files for each worker"
        datadirA = args.directory_A

        files = []

        for nodeA in args.from_nodes:
            nodeA = int(nodeA.strip("[").strip("]").strip(","))
            files_in_this_node = glob.glob(f"/net/node{nodeA}{datadirA}/{args.data}")
            assert (
                len(files_in_this_node) > 0
            ), f"No files found at {f'/net/node{nodeA}{datadirA}/{args.data}'}"
            print(f"Node {nodeA}: {len(files_in_this_node)} files")
            files += files_in_this_node

        files = sorted(files)
        print(len(files), f"All files found in nodes {args.from_nodes}")

        files_per_node = list(chunks(files, args.nfiles_per_worker))

        print(args.to_nodes)
        for n, nodeB in enumerate(args.to_nodes):
            nodeB = int(nodeB.strip("[").strip("]").strip(","))
            datadirB = f"{args.directory_B}"
            files_chunk = files_per_node[n]

            if args.remove_dir_first:
                remove_dir(nodeB, datadirB)

            mkdir_cmd = f"ssh node{nodeB} 'mkdir -p {datadirB}'"
            print(f"running {mkdir_cmd}")
            if not args.dry_run:
                os.system(mkdir_cmd)

            # assert os.path.isdir(datadirB)

            for fyl in files_chunk:
                mv_cmd = f"ssh node{nodeB} 'cp -r {fyl} {datadirB}'"
                print(f"running {mv_cmd}")
                if not args.dry_run:
                    os.system(mv_cmd)

    print("All done !")


if __name__ == "__main__":
    args = parser.parse_args()
    # print(args)
    main(args)
