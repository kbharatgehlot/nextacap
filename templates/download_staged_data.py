#!/usr/bin/env python3
"""
A script do download data from the LTA using the `html.txt` provided after the data staging process

.. note::

    The downloaded data is distributed across the given nodes, with the last node getting any remainder. Adjust `-nfiles_per_node` accordingly and of course make sure there is enough disk space on all the nodes used.

How to use this module

.. code-block:: python

    python3 /home/users/chege/theleap/leap/templates/download_staged_data.py --nodes_list ${nodeslist} --html_links_list ${htmltxt} --data_path ${datapath} --nfiles_per_node ${nfiles_per_node} --label ${label}

"""


import os
import asyncio
from argparse import ArgumentParser

parser = ArgumentParser(
    prog="download_staged_data",
    description="Download staged data from the LTA using the `html.txt` file provided after staging",
    usage="python3 /home/users/chege/theleap/leap/templates/download_staged_data.py --nodes_list ${nodeslist} --html_links_list ${htmltxt} --data_path ${datapath} --nfiles_per_node ${nfiles_per_node} --label ${label}",
    epilog="====================DISCLAIMER: Make sure you have enough disk space on each node====================",
)
parser.add_argument(
    "-i",
    "--html_links_list",
    required=False,
    type=str,
    help="html links list provided after the staging",
)
parser.add_argument(
    "--nodes_list",
    required=False,
    help="The list of nodes to which the data will be distributed",
)
parser.add_argument(
    "-d",
    "--data_path",
    type=str,
    required=False,
    help="The data path per node where the data will be placed",
)

parser.add_argument(
    "--label",
    type=str,
    required=False,
    default=None,
    help="An optional label to be appended to MS names just before the `.MS` extension",
)

parser.add_argument(
    "-n",
    "--nfiles_per_node",
    type=int,
    required=False,
    help="The number of files to be downloaded to each node. This is assumed to be  constant for the first n-1 nodes  in the nodeslist and the remainder goes to the last node",
)


def _chunks(lst, n):
    """Yield successive n-sized chunks of msfiles from a list of all msfiles."""
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def background(f):
    def wrapped(*args, **kwargs):
        return asyncio.get_event_loop().run_in_executor(None, f, *args, **kwargs)

    return wrapped


def _read_txt_to_list(txt):
    file_list = open(txt).readlines()
    file_list = [i.strip() for i in file_list]
    return file_list


def _getRealMSNameFromHtmlLink(MS_link):
    MSname = MS_link.split("/")[-1]
    MSname = MSname.split("MS")[0] + "MS"
    return MSname


def downloadMultipleMS(file_list):
    for i, MS_tar_link in enumerate(file_list):
        downloadSingleMS(MS_tar_link)
    return


# @background
def downloadSingleMS(MS_tar_link):
    """Download a single LOFAR Measurement Set (a subband) given its https link from the staging html file output

    Parameters
    ----------
    args : string
        The only required argument is --MS_tar_link : the ms tarball link
    """

    def _doDownload():
        print(f"Downloading {MS_tar_link} into {os.path.abspath('.')}/{MSname}")
        download_command = f"wget --no-check-certificate -c -O {MSname} {MS_tar_link}"
        untar_command = f"tar -xf {MSname}"

        os.system(download_command)
        os.system(untar_command)

        if label:
            finalMSname = MSname.replace(".MS", f"_{label}.MS")
            rename_MS_command = f"mv {MSname} {finalMSname}"
            os.system(rename_MS_command)

    MSname = _getRealMSNameFromHtmlLink(MS_tar_link)
    if label:
        finalMSname = MSname.replace(".MS", f"_{label}.MS")
        if os.path.isdir(finalMSname):
            print(finalMSname, "exists!")
        else:
            _doDownload()
    else:
        if os.path.isdir(MSname):
            print(MSname, "exists!")
        else:
            _doDownload()

    return


def main(args):
    file_list = _read_txt_to_list(args.html_links_list)
    msfiles_per_node = list(_chunks(file_list, args.nfiles_per_node))

    nodes_list = [str(i) for i in args.nodes_list.split(",")]

    global label
    label = args.label

    for _, (node, msfiles) in enumerate(zip(nodes_list, msfiles_per_node)):
        node_path = f"/net/node{node}/{args.data_path}"
        if not os.path.isdir(node_path):
            os.makedirs(node_path, exist_ok=True)
        os.chdir(node_path)

        downloadMultipleMS(msfiles)
        print(f"Node{node} done")


if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
