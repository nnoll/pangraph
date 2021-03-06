"""
builds a guide tree utilized by the pan-genome alignment
"""
import os, sys
import json
import builtins

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def register_args(parser):
    parser.add_argument("-p", "--prefix",
                        type=str,
                        nargs='?',
                        default=None,
                        help="prefix name of output files. if not given, will output all to stdout")
    parser.add_argument("input",
                        type=str,
                        nargs='?',
                        default="-",
                        help="pangraph json to serialize")

def open(path, *args, **kwargs):
    if path == '-':
        return sys.stdin
    return builtins.open(path, *args, **kwargs)

def main(args):
    '''
    Parameters
    ----------
    args : namespace
        arguments passed in via the command-line from pangraph
    Returns
    -------
    int
        returns 0 for success, 1 for general error
    '''
    with open(args.input, 'r') as fd:
        pg = json.load(fd)['tree']['graph']

    for n, g in enumerate(pg):
        s  = [SeqRecord(id=b["id"], seq=Seq(b["seq"])) for b in graph["blocks"]]
        if args.prefix:
            fa = f"{args.prefix}_{n:03d}.fa"
            with open(fa, 'w') as ofd:
                SeqIO.write(s, ofd, "fasta")
        else:
            print(s)
            if n < len(pg) - 1:
                print('=')
