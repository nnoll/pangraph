"""
builds a guide tree utilized by the pan-genome alignment
"""
import os
import io
import json
import subprocess as spawn
import tempfile

import numpy as np
import matplotlib.pylab as plt

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .tree  import Tree
from .utils import openany, parse_fasta

def register_args(parser):
    parser.add_argument("-d", "--dir",
                        metavar="directory",
                        type=str,
                        default=".",
                        help="directory for output file")
    parser.add_argument("-b", "--backend",
                        metavar="backend",
                        type=str,
                        nargs='?',
                        default="mash",
                        choices=['mash', 'random'],
                        help="backend used to estimate inter-sequence distance")
    parser.add_argument("input",
                        type=str,
                        nargs='*',
                        default="-",
                        help="fasta file(s) to cluster")

# TODO: make backends just a single function...

# ------------------------------------------------------------------------
# mash backend

# NOTE: mash takes '-' as filename if it is to read from stdin
def run_mash(inpath):
    stdout = spawn.check_output(f"mash triangle {inpath} 2>/dev/null", shell=True)
    return io.StringIO(stdout.decode("utf-8"))

def parse_mash(input):
    nrows = int(input.readline().strip())
    M     = np.zeros((nrows, nrows), dtype=float)
    names = []
    for i, line in enumerate(input):
        e = line.strip().split('\t')
        names.append(e[0])
        M[i,:i] = [float(x) for x in e[1:]]

    M  = M + M.T

    input.close()
    return M, np.array(names)

# ------------------------------------------------------------------------
# random matrix backend

def nop(input):
    return input

def random_matrix(input):
    np.random.seed(0)
    with openany(input) as fh:
        names = [rec.name for rec in parse_fasta(fh)]

    M = np.abs(np.random.randn(len(names),len(names)))
    M = M + M.T

    return M, np.array(names)


# ------------------------------------------------------------------------
# backends perform the pairwise distance approximation

class Backend(object):
    def __init__(self, run, parse):
        self.run   = run
        self.parse = parse

backends = {
    "mash"   : Backend(run=run_mash, parse=parse_mash),
    "random" : Backend(run=nop, parse=random_matrix),
    # add more backends here...
}

# ------------------------------------------------------------------------
# helpers

def asnwkstr(node, nwk, pdist, leafs):
    if node.is_leaf():
        return "%s:%.2f%s" % (leafs[node.id], pdist - node.dist, nwk)

def export(cls, names, fname):
    C = np.max(cls)
    with open(fname, "w+") as out:
        for c in range(C):
            out.write("\t".join(names[cls == (c+1)]) + "\n")

def escape(name):
    return name.replace(".", "-").replace("/", "#").replace("'", "").replace(",","")

# ------------------------------------------------------------------------
# Main point of entry

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

    backend = backends[args.backend]
    outdir  = args.dir.rstrip("/")

    tmp = None
    if isinstance(args.input, list):
        tmp = tempfile.NamedTemporaryFile('w', delete=True)
        for file in args.input:
            with openany(file, 'r') as fd:
                seqs = [SeqRecord(id = escape("_".join(s.description.split())), name = s.description, seq = s.seq) for s in SeqIO.parse(fd, 'fasta')]
                SeqIO.write(seqs, tmp, "fasta")
        input = tmp.name
    else:
        input = args.input

    dist, names = backend.parse(backend.run(input))
    tree = Tree.nj(dist, names)

    with openany(input, 'r') as fd:
        seqs = {s.id : s.seq for s in SeqIO.parse(fd, "fasta")}
    tree.attach(seqs)

    if tmp:
        tmp.close()

    # exports:
    np.savez(f"{outdir}/distmtx.npz", dist=dist, names=names)

    # for now we output both forms of the tree
    with open(f"{outdir}/guide.nwk", "w") as fd:
        tree.write_nwk(fd)

    with open(f"{outdir}/guide.json", "w") as fd:
        tree.write_json(fd)

    return 0
