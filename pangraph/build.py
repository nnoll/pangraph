"""
build a pangenome alignment from an annotated guide tree
"""
import os, sys
import builtins

from .utils import mkdir, log
from .tree import Tree

def open(path, *args, **kwargs):
    if path == '-':
        return sys.stdin
    return builtins.open(path, *args, **kwargs)

def register_args(parser):
    parser.add_argument("-d", "--dir",
                        metavar="directory",
                        type=str,
                        default=".",
                        help="directory used for output files")
    parser.add_argument("-l", "--len",
                        metavar="cutoff length",
                        type=int,
                        default=50,
                        help="minimum block size for nucleotides")
    parser.add_argument("-m", "--mu",
                        metavar="block cut chemical potential",
                        type=int,
                        default=100,
                        help="energy cost for cutting blocks (used during block merges)")
    parser.add_argument("-x", "--extensive",
                        default=False,
                        action='store_true',
                        help="boolean flag that toggles whether the energy associated to a cut grows proportionally to the number of sequences cut")
    parser.add_argument("-b", "--beta",
                        metavar="block diversity cost",
                        type=int,
                        default=22,
                        help="energy cost for mutations (used during block merges)")
    parser.add_argument("-w", "--window",
                        metavar="edge window",
                        type=int,
                        default=1000,
                        help="amount of sequence to align from for end repair")
    parser.add_argument("-e", "--extend",
                        metavar="edge extend",
                        type=int,
                        default=1000,
                        help="amount of sequence to extend for end repair")
    parser.add_argument("-c", "--circular",
                        default=False,
                        action='store_true',
                        help="toggle to consider a set of circular genomes")
    parser.add_argument("-s", "--statistics",
                        default=False,
                        action='store_true',
                        help="boolean flag that toggles whether the graph statistics are computed for intermediate graphs")
    parser.add_argument("-n", "--num",
                        type=int,
                        default=-1,
                        help="manually sets the tmp directory number. internal use only.")
    parser.add_argument("input",
                        type=str,
                        default="-",
                        help="input guide tree [json]")

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
    log("reading json")
    with open(args.input, 'r') as input:
        T = Tree.from_json(input)

    root = args.dir.rstrip('/')
    tmp  = f"{root}/tmp"
    if args.num == -1:
        i    = 0
        while os.path.isdir(tmp) and i < 64:
            i += 1
            tmp = f"{root}/tmp{i:03d}"
    else:
            tmp = f"{root}/tmp{args.num:03d}"
    mkdir(tmp)

    log("aligning")
    T.align(tmp, args.len, args.circular, args.mu, args.beta, args.extensive, args.window, args.extend, args.statistics)
    # TODO: when debugging phase is done, remove tmp directory

    graphs = T.collect()

    for i, g in enumerate(graphs):
        log(f"graph {i}: nseqs: {len(g.seqs)} nblks: {len(g.blks)}")

    for i, g in enumerate(graphs):
        with open(f"{root}/graph_{i:03d}.fa", 'w') as fd:
            g.write_fasta(fd)

    # NOTE: uncomment when done debugging
    # T.write_json(sys.stdout, no_seqs=True)

    return 0
