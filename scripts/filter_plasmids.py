#!/usr/bin/env python3
"""
script to filter plasmids and chromosomes from full genome assemblies
"""

import os
import sys
import gzip
import builtins
import argparse

from glob import glob

sys.path.insert(0, os.path.abspath('.')) # gross hack
from pangraph.utils import parse_fasta, breakpoint

def open(path, *args, **kwargs):
    if path.endswith('.gz'):
        return gzip.open(path, *args, **kwargs)
    else:
        return builtins.open(path, *args, **kwargs)

def main(dirs, plasmids=True):
    for d in dirs:
        in_dir = f"data/{d}/assemblies"
        if not os.path.exists(in_dir):
            print(f"{in_dir} doesn't exist. skipping...")
            continue

        if plasmids:
            out_dir = f"data/{d}-plasmid/assemblies"
        else:
            out_dir = f"data/{d}-chromosome/assemblies"

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        is_plasmid = lambda r: len(r.seq) < 1e6 and "plasmid" in r.meta

        for path in glob(f"{in_dir}/*.f?a*"):
            with open(path, 'rt') as fd, open(f"{out_dir}/{os.path.basename(path).replace('.gz', '')}", 'w') as wtr:
                for i, rec in enumerate(parse_fasta(fd)):
                    if is_plasmid(rec):
                        if plasmids:
                            wtr.write(str(rec))
                            wtr.write('\n')
                    else:
                        if not plasmids:
                            wtr.write(str(rec))
                            wtr.write('\n')

parser = argparse.ArgumentParser(description='seperate plasmids from chromosomes')
parser.add_argument('directories', metavar='dirs', nargs='+')
parser.add_argument('--chromosomes', default=False, action='store_true')

if __name__ == "__main__":
    args = parser.parse_args()
    main(args.directories, plasmids=not args.chromosomes)
