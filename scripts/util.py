import argparse, gzip
import csv
import numpy as np
from collections import defaultdict

def myopen(fname, mode='r'):
    if fname.endswith('.gz'):
        return gzip.open(fname, mode)
    else:
        return open(fname, mode)

def parse_paf(fname):
    assert fname.endswith(".paf")
    hits = []
    with myopen(fname) as fh:
        for line in fh:
            entries = line.strip().split()
            hit = {'query':{'name':entries[0], 'start':int(entries[2]), 'end':int(entries[3])},
                   'ref':  {'name':entries[5], 'start':int(entries[7]), 'end':int(entries[8])},
                   'aligned_bases':int(entries[9]), 'aligned_length':int(entries[10]),
                   'orientation':1 if entries[4]=='+' else -1, 'mapping_quality':int(entries[11])}
            for extra in entries[12:]:
                if extra.startswith('cg:'):
                    hit['cigar'] = extra.split(':')[-1]
            hits.append(hit)
    return hits

def parse_m8(fname):
    assert fname.endswith(".m8")
    hits = defaultdict(list)
    with open(fname) as fh:
        rdr = csv.reader(fh, delimiter = "\t", quotechar = '"')
        for row in rdr:
            hits[row[0]].append((row[11], row[1]))

        hits = dict(hits)
        for key in hits:
            hits[key] = sorted(hits[key], key = lambda x: x[0])[::-1]

    return hits

def parse_kmer(mtx):
    with open(mtx) as fh:
        nrows = int(fh.readline().strip())
        M = np.zeros((nrows, nrows), dtype=float)
        seq_names = []
        for li, line in enumerate(fh):
                e = line.strip().split()
                seq_names.append(e[0].split('/')[-1][:-3])
                M[li,:(li+1)] = [float(x) for x in e[1:]]

        Msym = M + M.T
        Msym -= np.eye(nrows)

    return Msym, np.array(seq_names)