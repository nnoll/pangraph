import os, sys
import csv
import gzip
import numpy as np

from io   import StringIO
from enum import IntEnum

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ------------------------------------------------------------------------
# Helper/debugging functions

def cdfplot(x, **kwargs):
    import matplotlib.pylab as plt

    plt.plot(sorted(x), np.linspace(0, 1, len(x)), **kwargs)

def breakpoint(msg):
    print(f"DEBUG: {msg}", file=sys.stderr)
    import ipdb; ipdb.set_trace()

# ------------------------------------------------------------------------
# Global Enums

# TODO: find a better place for this to live
class Strand(IntEnum):
    Plus  = +1
    Minus = -1
    Null  = 0

def Complement(S):
    if isinstance(S, int):
        S = Strand(S)

    if S == Strand.Plus:
        return Strand.Minus
    elif S == Strand.Minus:
        return Strand.Plus
    elif S == Strand.Null:
        return Strand.Null
    else:
        raise ValueError(f"expected type 'Strand', got '{type(S)}'")

# Table from https://en.wikipedia.org/wiki/Nucleic_acid_notation
wcpair = {
    'A' : 'T',
    'T' : 'A',
    'C' : 'G',
    'G' : 'C',
    'W' : 'W',
    'S' : 'S',
    'M' : 'K',
    'K' : 'M',
    'R' : 'Y',
    'Y' : 'R',
    'B' : 'V',
    'D' : 'H',
    'H' : 'D',
    'V' : 'B',
    'Z' : 'Z',
}

def rev_cmpl(seq):
    if isinstance(seq,Seq):
        return str(Seq.reverse_complement(seq))
    else:
        return str(Seq.reverse_complement(Seq(seq)))

# ------------------------------------------------------------------------
# errors

def panic(msg):
    SystemExit(f"panic: {msg}")

def debug(msg):
    print(f"bug: {msg}", file=sys.stderr)
    import ipdb; ipdb.set_trace()

def warn(msg):
    raise print(f"warn: {msg}", file=sys.stderr)

def tryprint(msg, verbose):
    if verbose:
        print(msg)
    else:
        pass

def log(msg, file=sys.stderr):
    print(msg, file=file)

# ------------------------------------------------------------------------
# simple conversions

def as_array(x):
    return np.array(list(x))

def as_string(x):
    try:
        return x.view(f'U{x.size}')[0]
    except:
        return "".join(str(c) for c in x)

def flatten(x):
    return np.ndarray.flatten(x[:])

def cat(*args):
    return np.concatenate(tuple(arg for arg in args))

def as_record(seq, name):
    return SeqRecord(Seq(seq), id=name, name=name, description="")

# ------------------------------------------------------------------------
# file handling

# equivalent to mkdir -p
def mkdir(path):
    if os.path.exists(path):
        return os.path.isdir(path)
    try:
        os.makedirs(path)
        return 1
    except:
        return 0

def openany(path, mode='r'):
    if path.endswith('.gz'):
        if not mode.endswith('b'):
            mode += 't'
        return gzip.open(path, mode)
    else:
        return open(path, mode)

# ------------------------------------------------------------------------
# misc

def new_strand(s, t):
    if not isinstance(s, Strand) or not isinstance(t, Strand):
        raise TypeError(f"Expected an enum! Recieved {type(t)} and {type(s)}")

    if s != t:
        return Strand.Minus
    else:
        return Strand.Plus

def getnwk(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "%s:%.8f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.8f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getnwk(node.get_left(), newick, node.dist, leaf_names)
        newick = getnwk(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick

def as_str(s):
    if isinstance(s, bytes):
        return s.decode('utf-8')
    return s

# ------------------------------------------------------------------------
# parsers

def parse_fasta(fh):
    class Record:
        def __init__(self, name=None, meta=None, seq=None):
            self.seq  = seq
            self.name = name
            self.meta = meta

        def __str__(self):
            NL = '\n'
            nc = 80
            return f">{self.name} {self.meta}\n{NL.join([self.seq[i:(i+nc)] for i in range(0, len(self.seq), nc)])}"

        def __repr__(self):
            return str(self)

    header = as_str(fh.readline())
    while header != "" and header[0] == ">":
        name = header[1:].split()
        seq  = StringIO()
        for line in fh:
            line = as_str(line)
            if line == "" or line[0] == ">":
                break
            seq.write(line[:-1])

        header = as_str(line)
        yield Record(name=name[0], meta=" ".join(name[1:]), seq=seq.getvalue())
        seq.close()

def parse_paf(fh):
    hits = []
    for line in fh:
        row = line.strip().split()
        hit = {'qry': {'name'    : row[0],
                       'len'     : int(row[1]),
                       'start'   : int(row[2]),
                       'end'     : int(row[3])},
               'ref': {'name'    : row[5],
                       'len'     : int(row[6]),
                       'start'   : int(row[7]),
                       'end'     : int(row[8])},
               'aligned_bases'   : int(row[9]),
               'aligned_length'  : int(row[10]),
               'mapping_quality' : int(row[11]),
               'orientation'     : Strand.Plus if row[4]=='+' else Strand.Minus}
        for xtra in row[12:]:
            if xtra.startswith('cg:'):
                hit['cigar'] = xtra.split(':')[-1]
            elif xtra.startswith('de:f'):
                hit['divergence'] = float(xtra.split(':')[-1])
            elif xtra.startswith('AS:i'):
                hit['align_score'] = int(xtra.split(":")[2]) / hit['aligned_length']

        hits.append(hit)
    return hits

def parse_cigar(aln, qryseq, refseq, cutoff=500):
    from cigar import Cigar

    aln = Cigar(aln)

    lq, rq = 0, 0
    lr, rr = 0, 0
    refs = []
    qrys = []
    blks = []

    R, Q = {}, {}
    blkseq = ""
    blkpos = 0
    refmap = [(rr, blkpos-rr)]
    qrymap = [(rq, blkpos-rq)]

    def push(qval=None, rval=None):
        nonlocal R, Q, blkseq, blkpos, refmap, qrymap
        assert not (qval is None and rval is None)

        def f(xs, x):
            if x is None:
                xs.append(None)
                return True
            else:
                l, r = zip(x)
                if l < r:
                    xs.append(x)
                    return True
                return False

        hasq = f(qrys, qval)
        hasr = f(refs, rval)

        if hasq or hasr:
            assert len(qrys) == len(refs)
            assert len(blkseq) > 0, "empty seq"
            blks.append((np.array(list(blkseq)), (Q, np.array(qrymap).T), (R, np.array(refmap).T)))

        R, Q = {}, {}
        blkseq = ""
        blkpos = 0
        refmap = [(rr, blkpos-rr)]
        qrymap = [(rq, blkpos-rq)]

    def recordbp():
        nonlocal blkpos, refmap, qrymap

        blkpos = len(blkseq)
        refmap.append((rr, blkpos-rr))
        qrymap.append((rq, blkpos-rq))

    for l, t in aln.items():
        if t in ['S', 'H']:
            if l >= cutoff:
                print(aln)
                import ipdb; ipdb.set_trace()

                push((lq, rq), (lr, rr))

                blkseq = qryseq[rq:rq+l]
                # TODO: Think through soft/hard clips
                # if t == 'S':
                rq += l
                recordbp()

                push((rq-l, rq), None)
                lq = rq
                lr = rr
            else:
                rq += l
                recordbp()

        elif t == 'M':
            rs = np.array(list(refseq[rr:rr+l]))
            qs = np.array(list(qryseq[rq:rq+l]))
            if len(rs) != len(qs):
                import ipdb; ipdb.set_trace()
            diff = np.where(np.array(rs != qs))[0]
            for i in diff:
                Q[i+blkpos] = qs[i]
            blkseq += refseq[rr:rr+l]

            rq += l
            rr += l

            recordbp()

        elif t == 'D':
            if l >= cutoff:
                push((lq, rq), (lr, rr))

                blkseq = refseq[rr:rr+l]
                rr += l
                recordbp()

                push(None, (rr-l, rr))
                lr = rr
                lq = rq
            else:
                for i in range(l):
                    Q[i+blkpos] = '-'
                blkseq += refseq[rr:rr+l]

                rr += l
                recordbp()

        elif t == 'I':
            if l >= cutoff:
                push((lq, rq), (lr, rr))

                blkseq = qryseq[rq:rq+l]
                rq += l
                recordbp()

                push((rq-l, rq), None)
                lq = rq
                lr = rr
            else:
                for i in range(l):
                    R[i+blkpos] = '-'
                blkseq += qryseq[rq:rq+l]

                rq += l
                recordbp()

    push((lq, rq), (lr, rr))
    assert len(qrys) == len(refs) and len(qrys) == len(blks)

    return qrys, refs, blks
