import numpy as np
import numpy.random as rng

from collections import defaultdict, Counter
from .utils import parse_cigar, wcpair, as_array, as_string

# ------------------------------------------------------------------------
# Helper functions

RS = rng.RandomState(0)
def randomid():
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    randomid.N += 1
    name = "".join([alphabet[i] for i in RS.randint(0, len(alphabet), 10)])

    return name
randomid.N = 0

# ------------------------------------------------------------------------
# Block class

class Block(object):
    """docstring for Block"""

    def __init__(self, gen=True):
        super(Block, self).__init__()
        self.id   = randomid() if gen else 0
        self.seq  = None
        self.muts = {}

    def __str__(self):
        return str(self.id)

    def __repr__(self):
        return str(self)

    # ------------------
    # properties

    @property
    def sequence(self):
        return "".join(n for n in self.seq if n != '-')

    @property
    def length(self):
        return len(self.sequence)

    @property
    def depth(self):
        return len(self.muts)

    @property
    def isolates(self):
        return dict(Counter([k[0] for k in self.muts]))

    # ------------------
    # static methods

    @classmethod
    def from_seq(cls, name, seq):
        new_blk      = cls()
        new_blk.seq  = as_array(seq)
        new_blk.muts = {(name, 0):{}}

        return new_blk

    @classmethod
    def from_dict(cls, d):
        def unpack(key):
            t = tuple(key.split("?###?"))
            return (t[0], int(t[1]))

        B      = Block()
        B.id   = d['id']
        B.seq  = as_array(d['seq'])
        B.muts = {unpack(k):{int(x):sub for x, sub in v.items()} for k, v in d['muts'].items()}

        return B

    @classmethod
    def cat(cls, blks):
        nblk = Block()
        assert all([blks[0].muts.keys() == b.muts.keys() for b in blks[1:]])

        nblk.seq  = np.concatenate([b.seq for b in blks])
        nblk.muts = { blks[0]: { s:dict(x) for s,x in blks[0].muts.items() }}
        offset    = len(blks[0])
        for b in blks[1:]:
            nblk.muts[b] = {}
            for s in b.muts.keys():
                nblk.muts[b][s] = {(x+offset):snp for x,snp in b.muts[s].items()}
            offset += len(b)

        return nblk

    @classmethod
    def from_aln(cls, aln, debug=False):
        def updatemuts(blk, xtramuts, xmap, omuts, ival):
            seq = blk.seq
            # Iterate over all sequences in the block
            isomap = {}
            for iso, muts in omuts.items():
                # Shift position by indel #'s
                opos = as_array(m for m in muts.keys() if m >= ival[0] and m < ival[1])
                npos = opos + xmap[1][np.searchsorted(xmap[0], opos, side='right')]

                assert (np.max(npos) if len(npos) > 0 else 0) < len(seq)

                newmuts = {np : muts[op] for np, op in zip(npos, opos)}
                for p, n in xtramuts.items():
                    if p in newmuts:
                        if newmuts[p] == seq[p]:
                            newmuts.pop(p)
                    else:
                        newmuts[p] = n

                isomap[iso] = blk.push(iso, newmuts)

            return blk, isomap

        qrys, refs, blks = parse_cigar(aln['cigar'], aln['qry_seq'], aln['ref_seq'])

        # Iterate over all merged blocks and merge their sequences + mutations.
        newblks = []
        isomap  = defaultdict(lambda: {})
        for i, blk in enumerate(blks):
            newblk      = cls()
            newblk.muts = {}

            newblk.seq, qry, ref = blk
            Q, qrymap = qry
            R, refmap = ref

            if qrys[i] is not None:
                newblk, isomap[newblk.id][aln['qry_name']] = updatemuts(newblk, Q, qrymap, aln['qry_cluster'], qrys[i])
            if refs[i] is not None:
                newblk, isomap[newblk.id][aln['ref_name']] = updatemuts(newblk, R, refmap, aln['ref_cluster'], refs[i])

            newblk.update_consensus()
            newblks.append(newblk)

        isomap = dict(isomap)

        qryblks = [nb for i, nb in enumerate(newblks) if qrys[i] is not None]
        if aln['orientation'] == -1:
            qryblks = qryblks[::-1]

        refblks    = [nb for i, nb in enumerate(newblks) if refs[i] is not None]
        sharedblks = [nb for i, nb in enumerate(newblks) if refs[i] is not None and qrys[i] is not None]

        return newblks, qryblks, refblks, sharedblks, isomap

    # --------------
    # methods

    def copy(self, keepname=False):
        b = Block()
        if keepname:
            b.name = self.name
        b.seq  = np.copy(self.seq)
        b.muts = {s:dict(self.muts[s]) for s in self.muts}

        return b

    def marginalize(self, *isolates):
        isolates = set(isolates)
        self.muts = {tag:m for tag, m in self.muts.items() if tag[0] in isolates}
        return self

    def extract(self, iso, num, strip_gaps=True, verbose=False):
        tag = (iso, num)
        tmp = np.copy(self.seq)
        for key, muts in self.muts.items():
            if not isinstance(key, tuple):
                continue
            for x in muts.keys():
                if x >= len(tmp):
                    import ipdb; ipdb.set_trace()

        for p, s in self.muts[tag].items():
            tmp[p] = s

        assert len(tmp) > 0, "empty sequence"

        if strip_gaps:
            return as_string(tmp[tmp != '-'])
        else:
            return as_string(tmp)

    def len_of(self, iso, num):
        tag    = (iso, num)
        length = len(self.seq)
        gaplen = sum(1 for s in self.muts[tag].values() if s == '-')
        return length - gaplen

    def is_empty(self, iso, num, strip_gaps=True):
        tag = (iso, num)
        seq = np.copy(self.seq)
        # NOTE: This is a hack. Need to investigate the error that arises.
        if tag not in self.muts:
            return True

        for p, s in self.muts[tag].items():
            seq[p] = s

        return len(seq) == 0 or all(nuc == "-" for nuc in seq)

    def rev_cmpl(self):
        from Bio import Seq

        nblk     = Block()
        nblk.seq = as_array(Seq.reverse_complement("".join(self.seq)))
        L        = len(self.seq)-1
        for s in self.muts:
            nblk.muts[s] = {L-p: wcpair.get(c,c) for p,c in self.muts[s].items()}

        return nblk

    def update_consensus(self):
        alleles = defaultdict(set)
        for tag, muts in self.muts.items():
            for locus, snp in muts.items():
                alleles[(locus,snp)].add(tag)

        all_tags = set(self.muts.keys())
        # NOTE: deal with case where site is not biallelic
        for (locus, snp), tags in alleles.items():
            if len(tags)/len(self.muts) > 0.5:
                new_snp = self.seq[locus]
                self.seq[locus] = snp

                # remove mutations in tag set
                for tag in tags:
                    self.muts[tag].pop(locus)

                # add mutation in inverse set
                for tag in all_tags - tags:
                    self.muts[tag][locus] = new_snp

    def push(self, iso, muts):
        tag = iso if isinstance(iso, tuple) else (iso, 0)

        while tag in self.muts:
            tag = (tag[0], tag[1]+1)

        self.muts[tag] = muts

        return tag

    def has(self, iso):
        for tag in self.muts.keys():
            if tag[0] == iso:
                return True
        return False

    def to_dict(self):
        def pack(key):
            return f"{key[0]}?###?{key[1]}"

        def fix(d):
            return {int(k):v for k, v in d.items()}

        return {'id'   : self.id,
                'seq'  : "".join(str(n) for n in self.seq),
                'muts' : {pack(k) : fix(v) for k, v in self.muts.items()}}

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, val):
        if isinstance(val, slice):
            b      = Block()
            offset = np.array([i for i,c in enumerate(self.seq) if c != '-'])
            start  = val.start or 0
            stop   = val.stop  or len(self.sequence)

            start, stop = offset[start], offset[stop-1]+1
            b.seq  = self.seq[start:stop]
            for s, _ in self.muts.items():
                b.muts[s] = {p-start:c for p,c in self.muts[s].items() if p>=start and p<stop}
            return b
        else:
            raise ValueError("item access not supported")
