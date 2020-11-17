import io, os, sys
import json
import numpy as np
import pprint
import subprocess
import tempfile

from io          import StringIO
from glob        import glob
from collections import defaultdict, Counter
from itertools   import chain
from copy        import deepcopy

from Bio           import AlignIO, SeqIO, Phylo
from Bio.Seq       import Seq
from Bio.SeqRecord import SeqRecord

from scipy.stats import entropy

from .         import suffix
from .block    import Block
from .sequence import Node, Path
from .utils    import Strand, as_string, parse_paf, panic, as_record, new_strand, breakpoint

# ------------------------------------------------------------------------
# globals

pp = pprint.PrettyPrinter(indent=4)

# ------------------------------------------------------------------------
# utility

def alignment_entropy(rdr):
    try:
        aln = np.array([list(rec) for rec in AlignIO.read(rdr, 'fasta')], np.character).view(np.uint8)
        S   = sum(entropy(np.bincount(aln[:,i])/aln.shape[0]) for i in range(aln.shape[1]))
        return S/aln.shape[1]
    except Exception as msg:
        print(f"ERROR: {msg}")
        return None

# ------------------------------------------------------------------------
# Junction class
# simple struct

class Junction(object):
    def __init__(self, left, right):
        self.left  = left
        self.right = right

    def __eq__(self, other):
        if self.data == other.data:
            return True
        elif self.data == other.reverse().data:
            return False
        else:
            return False

    def __hash__(self):
        return hash(frozenset([self.data, self.reverse().data]))

    def __str__(self):
        return f"({self.left}, {self.right})"

    def __repr__(self):
        return str(self)

    @property
    def data(self):
        return ((self.left.blk.id, self.left.strand), (self.right.blk.id, self.right.strand))

    @property
    def right_id(self):
        return self.right.blk.id

    @property
    def left_id(self):
        return self.left.blk.id

    @property
    def left_blk(self):
        return (self.left.blk.id, self.left.strand)

    @property
    def right_blk(self):
        return (self.right.blk.id, self.right.strand)

    def reverse(self):
        return Junction(
            Node(self.right.blk, self.right.num, Strand(-1*self.right.strand)),
            Node(self.left.blk,  self.left.num,  Strand(-1*self.left.strand)),
        )

def rev_blk(b):
    return (b[0], Strand(-1*b[1]))

# ------------------------------------------------------------------------
# Graph class

class Graph(object):
    """docstring for Graph"""

    def __init__(self):
        self.name = ""   # The name of graph. Will be used as file basename in exports
        self.blks = {}   # All blocks/alignments
        self.seqs = {}   # All sequences (as list of blocks)
        self.sfxt = None # Suffix tree of block records
        self.dmtx = None # Graph distance matrix

    # --- Class methods ---

    @classmethod
    def from_seq(cls, name, seq, circular):
        newg = cls()
        blk  = Block.from_seq(name, seq)
        newg.name = name
        newg.blks = {blk.id : blk}
        newg.seqs = {name : Path(name, Node(blk, 0, Strand.Plus), 0, circular)}

        return newg

    @classmethod
    def from_dict(cls, d):
        G = Graph()
        G.name = d['name']
        G.blks = [Block.from_dict(b) for b in d['blocks']]
        G.blks = {b.id : b for b in G.blks}
        G.seqs = [Path.from_dict(seq, G.blks) for seq in d['seqs']]
        G.seqs = {p.name : p for p in G.seqs}
        G.sfxt = None
        G.dmtx = None
        if d['suffix'] is not None:
            G.compile_suffix()
            G.dmtx = d['distmtx']

        return G

    @classmethod
    def connected_components(cls, G):
        # -----------------------------
        # internal functions
        def overlaps(s1, s2):
            return len(s1.intersection(s2)) > 0
        def component(graph, name):
            cc = Graph()
            cc.blks = {id:G.blks.pop(id) for id in graph}
            cc.seqs = {nm:G.seqs.pop(nm) for nm in name}
            cc.sfxt = None
            cc.dmtx = None
            return cc

        # -----------------------------
        # main body
        graphs, names = [], []
        for name, path in G.seqs.items():
            blks = set([b.id for b in path.blocks()])
            gi   = [i for i, g in enumerate(graphs) if overlaps(blks, g)]
            if len(gi) == 0:
                graphs.append(blks)
                names.append(set([name]))
                continue

            graphs[gi[0]] = graphs[gi[0]].union(blks, *(graphs.pop(i) for i in gi[:0:-1]))
            names[gi[0]]  = names[gi[0]].union(set([name]), *(names.pop(i) for i in gi[:0:-1]))

        return [component(graph, name) for graph, name in zip(graphs, names)]

    @classmethod
    def fuse(cls, g1, g2):
        ng = Graph()
        combine = lambda d1, d2: {**d1, **d2}
        ng.blks = combine(g1.blks, g2.blks)
        ng.seqs = combine(g1.seqs, g2.seqs)

        return ng

    # ---------------
    # methods

    def marginalize(self, *isolates):
        G = Graph()
        isolates = set(isolates)
        G.seqs   = {iso:deepcopy(path) for iso, path in self.seqs.items() if iso in isolates}
        G.blks   = {blk.id:blk.marginalize(*isolates) for blk in set.union(*[p.blocks() for p in G.seqs.values()])}
        G.remove_transitives()
        return G

    def union(self, qpath, rpath, out, cutoff=0, alpha=10, beta=2, extensive=False, edge_window=1000, edge_extend=2500):
        from seqanpy import align_global as align

        # ----------------------------------
        # internal functions

        def energy(hit):
            l    = hit["aligned_bases"]
            if l <= cutoff:
                return l

            num  = lambda k: len(self.blks[hit[k]["name"]].muts)
            cuts = lambda k: (hit[k]['start'] > cutoff) + ((hit[k]['len']-hit[k]['end']) > cutoff)

            if extensive:
                delP = num('qry')*cuts('qry') + num('ref')*cuts('ref')
            else:
                delP = cuts('qry') + cuts('ref')
            dmut = hit["aligned_length"] * hit["divergence"]

            return -l + alpha*delP + beta*dmut

        def accepted(hit):
            return energy(hit) < 0

        if cutoff <= 0:
            def proc(hit):
                return hit
        else:
            def proc(hit):
                # -----------------------
                # load in sequences

                with open(f"{qpath}.fa", 'r') as fd:
                    qfa = {s.id:str(s.seq) for s in SeqIO.parse(fd, 'fasta')}

                if qpath == rpath:
                    rfa = qfa
                else:
                    with open(f"{rpath}.fa", 'r') as fd:
                        rfa = {s.id:str(s.seq) for s in SeqIO.parse(fd, 'fasta')}

                # -----------------------
                # internal functions

                def to_cigar(aln):
                    cigar = ""
                    s1, s2 = np.fromstring(aln[0], dtype=np.int8), np.fromstring(aln[1], dtype=np.int8)
                    M, I, D = 0, 0, 0
                    for (c1, c2) in zip(s1, s2):
                        if c1 == ord("-") and c2 == ord("-"):
                            breakpoint("panic")
                        elif c1 == ord("-"):
                            if I > 0:
                                cigar += f"{I}I"
                                I = 0
                            elif M > 0:
                                cigar += f"{M}M"
                                M = 0
                            D += 1
                        elif c2 == ord("-"):
                            if D > 0:
                                cigar += f"{D}D"
                                D = 0
                            elif M > 0:
                                cigar += f"{M}M"
                                M = 0
                            I += 1
                        else:
                            if D > 0:
                                cigar += f"{D}D"
                                D = 0
                            elif I > 0:
                                cigar += f"{I}I"
                                I = 0
                            M += 1
                    if I > 0:
                        cigar += f"{I}I"
                        I = 0
                    elif M > 0:
                        cigar += f"{M}M"
                        M = 0
                    elif D > 0:
                        cigar += f"{D}D"
                        M = 0

                    return cigar

                def revcmpl_if(s, cond):
                    if cond:
                        return str(Seq.reverse_complement(Seq(s)))
                    else:
                        return s

                def get_seqs():
                    return qfa[hit['qry']['name']], rfa[hit['ref']['name']]

                # -----------------------
                # body

                dS_q = hit['qry']['start']
                dE_q = hit['qry']['len'] - hit['qry']['end']
                dS_r = hit['ref']['start']
                dE_r = hit['ref']['len'] - hit['ref']['end']

                # Left side of match
                if 0 < dS_q <= cutoff and (dS_r > cutoff or dS_r == 0):
                    hit['cigar'] = f"{dS_q}I" + hit['cigar']
                    hit['qry']['start'] = 0
                elif 0 < dS_r <= cutoff and (dS_q > cutoff or dS_q == 0):
                    hit['cigar'] = f"{dS_r}D" + hit['cigar']
                    hit['ref']['start'] = 0
                elif 0 < dS_q <= cutoff and 0 < dS_r <= cutoff:
                    qseq, rseq = get_seqs()
                    aln = align(revcmpl_if(qseq, hit['orientation']==Strand.Minus)[0:dS_q], rseq[0:dS_r])[1:]

                    hit['cigar'] = to_cigar(aln) + hit['cigar']
                    hit['qry']['start'] = 0
                    hit['ref']['start'] = 0
                    hit['aligned_bases'] += len(aln[0])

                # Right side of match
                if 0 < dE_q <= cutoff and (dE_r > cutoff or dE_r == 0):
                    hit['cigar'] += f"{dE_q}I"
                    hit['qry']['end'] = hit['qry']['len']
                elif 0 < dE_r <= cutoff and (dE_q > cutoff or dE_q == 0):
                    hit['cigar'] += f"{dE_r}D"
                    hit['ref']['end'] = hit['ref']['len']
                elif 0 < dE_q <= cutoff and 0 < dE_r <= cutoff:
                    qseq, rseq = get_seqs()
                    aln = align(revcmpl_if(qseq, hit['orientation']==Strand.Minus)[-dE_q:], rseq[-dE_r:])[1:]

                    hit['cigar'] = hit['cigar'] + to_cigar(aln)
                    hit['qry']['end'] = hit['qry']['len']
                    hit['ref']['end'] = hit['ref']['len']
                    hit['aligned_bases'] += len(aln[0])

                return hit

        # ----------------------------------
        # body

        os.system(f"minimap2 -t 2 -x asm20 -m 10 -n 2 -s 30 -D -c {rpath}.fa {qpath}.fa 1>{out}.paf 2>log")

        with open(f"{out}.paf") as fd:
            paf = parse_paf(fd)
        paf.sort(key = lambda x: energy(x))

        merged_blks = set()
        if len(paf) == 0:
            return self, False

        merged = False
        for hit in paf:
            if hit['qry']['name'] in merged_blks \
            or hit['ref']['name'] in merged_blks \
            or(hit['ref']['name'] <= hit['qry']['name'] and qpath == rpath) \
            or not accepted(hit):
                continue

            merged = True
            self.merge(proc(hit), edge_window, edge_extend)
            merged_blks.add(hit['ref']['name'])
            merged_blks.add(hit['qry']['name'])

        for b in self.blks.values():
            L = len(b.seq)
            for mut in b.muts.values():
                for x in mut.keys():
                    if x >= L:
                        import ipdb; ipdb.set_trace()

        self.remove_transitives()

        for path in self.seqs.values():
            path.rm_nil_blks()

        return self, merged

    # a junction is a pair of adjacent blocks.
    def junctions(self):
        junctions = defaultdict(list)
        for iso, path in self.seqs.items():
            if len(path.nodes) == 1:
                continue

            for i, n in enumerate(path.nodes):
                if path.circular or i > 0:
                    j = Junction(path.nodes[i-1], n)
                    junctions[j].append(iso)
        return { k:dict(Counter(v)) for k, v in junctions.items() }

    def remove_transitives(self):
        js = self.junctions()
        transitives = []
        for j, isos in js.items():
            left_eq_right = self.blks[j.left.blk.id].isolates == self.blks[j.right.blk.id].isolates
            left_eq_isos  = isos == self.blks[j.left.blk.id].isolates
            if left_eq_right and left_eq_isos:
                transitives.append(j)

        chains = {}
        for j in transitives:
            if j.left_id in chains and j.right_id in chains:
                c1, c2 = chains[j.left_id], chains[j.right_id]
                if c1 == c2:
                    continue

                if j.left_blk==c1[-1] and j.right_blk==c2[0]:
                    new_chain = c1 + c2
                elif j.left_blk==c1[-1] and rev_blk(j.right_blk)==c2[-1]:
                    new_chain = c1 + [rev_blk(b) for b in c2[::-1]]
                elif rev_blk(j.left_blk)==c1[0] and j.right_blk==c2[0]:
                    new_chain = [rev_blk(b) for b in c1[::-1]] + c2
                elif rev_blk(j.left_blk)==c1[0] and rev_blk(j.right_blk)==c2[-1]:
                    new_chain = c2 + c1
                else:
                    breakpoint("case not covered")

                for b, _ in new_chain:
                    chains[b] = new_chain

            elif j.left_id in chains:
                c = chains[j.left_id]
                if j.left_blk == c[-1]:
                    c.append(j.right_blk)
                elif rev_blk(j.left_blk) == c[0]:
                    c.insert(0, rev_blk(j.right_blk))
                else:
                    breakpoint("chains should be linear")
                chains[j.right_id] = c
            elif j.right_id in chains:
                c = chains[j.right_id]
                if j.right_blk == c[-1]:
                    c.append(rev_blk(j.left_blk))
                elif j.right_blk == c[0]:
                    c.insert(0, j.left_blk)
                else:
                    breakpoint("chains should be linear")
                chains[j.left_id] = c
            else:
                chains[j.left_id]  = [j.left_blk, j.right_blk]
                chains[j.right_id] = chains[j.left_id]

        chains = list({id(c):c for c in chains.values()}.values())

        for i, c in enumerate(chains):
            new_blk = Block.cat([self.blks[b] if s == Strand.Plus else self.blks[b].rev_cmpl() for b, s in c])

            # TODO: check that isos is constant along the chain
            for iso in self.blks[c[0][0]].isolates.keys():
                self.seqs[iso].merge(c[0], c[-1], new_blk)

            new_blk.muts = {key:val for key,val in new_blk.muts.items() if isinstance(key,tuple)}
            # if set(self.blks[c[0][0]].isolates.keys()) != set(new_blk.isolates.keys()):
            #     breakpoint("bad mutation keys")

            self.blks[new_blk.id] = new_blk

            for b, _ in c:
                self.blks.pop(b)

    def prune_blks(self):
        blks = set()
        for path in self.seqs.values():
            blks.update(path.blocks())
        self.blks = {b.id:self.blks[b.id] for b in blks}

    def merge(self, hit, window, extend):
        old_ref = self.blks[hit['ref']['name']]
        old_qry = self.blks[hit['qry']['name']]

        # As we slice here, we DONT need to remember the starting position.
        # This is why in from_aln(aln) we set the start index to 0
        ref = old_ref[hit['ref']['start']:hit['ref']['end']]
        qry = old_qry[hit['qry']['start']:hit['qry']['end']]

        if hit["orientation"] == Strand.Minus:
            qry = qry.rev_cmpl()

        aln = {"ref_seq"     : ref.sequence,
               "qry_seq"     : qry.sequence,
               "cigar"       : hit["cigar"],
               "ref_cluster" : ref.muts,
               "qry_cluster" : qry.muts,
               "ref_start"   : hit["ref"]["start"],
               "ref_name"    : hit["ref"]["name"],
               "qry_start"   : hit["qry"]["start"],
               "qry_name"    : hit["qry"]["name"],
               "orientation" : hit["orientation"]}

        merged_blks, new_qrys, new_refs, shared_blks, blk_map = Block.from_aln(aln)

        for merged_blk in merged_blks:
            self.blks[merged_blk.id] = merged_blk

        def update(blk, add_blks, hit, strand):
            new_blks = []

            # The convention for the tuples are (block, strand orientation, merged)
            if hit['start'] > 0:
                left = blk[0:hit['start']]
                self.blks[left.id] = left
                new_blks.append((left, Strand.Plus, False))

            for b in add_blks:
                new_blks.append((b, strand, True))

            if hit['end'] < len(blk):
                right = blk[hit['end']:]
                self.blks[right.id] = right
                new_blks.append((right, Strand.Plus, False))

            for tag in blk.muts.keys():
                path = self.seqs[tag[0]]
                path.replace(blk, tag, new_blks, blk_map)

            return new_blks

        # NOTE: This has repeated blocks...
        new_blocks = []
        new_blocks.extend(update(old_ref, new_refs, hit['ref'], Strand.Plus))
        new_blocks.extend(update(old_qry, new_qrys, hit['qry'], hit['orientation']))

        # lblks_set_x, rblks_set_x = set(), set()
        # lblks_set_s, rblks_set_s = set(), set()
        # first    = True
        # num_seqs = 0
        # for tag in shared_blks[0].muts.keys():
        #     pos    = [self.seqs[tag[0]].position_of(b, tag[1]) for b in shared_blks]
        #     strand = [self.seqs[tag[0]].orientation_of(b, tag[1]) for b in shared_blks]
        #     beg, end = pos[0], pos[-1]
        #     if strand[0] == Strand.Plus:
        #         lwindow = min(window, shared_blks[0].len_of(*tag))
        #         rwindow = min(window, shared_blks[-1].len_of(*tag))

        #         lblks_x = self.seqs[tag[0]][beg[0]-extend:beg[0]+lwindow]
        #         rblks_x = self.seqs[tag[0]][end[1]-rwindow:end[1]+extend]
        #     elif strand[0] == Strand.Minus:
        #         lwindow = min(window, shared_blks[-1].len_of(*tag))
        #         rwindow = min(window, shared_blks[0].len_of(*tag))

        #         rblks_x = self.seqs[tag[0]][beg[0]-extend:beg[0]+rwindow]
        #         lblks_x = self.seqs[tag[0]][end[1]-lwindow:end[1]+extend]
        #     else:
        #         raise ValueError("unrecognized strand polarity")

        #     if first:
        #         lblks_set_x = set([b.id for b in lblks_x])
        #         rblks_set_x = set([b.id for b in rblks_x])

        #         lblks_set_s = set([b.id for b in lblks_x])
        #         rblks_set_s = set([b.id for b in rblks_x])

        #         first = False
        #     else:
        #         lblks_set_x.intersection_update(set([b.id for b in lblks_x]))
        #         rblks_set_x.intersection_update(set([b.id for b in rblks_x]))

        #         lblks_set_s.update(set([b.id for b in lblks_x]))
        #         rblks_set_s.update(set([b.id for b in rblks_x]))
        #     num_seqs += 1

        # def emit(side):
        #     if side == 'left':
        #         delta  = len(lblks_set_s)-len(lblks_set_x)
        #     elif side == 'right':
        #         delta = len(lblks_set_s)-len(rblks_set_x)
        #     else:
        #         raise ValueError(f"unrecognized argument '{side}' for side")

        #     if delta > 0 and num_seqs > 1:
        #         print(f">LEN={delta}", end=';')
        #         try:
        #             fd, path = tempfile.mkstemp()
        #             with os.fdopen(fd, 'w') as tmp:
        #                 for i, tag in enumerate(merged_blks[0].muts.keys()):
        #                     pos    = [self.seqs[tag[0]].position_of(b, tag[1]) for b in shared_blks]
        #                     strand = [self.seqs[tag[0]].orientation_of(b, tag[1]) for b in shared_blks]
        #                     beg, end = pos[0], pos[-1]

        #                     if strand[0] == Strand.Plus:
        #                         if side == 'left':
        #                             left, right = beg[0]-extend,beg[0]+min(window,shared_blks[0].len_of(*tag))
        #                         elif side == 'right':
        #                             left, right = end[1]-min(window,shared_blks[-1].len_of(*tag)),end[1]+extend
        #                         else:
        #                             raise ValueError(f"unrecognized argument '{side}' for side")

        #                     elif strand[0] == Strand.Minus:
        #                         if side == 'left':
        #                             left, right = end[1]-min(window,shared_blks[-1].len_of(*tag)),end[1]+extend
        #                         elif side == 'right':
        #                             left, right = beg[0]-extend,beg[0]+min(window, shared_blks[0].len_of(*tag))
        #                         else:
        #                             raise ValueError(f"unrecognized argument '{side}' for side")

        #                     iso_blks = self.seqs[tag[0]][left:right]
        #                     tmp.write(f">isolate_{i:04d} {','.join(b.id for b in iso_blks)}\n")
        #                     s = self.seqs[tag[0]].sequence_range(left,right)
        #                     if len(s) > extend + window:
        #                         breakpoint(f"bad sequence slicing: {len(s)}")
        #                     tmp.write(s + '\n')

        #                 tmp.flush()

        #                 proc = subprocess.Popen(f"mafft --auto {path}",
        #                             stdout=subprocess.PIPE,
        #                             stderr=subprocess.PIPE,
        #                             shell=True)
        #                 out, err = proc.communicate()
        #                 print(f"ALIGNMENT={out}", end=";")
        #                 rdr = StringIO(out.decode('utf-8'))
        #                 print(f"SCORE={alignment_entropy(rdr)}", end=";")
        #                 rdr.close()
        #                 print("\n", end="")
        #         finally:
        #             os.remove(path)
        #     else:
        #         print(f">NO MATCH")

        # emit('left')
        # emit('right')
        self.prune_blks()

        return [b[0] for b in new_blocks]

    def extract(self, name, strip_gaps=True, verbose=False):
        seq = self.seqs[name].sequence()
        if strip_gaps:
            seq = seq.replace('-', '')
        return seq

    def compress_ratio(self, extensive=False, name=None):
        unc = 0
        if name is None:
            for n in self.seqs:
                seq  = self.extract(n)
                unc += len(seq)
            cmp = np.sum([len(x) for x in self.blks.values()])
        else:
            cmp = np.sum([len(x) for x in self.blks.values() if name in x.muts])

        return unc/cmp/len(self.seqs) if not extensive else unc/cmp

    def contains(self, other):
        return set(other.seqs.keys()).issubset(set(self.seqs.keys()))

    def pairwise_distance(self, weighted=False):
        strings     = {iso: [(n.blk.id,n.strand) for n in path.nodes] for iso,path in self.seqs.items()}
        suffix_tree = suffix.Tree(strings)
        isos        = sorted(list(self.seqs.keys()))
        N           = len(self.seqs)

        if not weighted:
            def kernel():
                D = np.zeros((N,N))
                for i, iso1 in enumerate(isos):
                    for j, iso2 in enumerate(isos[:i]):
                        matches    = suffix_tree.matches(iso1, iso2)
                        num_events = len(matches)
                        if num_events > 0:
                            D[i,j] = num_events
                        else:
                            D[i,j] = np.infty
                        D[j,i] = D[i,j]
                return D
        else:
            def kernel():
                D = np.zeros((N,N))
                for i, iso1 in enumerate(isos):
                    L1 = len(self.seqs[iso1])
                    for j, iso2 in enumerate(isos[:i]):
                        L2 = len(self.seqs[iso2])

                        matches = suffix_tree.matches(iso1, iso2)
                        length  = lambda name: sum(self.blks[b[0]].len_of(name, 0) for m in matches for b in m)
                        f1, f2  = length(iso1)/L1, length(iso2)/L2

                        D[i,j] = np.sqrt(f1*f2)
                        D[j,i] = D[i,j]
                return D

        return (kernel(), isos)

    def to_json(self, wtr, minlen=500):
        J = {}
        cleaned_seqs = {s:[b for b in self.seqs[s] if len(self.blks[b[0]])>minlen]
                             for s in self.seqs}
        relevant_blocks = set()
        for s in cleaned_seqs.values():
            relevant_blocks.update([b[0] for b in s])

        J['Isolate_names'] = list(cleaned_seqs.keys())
        J['Plasmids']      = [[x for x in cleaned_seqs[s]] for s in J['Isolate_names']]

        # Build node substructure
        nodes = {}
        for b in relevant_blocks:
            aln = { J["Isolate_names"].index(iso) :
                    self.blks[b].extract(iso, num, strip_gaps=False) for iso, num in self.blks[b].muts }
            nodes[b] = {"ID"        : b,
                        "Genomes"   : {"Consensus" : self.blks[b].sequence,
                                       "Alignment" : aln },
                        "Out_Edges" : [],
                        "In_Edges"  : []}

        # Gather edges (node/node junctions) for each node
        edges = {}
        for pname, p in zip(range(len(J["Isolate_names"])), J['Plasmids']):
            for i in range(len(p)-1):
                e = (p[i], p[i+1])
                if e in edges:
                    edges[e]["Isolates"].append(pname)
                else:
                    edges[e] = {"Source" : e[0], "Target" : e[1], "Isolates" : [pname]}
            e = (p[-1], p[0])
            if e in edges:
                edges[e]["Isolates"].append(pname)
            else:
                edges[e] = {"Source" : e[0], "Target" : e[1], "Isolates" : [pname]}

        for e in edges:
            nodes[e[0][0]]["Out_Edges"].append(edges[e])
            nodes[e[1][0]]["In_Edges"].append(edges[e])

        J["Nodes"] = nodes
        json.dump(J, wtr)

    def to_dict(self):
        return {'name'   : self.name,
                'seqs'   : [s.to_dict() for s in self.seqs.values()],
                'blocks' : [b.to_dict() for b in self.blks.values()],
                'suffix' : None if self.sfxt is None else "compiled",
                'distmtx': self.dmtx}

    def write_fasta(self, wtr):
        SeqIO.write(sorted([ as_record(c.sequence, c.id) for c in self.blks.values() ],
            key=lambda x: len(x), reverse=True), wtr, format='fasta')
