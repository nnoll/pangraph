import sys
import numpy as np

from collections import ChainMap

from .utils import Strand, log, breakpoint, new_strand, rev_cmpl

# ------------------------------------------------------------------------
# Node class: one visit along a path

class Node(object):
    """docstring for Node"""

    def __init__(self, blk, num, strand):
        self.blk    = blk
        self.num    = num
        self.strand = strand

    def __str__(self):
        return f"({self.blk}, {self.num}, {self.strand})"

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return self.blk.id == other.blk.id and self.strand == other.strand

    def __hash__(self):
        return hash((self.blk.id, self.strand))

    @classmethod
    def from_dict(cls, d, blks):
        return Node(blks[d['id']], d['num'], Strand(d['strand']))

    def to_dict(self):
        return {'id': self.blk.id, 'num': self.num, 'strand': int(self.strand)}

    def length(self, name):
        return self.blk.len_of(name, self.num)

# ------------------------------------------------------------------------
# Path class: (circular) list of blocks

class Path(object):
    """docstring for Path"""

    def __init__(self, name, nodes, offset, circular):
        super(Path, self).__init__()
        self.name     = name
        self.nodes    = nodes if isinstance(nodes, list) else [nodes]
        self.offset   = offset
        self.position = np.cumsum([0] + [n.length(name) for n in self.nodes])
        self.circular = circular

        # if offset > 0 and not circular:
        #     raise ValueError("sequence path cannot have non-zero offset if it corresponds to linear genome")

    def __str__(self):
        return f"{self.name}: {[str(n) for n in self.nodes]}"

    def __repr__(self):
        return str(self)

    @classmethod
    def from_dict(cls, d, blks):
        if 'circular' in d:
            return Path(d['name'], [Node.from_dict(n, blks) for n in d['nodes']], d['offset'], d['circular'])
        return Path(d['name'], [Node.from_dict(n, blks) for n in d['nodes']], d['offset'], False)

    def to_dict(self):
        return {'name': self.name, 'offset': self.offset, 'nodes': [n.to_dict() for n in self.nodes], 'circular': self.circular}

    def blocks(self):
        return set([n.blk for n in self.nodes])

    def sequence(self, verbose=False):
        seq = ""
        for n in self.nodes:
            s = n.blk.extract(self.name, n.num, strip_gaps=True, verbose=verbose)
            if n.strand == Strand.Plus:
                seq += s
            else:
                seq += rev_cmpl(s)

        if self.offset != 0:
            # if not self.circular:
            #     raise ValueError("invalid sequence path: non-zero offset for linear genome")
            seq = seq[self.offset:] + seq[:self.offset]

        return seq

    def rm_nil_blks(self):
        good, popped = [], set()
        for i, n in enumerate(self.nodes):
            if n.blk.id in popped:
                continue

            if n.blk.is_empty(self.name, n.num):
                if (self.name, n.num) not in n.blk.muts:
                    breakpoint("malformed mutation bookkeeping!")
                n.blk.muts.pop((self.name, n.num))
            else:
                good.append(i)

            if not n.blk.has(self.name):
                popped.add(n.blk.id)

        self.nodes    = [self.nodes[i] for i in good]
        self.position = np.cumsum([0] + [n.length(self.name) for n in self.nodes])

    def merge(self, start, stop, new):
        N = 0
        while True:
            ids = [n.blk.id for n in self.nodes]
            try:
                i, j = ids.index(start[0]), ids.index(stop[0])

                if self.nodes[i].strand == start[1]:
                    beg, end, s = i, j, Strand.Plus
                else:
                    beg, end, s = j, i, Strand.Minus

                key = (self.name,N)

                if beg < end:
                    # print(f"----> case 1: {self.nodes[beg:end+1]} ({beg}, {end})")
                    # s0  = "".join(n.blk.extract(self.name, n.num) for n in self.nodes[beg:end+1])
                    val = dict(ChainMap(*[new.muts[n.blk][(self.name,n.num)] for n in self.nodes[beg:end+1]]))
                    new.muts.update({key:val})
                    self.nodes = self.nodes[:beg] + [Node(new, N, s)] + self.nodes[end+1:]

                    # s1 = new.extract(self.name,N)
                    # if s0 != s1:
                    #     breakpoint("bad fwd-ordered mutations")
                else:
                    # print(f"----> case 2: {self.nodes[beg:] + self.nodes[:end+1]} ({beg}, {end})")
                    # s0 = "".join(n.blk.extract(self.name, n.num) for n in self.nodes[beg:] + self.nodes[:end+1])
                    # if not self.circular:
                    #     raise ValueError("attempted to rotate non-circular sequence")
                    self.offset += sum(n.blk.len_of(self.name, N) for n in self.nodes[beg:])
                    val = dict(ChainMap(*[new.muts[n.blk][(self.name,n.num)] for n in self.nodes[beg:] + self.nodes[:end+1]]))
                    new.muts.update({key:val})
                    self.nodes   = [Node(new, N, s)] + self.nodes[end+1:beg]

                    # s1 = new.extract(self.name, N)
                    # if s0 != s1:
                    #     breakpoint("bad rev-ordered mutations")
                self.position  = np.cumsum([0] + [n.length(self.name) for n in self.nodes])
                N += 1
            except ValueError as err:
                print(f"Error: {err}")
                return

    def replace(self, blk, tag, new_blks, blk_map):
        new = []
        for n in self.nodes:
            if n.blk.id == blk.id and n.num == tag[1]:
                os  = n.strand
                mk  = lambda b,ns,merged: Node(b, blk_map[b.id][blk.id][tag][1], new_strand(os, ns)) if merged else Node(b, n.num, new_strand(os, ns))
                tmp = [mk(blk,ns,flag) for blk, ns, flag in new_blks]
                if os == Strand.Minus:
                    tmp = tmp[::-1]
                new.extend(tmp)
            else:
                new.append(n)

        self.nodes    = new
        self.position = np.cumsum([0] + [n.length(self.name) for n in self.nodes])

    def position_of(self, blk, num):
        index = { n.num:i for i, n in enumerate(self.nodes) if n.blk == blk }
        if not num in index:
            return None
        return (self.position[index[num]], self.position[index[num]+1])

    def orientation_of(self, blk, num):
        orientation = { n.num:n.strand for i, n in enumerate(self.nodes) if n.blk == blk }
        if not num in orientation:
            return None
        return orientation[num]

    # TODO: pull out common functionality into a helper function
    # TODO: merge with other sequence function
    def sequence_range(self, start=None, stop=None):
        beg = start or 0
        end = stop or self.position[-1]
        l, r = "", ""
        if beg < 0:
            if len(self.nodes) > 1 and self.circular:
                l = self.sequence_range(self.position[-1]+beg,self.position[-1])
            beg = 0
        if end > self.position[-1]:
            if len(self.nodes) > 1 and self.circular:
                r = self.sequence_range(0,end-self.position[-1])
            end = self.position[-1]
        if beg > end:
            beg, end = end, beg

        i = np.searchsorted(self.position, beg, side='right') - 1
        j = np.searchsorted(self.position, end, side='left')
        m = ""
        if i < j:
            if j >= len(self.position):
                breakpoint("what?")
            if i == j - 1:
                m = self.nodes[i].blk.extract(self.name, self.nodes[i].num)[(beg-self.position[i]):(end-self.position[i])]
            else:
                m = self.nodes[i].blk.extract(self.name, self.nodes[i].num)[(beg-self.position[i]):]
                for n in self.nodes[i+1:j-1]:
                    m += n.blk.extract(self.name, n.num)
                n  = self.nodes[j-1]
                b  = n.blk.extract(self.name, n.num)
                m += b[0:(end-self.position[j-1])]
        return l + m + r

    def __getitem__(self, index):
        if isinstance(index, slice):
            beg = index.start or 0
            end = index.stop or self.position[-1]
            l, r = [], []
            if beg < 0:
                if len(self.nodes) > 1 and self.circular:
                    l = self[(self.position[-1]+beg):self.position[-1]]
                beg = 0
            if end > self.position[-1] and self.circular:
                if len(self.nodes) > 1:
                    r = self[0:(end-self.position[-1])]
                end = self.position[-1]
            if beg > end:
                beg, end = end, beg

            i = np.searchsorted(self.position, beg, side='right') - 1
            j = np.searchsorted(self.position, end, side='left')
            if i > j:
                breakpoint(f"not sorted, {beg}-{end}")
            return l + [n.blk for n in self.nodes[i:j]] + r
        elif isinstance(index, int):
            i = np.searchsorted(self.position, index, side='right') - 1
            return self.nodes[i].blk
        else:
            raise ValueError(f"type '{type(index)}' not supported as index")

    def __len__(self):
        return self.position[-1]
