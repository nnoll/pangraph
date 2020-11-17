import os, sys
import json

from math import inf
from copy import deepcopy

import numpy as np
import matplotlib.pylab as plt

from Bio.Seq import Seq

from .utils import Strand, log, flatten, panic, breakpoint, rev_cmpl
from .graph import Graph

# ------------------------------------------------------------------------
# Global variables

MAXSELFMAPS = 25

# ------------------------------------------------------------------------
# Helper functions

def nodiag(mtx):
    return mtx-np.diag(np.diag(mtx))

def parse(mtx):
    with open(mtx) as fh:
        nrows = int(fh.readline().strip())
        M, r  = np.zeros((nrows, nrows), dtype=float), 0

        del_idxs  = []
        seq_names = []
        for li, line in enumerate(fh):
            e = line.strip().split()
            n = e[0].split('/')[-1][:-3]
            if n not in seq_names:
                seq_names.append(n)
                M[li,:(li+1)] = [float(x) for x in e[1:]]
            else:
                del_idxs.append(li)

    M = np.delete(M, del_idxs, axis=0)
    M = np.delete(M, del_idxs, axis=1)

    # Symmetrize
    M = 1 - M
    M = nodiag(M + M.T)/2

    return M, seq_names

def to_list(dmtx):
    assert len(dmtx.shape) == 2 and dmtx.shape[0] == dmtx.shape[1], "expected a square matrix"

    dlst = []
    for n in range(dmtx.shape[0]):
        dlst.append(list(dmtx[n,:(n+1)]))

    return dlst

# ------------------------------------------------------------------------
# Clade and Tree classes

class Clade(object):
    # ---------------------------------
    # Internal functions

    def __init__(self, name, parent, dist, children=[]):
        self.name   = name
        self.dist   = dist
        self.parent = parent
        self.child  = children
        self.fapath = ""
        self.graph  = None

    def __str__(self):
        if self.dist is None:
            return f"{self.name} :: Unknown"
        else:
            return f"{self.name} :: {self.dist:.4f}"

    def __repr__(self):
        return self.__str__()

    # ---------------------------------
    # Static functions

    @classmethod
    def from_dict(cls, d, parent):
        N = Clade(d['name'], parent, d['dist'])
        N.child = [Clade.from_dict(child, N) for child in d['child']]
        N.fapath = d['fapath']
        if d['graph'] is not None:
            if isinstance(d['graph'], list):
                N.graph  = [Graph.from_dict(elt) for elt in d['graph']]
            else:
                N.graph  = Graph.from_dict(d['graph'])
        else:
            N.graph = None

        return N

    # ---------------------------------
    # Class methods

    def is_leaf(self):
        return len(self.child) == 0

    def postorder(self):
        for child in self.child:
            for it in child.postorder():
                yield it
        yield self

    def children_graphs(self, func=lambda n: n.graph):
        gs = []
        for child in self.child:
            if child.graph:
                gs.append(func(child))
            else:
                gs.extend(child.children_graphs(func))
        return gs

    def new_parent(self, parent, dist):
        self.parent = parent
        self.dist   = dist if dist > 0 else 0

    def to_nwk(self, wtr):
        if not self.is_leaf():
            wtr.write("(")
            for i, child in enumerate(self.child):
                if i > 0:
                    wtr.write(",")
                child.to_nwk(wtr)
            wtr.write(")")

        wtr.write(self.name)
        wtr.write(":")
        wtr.write(f"{self.dist:.6f}")

    def to_json(self):
        serialize = lambda gs: [g.to_dict() for g in gs] if isinstance(gs, list) else [gs.to_dict()]
        return {'name'     : self.name,
                'dist'     : self.dist,
                'child'    : [ child.to_json() for child in self.child ],
                'fapath'   : self.fapath,
                'graph'    : serialize(self.graph) if self.graph is not None else None }

    def set_level(self, level):
        for c in self.child:
            c.set_level(level+1)
        self.level = level

class Tree(object):
    # ------------------- 
    # Class constructor
    def __init__(self, bare=False):
        self.root   = Clade("ROOT", None, 0) if not bare else None
        self.seqs   = None
        self.leaves = None

    # ------------------- 
    # Static methods

    # Loading from json
    @classmethod
    def from_json(cls, rdr):
        data   = json.load(rdr)
        T      = Tree(bare=True)
        T.root = Clade.from_dict(data['tree'], None)

        leafs  = {n.name: n for n in T.get_leafs()}
        if data['seqs'] is not None:
            T.seqs = {leafs[k]:Seq(v) for k,v in data['seqs'].items()}

        return T

    # our own neighbor joining
    # Biopython implementation is WAY too slow.
    @classmethod
    def nj(cls, mtx, names, verbose=False):
        # -----------------------------
        # internal functions

        def q(D):
            n = D.shape[0]
            Q = (n-2)*D - (np.sum(D,axis=0,keepdims=True) + np.sum(D,axis=1,keepdims=True))
            np.fill_diagonal(Q, np.inf)
            return Q

        def minpair(q):
            i, j = np.unravel_index(np.argmin(q), q.shape)
            qmin = q[i, j]
            if i > j:
                i, j = j, i
            return (i, j), qmin

        def pairdists(D, i, j):
            n  = D.shape[0]
            d1 = .5*D[i,j] + 1/(2*(n-2)) * (np.sum(D[i,:], axis=0) - np.sum(D[j,:], axis=0))
            d2 = D[i,j] - d1

            # remove negative branches while keeping total fixed
            if d1 < 0:
                d2 -= d1
                d1  = 0
            if d2 < 0:
                d1 -= d2
                d2  = 0

            dnew = .5*(D[i,:] + D[j,:] - D[i, j])
            return d1, d2, dnew

        def join(D, debug=False):
            nonlocal idx
            Q = q(D)
            (i, j), qmin = minpair(Q)
            if debug:
                q0min = min(flatten(Q[:]))
                assert abs(qmin-q0min) < 1e-2, f"minimum not found correctly. returned {qmin}, expected {q0min}"
                print(f"{D}\n--> Joining {i} and {j}. d={D[i,j]}")

            node   = Clade(f"NODE_{idx:05d}", T.root, None, [T.root.child[i], T.root.child[j]])

            d1, d2, dnew = pairdists(D, i, j)
            node.child[0].new_parent(node, d1)
            node.child[1].new_parent(node, d2)

            D[i, :] = dnew
            D[:, i] = dnew
            D[i, i] = 0
            D = np.delete(D, j, axis=0)
            D = np.delete(D, j, axis=1)
            T.root.child[i] = node
            T.root.child.pop(j)

            idx = idx + 1

            return D

        # -----------------------------
        # body
        assert len(names) == len(set(names)), "non-unique names found"

        T = Tree()
        for name in names:
            T.root.child.append(Clade(name, T.root, None, children=[]))
        idx = 0

        while mtx.shape[0] > 2:
            if verbose:
                print(f"--> Matrix size={mtx.shape[0]}. Number of root child={len(T.root.child)}")
            mtx = join(mtx)

        assert mtx.shape[0] == 2
        d = mtx[0, 1]
        T.root.child[0].dist = d/2
        T.root.child[1].dist = d/2

        return T

    # ------------------- 
    # methods 

    def postorder(self):
        return self.root.postorder()

    def preterminals(self):
        for n in self.postorder():
            if n.is_leaf():
                continue

            pre_terminal = True
            for c in n.child:
                if not c.is_leaf():
                    pre_terminal = False
                    break
            if pre_terminal:
                yield n

    def node(self, name):
        for n in self.postorder():
            if n.name == name:
                return n
        return None

    def get_leafs(self):
        if self.leaves is None:
            self.leaves = [node for node in self.postorder() if node.is_leaf()]

        return self.leaves

    def num_leafs(self):
        if self.leaves is None:
            self.leaves = [node for node in self.postorder() if node.is_leaf()]

        return len(self.leaves)

    def attach(self, seqs):
        leafs = {n.name: n for n in self.get_leafs()}
        self.seqs = {leafs[name]:seq for name,seq in seqs.items()}

    def align(self, tmpdir, min_blk_len, circular, mu, beta, extensive, edge_window, edge_extend, log_stats=False, verbose=False):
        self.root.set_level(0) # NOTE: for debug logging
        stats = {}
        # ---------------------------------------------
        # internal functions
        # debugging function that will check reconstructed sequence against known real one.
        def check(seqs, G, verbose=False):
            nerror = 0
            uncompressed_length = 0
            for n in self.get_leafs():
                if n.name not in G.seqs:
                    continue

                seq  = seqs[n]
                orig = str(seq[:]).upper()
                rec  = G.extract(n.name)
                uncompressed_length += len(orig)
                if orig != rec:
                    nerror += 1

                    with open("test.fa", "w+") as out:
                        out.write(f">original\n{orig}\n")
                        out.write(f">reconstructed\n{rec}")

                    for i in range(len(orig)//100):
                        if (orig[i*100:(i+1)*100] != rec[i*100:(i+1)*100]):
                            log("-----------------")
                            log(f"original: {i} {orig[i*100:(i+1)*100]}")
                            log(f"graph:    {i} {rec[i*100:(i+1)*100]}")
                            breakpoint("inconsistency")

                            diffs = [i for i in range(len(rec)) if rec[i] != orig[i]]
                            pos   = [0]
                            seq   = G.seqs[n.name]
                            for nn in seq.nodes:
                                pos.append(pos[-1] + len(G.blks[nn.blk.id].extract(n.name, nn.num)))
                            pos = pos[1:]

                            testseqs = []
                            for nn in G.seqs[n.name].nodes:
                                if nn.strand == Strand.Plus:
                                    testseqs.append("".join(G.blks[nn.blk.id].extract(n.name, nn.num)))
                                else:
                                    testseqs.append("".join(rev_cmpl(G.blks[nn.blk.id].extract(n.name, nn.num))))

            if nerror == 0:
                log("all sequences correctly reconstructed")
                tlen = np.sum([len(x) for x in G.blks.values()])
                log(f"--- total graph length: {tlen}")
                log(f"--- total input sequence: {uncompressed_length}")
                log(f"--- compression: {uncompressed_length/tlen:1.2f}")
            else:
                raise ValueError("bad sequence reconstruction")

        def merge(node1, node2):
            if node1 != node2:
                graph1, fapath1 = node1.graph, node1.fapath
                graph2, fapath2 = node2.graph, node2.fapath
                graph    = Graph.fuse(graph1, graph2)
                graph, _ = graph.union(fapath1, fapath2, f"{tmpdir}/{n.name}", min_blk_len, mu, beta, extensive, edge_window, edge_extend)
            else:
                graph = node1.graph

            for i in range(MAXSELFMAPS):
                log(f"----> merge round {i}")
                # check(self.seqs, graph)
                itr = f"{tmpdir}/{n.name}_iter_{i}"
                with open(f"{itr}.fa", 'w') as fd:
                    graph.write_fasta(fd)
                graph, contin = graph.union(itr, itr, f"{tmpdir}/{n.name}_iter_{i}", min_blk_len, mu, beta, extensive, edge_window, edge_extend)
                if not contin:
                    break

            check(self.seqs, graph)
            return graph

        # --------------------------------------------
        # body

        if self.num_leafs() == 1:
            return Graph()

        # fix trivial graphs onto the leafs
        for i, n in enumerate(self.get_leafs()):
            seq      = self.seqs[n]
            n.graph  = Graph.from_seq(n.name, str(seq).upper(), circular)
            n.fapath = f"{tmpdir}/{n.name}"
            with open(f"{n.fapath}.fa", 'w') as fd:
                n.graph.write_fasta(fd)

        for n in self.postorder():
            if n.is_leaf():
                continue
            n.fapath = f"{tmpdir}/{n.name}"
            log(f"fusing {n.child[0].name} with {n.child[1].name} @ {n.name}")
            n.graph = merge(*n.child)
            # delete references to children graphs for cleanup
            for c in n.child:
                if log_stats:
                    stats[c.name] = {
                        'length' : [b.length for b in c.graph.blks.values()],
                        'depth'  : [b.depth  for b in c.graph.blks.values()],
                    }
                c.graph = None

            with open(f"{n.fapath}.fa", 'w') as fd:
                n.graph.write_fasta(fd)

                log((f"--> compression ratio: "
                       f"{n.graph.compress_ratio()}"))
                log((f"--> number of blocks: "
                       f"{len(n.graph.blks)}"))
                log((f"--> number of members: "
                       f"{len(n.graph.seqs)}"))

    def collect(self):
        if not self.root.graph:
            return None
        self.root.graph = Graph.connected_components(self.root.graph)
        return self.root.graph

    def write_nwk(self, wtr):
        self.root.to_nwk(wtr)
        wtr.write(";")

    def write_json(self, wtr, no_seqs=False):
        data = {'tree' : self.root.to_json(),
                'seqs' : None if no_seqs else {k.name:str(v) for k,v in self.seqs.items()}}
        wtr.write(json.dumps(data))
