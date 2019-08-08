import numpy as np
from Bio import Seq, SeqIO, SeqRecord
from block import Block

plus_strand = 1
minus_strand = -1

class Graph(object):
    """docstring for Graph"""
    def __init__(self):
        super(Graph, self).__init__()
        self.blocks = {}
        self.sequences = {}

    @classmethod
    def from_sequence(cls, sequence_name, sequence):
        new_graph = cls()
        b = Block.from_sequence(sequence_name, sequence)
        new_graph.blocks = {b.name:b}
        new_graph.sequences = {sequence_name: [(b.name, plus_strand)]}
        return new_graph

    @classmethod
    def fuse(cls, graph1, graph2):
        new_graph = Graph()
        new_graph.blocks = {}
        new_graph.blocks.update(graph1.blocks)
        new_graph.blocks.update(graph2.blocks)
        new_graph.sequences = {s:list(b) for s,b in list(graph1.sequences.items())+list(graph2.sequences.items())}
        return new_graph

    def merge_hit(self,hit):
        orig_ref_block = self.blocks[hit['ref']['name']]
        orig_seq_block = self.blocks[hit['query']['name']]
        ref_block = orig_ref_block[hit['ref']['start']:hit['ref']['end']]
        seq_block = orig_seq_block[hit['query']['start']:hit['query']['end']]

        if hit["orientation"]==minus_strand:
            seq_block = seq_block.reverse_complement()

        aln = {"ref_seq":"".join(ref_block.consensus),
               "query_seq":"".join(seq_block.consensus),
               "cigar":hit["cigar"], "ref_cluster":ref_block.sequences,
               "query_cluster":seq_block.sequences, 'ref_start':0}

        merged_block = Block.from_cluster_alignment(aln)
        self.blocks[merged_block.name] = merged_block

        for tmp_block, subhit, strand in [(orig_ref_block, hit['ref'], plus_strand),
                                          (orig_seq_block, hit['query'], hit["orientation"])]:
            new_blocks = []
            if subhit['start']:
                left = tmp_block[0:subhit['start']]
                new_blocks.append((left.name, plus_strand))
                self.blocks[left.name] = left
            new_blocks.append((merged_block.name, strand))
            if subhit['end']<len(tmp_block):
                right = tmp_block[subhit['end']:]
                new_blocks.append((right.name, plus_strand))
                self.blocks[right.name] = right

            for s in tmp_block.sequences:
                new_block_sequence = []
                for bi,b in enumerate(self.sequences[s]):
                    if b[0]==tmp_block.name:
                        orig_strand = b[1]
                        if orig_strand==plus_strand:
                            tmp_new_blocks = [(x, orig_strand*y) for x,y in new_blocks]
                        else:
                            tmp_new_blocks = [(x, orig_strand*y) for x,y in new_blocks][::-1]
                        new_block_sequence = self.sequences[s][:bi] + tmp_new_blocks + self.sequences[s][bi+1:]
                        print(b[0], bi, len(new_block_sequence), len(self.sequences[s]))
                        break
                if not len(new_block_sequence):
                    import ipdb; ipdb.set_trace()
                self.sequences[s] = new_block_sequence

        print('before', self.blocks.keys())
        remaining_blocks = set()
        for s in self.sequences:
            remaining_blocks.update([s[0] for s in self.sequences[s]])
        self.blocks = {b:self.blocks[b] for b in remaining_blocks}
        print('after', self.blocks.keys())


    def extract(self, name):
        seq = ""
        for b,strand in self.sequences[name]:
            tmp_seq = self.blocks[b].extract(name)
            if strand==plus_strand:
                seq += tmp_seq
            else:
                seq += Seq.reverse_complement(tmp_seq)

        return seq


    def to_fasta(self, fname):
        SeqIO.write([SeqRecord.SeqRecord(seq=Seq.Seq("".join(c.consensus)), id=c.name, description='')
                    for c in self.blocks.values()], fname, format='fasta')


    def to_json(self, fname):
        J = {}
        J['Isolate_names'] = list(self.sequences.keys())
        J['Plasmids'] = [[x[0] for x in self.sequences[s]] for s in J['Isolate_names']]
        nodes = {}
        for b in self.blocks:
            nodes[b] = {"ID":b,
                        "Genomes":{"Consensus":''.join(self.blocks[b].consensus),
                        "Alignment":{J["Isolate_names"].index(s):self.blocks[b].extract(s) for s in self.blocks[b].sequences}},
                        "Out_Edges":[], "In_Edges":[]}

        edges = {}
        for pname, p in zip(range(len(J["Isolate_names"])), J['Plasmids']):
            for i in range(len(p)-1):
                e = (p[i],p[i+1])
                if e in edges:
                    edges[e]["Isolates"].append(pname)
                else:
                    edges[e] = {"Source":e[0], "Target":e[1], "Isolates":[pname]}
            e = (p[-1],p[0])
            if e in edges:
                edges[e]["Isolates"].append(pname)
            else:
                edges[e] = {"Source":e[0], "Target":e[1], "Isolates":[pname]}
        for e in edges:
            nodes[e[0]]["Out_Edges"].append(edges[e])
            nodes[e[1]]["In_Edges"].append(edges[e])
        J["Nodes"] = nodes
        import json
        with open(fname, 'w') as fh:
            json.dump(J, fh)


if __name__ == '__main__':
    g1 = Graph.from_sequence("S1", "xxxxxxACACACyyyy")
    g2 = Graph.from_sequence("S2", "zzzGTCTGTqqq")

    g = Graph.fuse(g1,g2)

    # this defines the homology hit format, meant to be parsed from paf
    hit = {'ref': {'name':g.sequences['S1'][0][0], 'start':6, 'end':12},
           'query': {'name':g.sequences['S2'][0][0], 'start':3, 'end':9},
           'cigar':'6M', 'orientation':minus_strand}

    g.merge_hit(hit)

    print(g.extract('S1'))
    print(g.extract('S2'))
