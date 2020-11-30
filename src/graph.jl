module PanGraph

include("block.jl")
include("path.jl")
include("align.jl")

export Graph, write

# ------------------------------------------------------------------------
# graph data structure

struct Graph
    blocks::Dict{String,Block}
    sequence::Dict{String,Path}
    # TODO: add edge data structure
end

# ------------------------------------------------------------------------
# serialization

function write_fasta(io, G::Graph; numcols=80)
    NL = "\n"
    for block in values(G.blocks)
        write(io, "f>{block.uuid}{NL}")
        write(io, join([block.sequence[1+numcols*i:numcols*(i+1)] for i in 1:ceil(length(block),numcols)], NL))
    end
end

function write(io, G::Graph; fmt=:fasta)
    @match fmt begin
        :fasta || :fa => write_fasta(io, G)
        _ => error(f"{format} not a recognized output format")
    end
end

end
