using Viznet
using AbstractTrees
using OMEinsum: NestedEinsum

function AbstractTrees.children(ne::NestedEinsum)
    d = Dict()
    for (k,item) in enumerate(ne.args)
        d[k] = item isa Integer ? join(OMEinsum.getixs(ne.eins)[k]) : item
    end
    d
end

function AbstractTrees.printnode(io::IO, x::String)
    print(io, x)
end
function AbstractTrees.printnode(io::IO, x::NestedEinsum)
    print(io, x.eins)
end

function Base.show(io::IO, e::EinCode{ixs, iy}) where {ixs, iy}
    s = join([_join(ix) for ix in ixs], ", ") * " -> " * _join(iy)
    print(io, s)
end
function Base.show(io::IO, e::NestedEinsum)
    print_tree(io, e)
end
Base.show(io::IO, ::MIME"text/plain", e::NestedEinsum) = show(io, e)
Base.show(io::IO, ::MIME"text/plain", e::EinCode) = show(io, e)
_join(ix::NTuple{0}) = ""
_join(ix::NTuple{N,Char}) where N = join(ix, "")
_join(ix::NTuple{N,Int}) where N = join(ix, "∘")