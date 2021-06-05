using TupleTools

export bounding_contract

Base.isnan(x::Tropical) = isnan(x.n)
function backward_tropical(@nospecialize(ixs), @nospecialize(xs), @nospecialize(iy), @nospecialize(y), @nospecialize(ymask), i, size_dict)
    nixs = TupleTools.insertat(ixs, i, (iy,))
    nxs  = TupleTools.insertat( xs, i, (OMEinsum.asarray(inv.(y) .* ymask),))
    niy = ixs[i]
    return inv.(einsum(EinCode(nixs, niy), nxs, size_dict)) .== xs[i]
end

struct CacheTree{T}
    content::AbstractArray{T}
    siblings::Vector{CacheTree{T}}
end
function cached_einsum(code::Int, @nospecialize(xs), size_dict)
    y = xs[code]
    CacheTree(y, CacheTree{eltype(y)}[])
end
function cached_einsum(code::NestedEinsum, @nospecialize(xs), size_dict)
    caches = [cached_einsum(arg, xs, size_dict) for arg in code.args]
    y = einsum(code.eins, (getfield.(caches, :content)...,), size_dict)
    CacheTree(y, caches)
end

function generate_masktree(code::Int, cache, mask, size_dict)
    CacheTree(mask, CacheTree{Bool}[])
end
function generate_masktree(code::NestedEinsum, cache, mask, size_dict)
    submasks = map(1:length(code.args)) do i
        backward_tropical(OMEinsum.getixs(code.eins), (getfield.(cache.siblings, :content)...,), OMEinsum.getiy(code.eins), cache.content, mask, i, size_dict)
    end
    return CacheTree(mask, generate_masktree.(code.args, cache.siblings, submasks, Ref(size_dict)))
end

function masked_einsum(code::Int, @nospecialize(xs), masks, size_dict)
    y = copy(xs[code])
    y[OMEinsum.asarray(.!masks.content)] .= Ref(zero(eltype(y))); y
end
function masked_einsum(code::NestedEinsum, @nospecialize(xs), masks, size_dict)
    xs = [masked_einsum(arg, xs, mask, size_dict) for (arg, mask) in zip(code.args, masks.siblings)]
    y = einsum(code.eins, (xs...,), size_dict)
    y[OMEinsum.asarray(.!masks.content)] .= Ref(zero(eltype(y))); y
end

function bounding_contract(@nospecialize(code::EinCode), @nospecialize(xsa), ymask, @nospecialize(xsb); size_info=nothing)
    bounding_contract(NestedEinsum((1:length(xsa)), code), xsa, ymask, xsb; size_info=size_info)
end
function bounding_contract(code::NestedEinsum, @nospecialize(xsa), ymask, @nospecialize(xsb); size_info=nothing)
    size_dict = OMEinsum.get_size_dict(OMEinsum.getixs(Iterators.flatten(code)), xsa, size_info)
    # compute intermediate tensors
    c = cached_einsum(code, xsa, size_dict)
    # compute masks from cached tensors
    mt = generate_masktree(code, c, ymask, size_dict)
    # compute results with masks
    masked_einsum(code, xsb, mt, size_dict)
end