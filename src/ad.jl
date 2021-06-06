using TupleTools

export bounding_contract

Base.isnan(x::Tropical) = isnan(x.n)
function backward_tropical(mode, @nospecialize(ixs), @nospecialize(xs), @nospecialize(iy), @nospecialize(y), @nospecialize(ymask), i, size_dict)
    nixs = TupleTools.insertat(ixs, i, (iy,))
    nxs  = TupleTools.insertat( xs, i, (OMEinsum.asarray(inv.(y) .* ymask),))
    niy = ixs[i]
    A = inv.(einsum(EinCode(nixs, niy), nxs, size_dict))
    if mode == :all
        return A .== xs[i]
    elseif mode == :single  # wrong, need `B` matching `A`.
        mask = falses(size(A)...)
        found = false
        for j=1:length(A)
            if xs[i][j] == A[j] && !found
                #@show j, A, xs[i]
                mask[j] = true
                found = true
            else
                xs[i][j] = zero(eltype(xs[i]))
            end
        end
        return mask
    else
        error("unkown mode: $mod")
    end
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

function generate_masktree(code::Int, cache, mask, size_dict, mode=:all)
    CacheTree(mask, CacheTree{Bool}[])
end
function generate_masktree(code::NestedEinsum, cache, mask, size_dict, mode=:all)
    submasks = map(1:length(code.args)) do i
        backward_tropical(mode, OMEinsum.getixs(code.eins), (getfield.(cache.siblings, :content)...,), OMEinsum.getiy(code.eins), cache.content, mask, i, size_dict)
    end
    return CacheTree(mask, generate_masktree.(code.args, cache.siblings, submasks, Ref(size_dict), mode))
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
    mt = generate_masktree(code, c, ymask, size_dict, :all)
    # compute results with masks
    masked_einsum(code, xsb, mt, size_dict)
end

function mis_config_ad(@nospecialize(code::EinCode), @nospecialize(xsa), ymask; size_info=nothing)
    mis_config_ad(NestedEinsum((1:length(xsa)), code), xsa, ymask; size_info=size_info)
end

function mis_config_ad(code::NestedEinsum, @nospecialize(xsa), ymask; size_info=nothing)
    size_dict = OMEinsum.get_size_dict(OMEinsum.getixs(Iterators.flatten(code)), xsa, size_info)
    # compute intermediate tensors
    c = cached_einsum(code, xsa, size_dict)
    # compute masks from cached tensors
    mt = generate_masktree(code, c, ymask, size_dict, :single)
    c.content, read_config!(code, mt, Dict())
end

function read_config!(code::NestedEinsum, mt, out)
    for (arg, ix, sibling) in zip(code.args, OMEinsum.getixs(code.eins), mt.siblings)
        if arg isa Int
            assign = sibling.content
            if length(ix) == 1
                if !assign[1] && assign[2]
                    out[ix[1]] = 1
                elseif !assign[2] && assign[1]
                    out[ix[1]] = 0
                else
                    error("invalid assign $(assign)")
                end
            end
        else  # nested
            read_config!(arg, sibling, out)
        end
    end
    return out
end