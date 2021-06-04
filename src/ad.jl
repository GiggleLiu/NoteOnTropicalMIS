using ChainRulesCore
using TropicalNumbers
using ForwardDiff
using OMEinsum: NestedEinsum
using OMEinsum
using TupleTools

a = Float64.(rand(1:5, 10, 10))
A = TropicalF64.(a)
B = TropicalF64.(rand(1:5, 10, 10))
C = A * B
gC = ones(TropicalF64, 10, 10)

function loss(a)
    A = Tropical.(a)
    C = A * B
    prod(C).n
end
ForwardDiff.gradient(loss, a)

function backward_tropicalmatmul(A, B, C, CMASK)
    (inv.(inv.(C .* CMASK) * transpose(B)) .== A)
end

function backward_tropical(ixs, xs, iy, y, ymask, i; size_info)
    nixs = TupleTools.insertat(ixs, i, (iy,))
    nxs  = TupleTools.insertat( xs, i, (inv.(y .* ymask),))
    niy = ixs[i]
    return inv.(EinCode(nixs, niy)(nxs...; size_info=size_info)) .== xs[i]
end

struct CacheTree{T}
    content::AbstractArray{T}
    siblings::Vector{CacheTree{T}}
end
function cached_einsum(code::Int, xs; size_info=nothing)
    y = xs[code]
    CacheTree(y, CacheTree{eltype(y)}[])
end
function cached_einsum(code::NestedEinsum, xs; size_info=nothing)
    caches = [cached_einsum(arg, xs; size_info=size_info) for arg in code.args]
    y = code.eins(getfield.(caches, :content)...; size_info=size_info)
    CacheTree(y, caches)
end

function generate_masktree(code::Int, cache, mask; size_info=nothing)
    CacheTree(mask, CacheTree{Bool}[])
end
function generate_masktree(code::NestedEinsum, cache, mask; size_info=nothing)
    submasks = map(1:length(code.args)) do i
        backward_tropical(OMEinsum.getixs(code.eins), (getfield.(cache.siblings, :content)...,), OMEinsum.getiy(code.eins), cache.content, mask, i; size_info=size_info)
    end
    return CacheTree(mask, generate_masktree.(code.args, cache.siblings, submasks; size_info=size_info))
end

function masked_einsum(code::Int, xs, masks; size_info=nothing)
    y = xs[code]
    y[.!masks.content] .= Ref(zero(eltype(y))); y
end
function masked_einsum(code::NestedEinsum, xs, masks; size_info=nothing)
    xs = [masked_einsum(arg, xs, mask; size_info=size_info) for (arg, mask) in zip(code.args, masks.siblings)]
    y = code.eins(xs...; size_info=size_info)
    y[.!masks.content] .= Ref(zero(eltype(y))); y
end

function bounding_contract(code, xsa, ymask, xsb; size_info=nothing)
    # compute intermediate tensors
    c = cached_einsum(code, xsa; size_info=size_info)
    # compute masks from cached tensors
    mt = generate_masktree(code, c, ymask; size_info=size_info)
    # compute results with masks
    masked_einsum(code, xsb, mt; size_info=size_info)
end

using Test
@testset "cached einsum" begin
    xs = map(x->Tropical.(x), [randn(2,2), randn(2), randn(2,2), randn(2,2), randn(2,2)])
    code = ein"((ij,j),jk, kl), ii->kli"
    c = cached_einsum(code, xs)
    @test c.content == code(xs...)
    mt = generate_masktree(code, c, rand(Bool,2,2,2))
    @test mt isa CacheTree{Bool}
    y = masked_einsum(code, xs, mt)
    @test y isa AbstractArray
end

@testset "bounding contract" begin
    xs = map(x->Tropical.(x), [rand(1:5,2,2), rand(1:5,2), rand(1:5,2,2), rand(1:5,2,2), rand(1:5,2,2)])
    code = ein"((ij,j),jk, kl), ii->kli"
    y1 = code(xs...)
    y2 = bounding_contract(code, xs, BitArray(ones(Bool,2,2,2)), xs)
    @test y1 ≈ y2
end