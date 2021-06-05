using Test, OMEinsum, NoteOnTropicalMIS, TropicalNumbers, Random
using NoteOnTropicalMIS: cached_einsum, generate_masktree, masked_einsum, CacheTree

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
    for seed in 1:100
        Random.seed!(seed)
        xs = map(x->TropicalF64.(x), [rand(1:5,2,2), rand(1:5,2), rand(1:5,2,2), rand(1:5,2,2), rand(1:5,2,2)])
        code = ein"((ij,j),jk, kl), ii->kli"
        y1 = code(xs...)
        y2 = bounding_contract(code, xs, BitArray(ones(Bool,2,2,2)), xs)
        @test y1 ≈ y2
    end
    rawcode = random_regular_eincode(10, 3)
    optcode = OMEinsum.optimize_greedy(rawcode, uniformsize(rawcode, 2))
    xs = NoteOnTropicalMIS.generate_xs!(NoteOnTropicalMIS._auto_mistensor, TropicalF64(1.0), rawcode, Vector{Any}(undef, NoteOnTropicalMIS.ninput(rawcode)))
    y1 = rawcode(xs...)
    y2 = bounding_contract(rawcode, xs, BitArray(fill(true)), xs)
    @test y1 ≈ y2
    y1 = optcode(xs...)
    y2 = bounding_contract(optcode, xs, BitArray(fill(true)), xs)
    @test y1 ≈ y2
end

