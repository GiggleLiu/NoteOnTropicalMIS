using NoteOnTropicalMIS, Test, OMEinsum, OMEinsumContractionOrders
using Mods, Polynomials, TropicalNumbers
using LightGraphs, Random

@testset "bond and vertex tensor" begin
    @test misb(TropicalF64) == [TropicalF64(0) TropicalF64(0); TropicalF64(0) TropicalF64(-Inf)]
    @test misv(TropicalF64, 2.0) == [TropicalF64(0), TropicalF64(2.0)]
end

@testset "graph generator" begin
    g = diagonal_coupled_graph(trues(3, 3))
    @test ne(g) == 20
    g = diagonal_coupled_graph((x = trues(3, 3); x[2,2]=0; x))
    @test ne(g) == 12
    g = diagonal_coupled_eincode(trues(3, 3))
    @test length(NoteOnTropicalMIS.symbols(g)) == 9
end

@testset "independence_polynomial" begin
    Random.seed!(2)
    code = random_regular_eincode(10, 3)
    code = optimize_kahypar(code, uniformsize(code, 2), sc_target=4, max_group_size=5)
    p1 = independence_polynomial(Val(:fitting), code)
    p2 = independence_polynomial(Val(:polynomial), code)
    p3 = independence_polynomial(Val(:fft), code)
    p4 = independence_polynomial(Val(:finitefield), code)
    @test p1 ≈ p2
    @test p1 ≈ p3
    @test p1 ≈ p4
end

@testset "arithematics" begin
    for (a, b, c) in [
                    (TropicalF64(2), TropicalF64(8), TropicalF64(9)),
                    (CountingTropicalF64(2, 8), CountingTropicalF64(8, 9), CountingTropicalF64(9, 2)),
                    (Mod{17}(2), Mod{17}(8), Mod{17}(9)),
                    (Polynomial([0,1,2,3.0]), Polynomial([3,2.0]), Polynomial([1,7.0])),
                    (TropicalF64(5), TropicalF64(3), TropicalF64(-9)),
                    (CountingTropicalF64(5, 3), CountingTropicalF64(3, 9), CountingTropicalF64(-3, 2)),
                    (ConfigTropical{Float64,10,1}(5.0, BitVector(rand(Bool, 10))), ConfigTropical{Float64,10,1}(3.0, BitVector(rand(Bool, 10))), ConfigTropical{Float64,10,1}(-3.0, BitVector(rand(Bool, 10)))),
                    (CountingTropical(5.0, ConfigEnumerator([StaticBitVector(rand(Bool, 10)) for j=1:3])), CountingTropical(3.0, ConfigEnumerator([StaticBitVector(rand(Bool, 10)) for j=1:4])), CountingTropical(-3.0, ConfigEnumerator([StaticBitVector(rand(Bool, 10)) for j=1:5]))),
                    ]
        @test is_commutative_semiring(a, b, c)
    end
end