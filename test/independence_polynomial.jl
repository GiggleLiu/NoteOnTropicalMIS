using NoteOnTropicalMIS, Test, OMEinsum
using Mods, Polynomials, TropicalNumbers

@testset "independence_polynomial" begin
    code = random_regular_eincode(10, 3)
    code = optimize_greedy(code, uniformsize(code, 2))
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