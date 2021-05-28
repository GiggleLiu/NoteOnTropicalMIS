using NoteOnTropicalMIS, Test, OMEinsum

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