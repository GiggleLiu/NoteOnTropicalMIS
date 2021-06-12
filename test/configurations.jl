using NoteOnTropicalMIS, Test, LightGraphs
using OMEinsum

@testset "ConfigTropical" begin
    x = one(ConfigTropical{Float64, 5, 1})
    @test x.n == 0
    @test x.config == falses(5)
    x = zero(ConfigTropical{Float64, 5, 1})
    @test x.n == -Inf
    @test x.config == trues(5)
end

@testset "enumerating" begin
    rawcode = random_regular_eincode(10, 3)
    optcode = OMEinsum.optimize_greedy(rawcode, uniformsize(rawcode, 2))
    for code in [rawcode, optcode]
        res0 = mis_size(code)
        res1 = mis_count(code)
        res2 = mis_config(code; all=true, usemask=true)[]
        res3 = mis_config(code; all=false, usemask=false)[]
        res4 = mis_config(code; all=true, usemask=false)[]
        @test res0 == res2.n == res3.n == res4.n
        @test res1 == length(res2.c) == length(res4.c)
        @test res3.config ∈ res2.c.data
        @test res3.config ∈ res4.c.data
        res5 = mis_config(code; all=false, usemask=true)[]
        @test res5.n == res0
        @test res5.config ∈ res2.c.data
    end
end