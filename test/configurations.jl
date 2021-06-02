using NoteOnTropicalMIS, Test, LightGraphs

@testset "ConfigTropical" begin
    x = one(ConfigTropical{Float64, 5})
    @test x.n == 0
    @test x.config == falses(5)
    x = zero(ConfigTropical{Float64, 5})
    @test x.n == -Inf
    @test x.config == trues(5)
end

@testset "enumerating" begin
    @test NoteOnTropicalMIS.onehot(ConfigEnumerator{10}, 3) == ConfigEnumerator{10}(BitVector[[0, 0, 1, 0, 0, 0, 0, 0, 0, 0]])
    code = random_regular_eincode(10, 3)
    res1 = mis_count(code)[]
    res2 = mis_config(code; all=true)[]
    res3 = mis_config(code; all=false)[]
    @test res1.n == res2.n == res3.n
    @test res1.c == length(res2.c)
    @test res3.config ∈ res2.c.data
end