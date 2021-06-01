using NoteOnTropicalMIS, Test, LightGraphs

@testset "enumerating" begin
    @test NoteOnTropicalMIS.onehot(ConfigEnumerator{10}, 3) == ConfigEnumerator{10}(BitVector[[0, 0, 1, 0, 0, 0, 0, 0, 0, 0]])
    code = random_regular_eincode(10, 3)
    res1 = mis_count(code)[]
    res2 = mis_enumerate(code)[]
    @test res1.n == res2.n
    @test res1.c == length(res2.c)
end