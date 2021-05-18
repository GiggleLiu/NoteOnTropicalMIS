using Test, NoteOnTropicalMIS
using TropicalNumbers

@testset "bond and vertex tensor" begin
    @test misb(TropicalF64) == [TropicalF64(0) TropicalF64(0); TropicalF64(0) TropicalF64(-Inf)]
    @test misv(TropicalF64, 1, 2.0) == [TropicalF64(0), TropicalF64(2.0)]
end

edges = find_edges(getindex.(nodes, 2), sqrt(0.05))
@test length(edges) == 32
vizconfig(nodes, edges)
res2 = contract_a(TropicalF64, ones(16))

