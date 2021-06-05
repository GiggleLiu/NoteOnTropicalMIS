using Test, NoteOnTropicalMIS
using TropicalNumbers
using NoteOnTropicalMIS: unitdisk_graph, diagonal_coupled_graph

@testset "bond and vertex tensor" begin
    @test misb(TropicalF64) == [TropicalF64(0) TropicalF64(0); TropicalF64(0) TropicalF64(-Inf)]
    @test misv(TropicalF64, 1, 2.0) == [TropicalF64(0), TropicalF64(2.0)]
end

@testset "graph generator" begin
    g = diagonal_coupled_graph(trues(3, 3))
    @test ne(g) == 20
    g = diagonal_coupled_graph((x = trues(3, 3); x[2,2]=0; x))
    @test ne(g) == 12
    g = diagonal_coupled_eincode(trues(3, 3))
    @test NoteOnTropicalMIS.ninput(g) == 29
end