using Test, NoteOnTropicalMIS
using TropicalNumbers

@testset "bond and vertex tensor" begin
    @test misb(TropicalF64) == [TropicalF64(0) TropicalF64(0); TropicalF64(0) TropicalF64(-Inf)]
    @test misv(TropicalF64, 1, 2.0) == [TropicalF64(0), TropicalF64(2.0)]
end