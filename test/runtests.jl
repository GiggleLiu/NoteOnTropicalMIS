using NoteOnTropicalMIS
using Test

@testset "tensors" begin
    include("tensors.jl")
end

@testset "independence polynomial" begin
    include("independence_polynomial.jl")
end

@testset "configurations" begin
    include("configurations.jl")
end