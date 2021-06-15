module NoteOnTropicalMIS

using OMEinsumContractionOrders: OMEinsum
using Core: Argument
using TropicalGEMM, TropicalNumbers
using OMEinsum

project_relative_path(xs...) = normpath(joinpath(dirname(dirname(pathof(@__MODULE__))), xs...))

# patch
if isdefined(OMEinsum, :flatten)
    using OMEinsum: flatten
else
    flatten(code) = Iterators.flatten(code::Union{NestedEinsum,EinCode})
end

using TensorOperations, LinearAlgebra
function LinearAlgebra.permutedims!(C::Array{T,N}, A::StridedArray{T,N}, perm) where {T,N}
    if isbitstype(T)
        TensorOperations.tensorcopy!(A, ntuple(identity,N), C, perm)
    else
        invoke(permutedims!, Tuple{Any,AbstractArray,Any}, C, A, perm)
    end
end

include("arithematics.jl")
include("widgets.jl")
include("independence_polynomial.jl")
include("configurations.jl")
include("graphs.jl")
include("bounding.jl")
include("viz.jl")
include("demos.jl")

using Requires
function __init__()
    @require CUDA="052768ef-5323-5732-b1bb-66c8b64840ba" include("cuda.jl")
end

end
