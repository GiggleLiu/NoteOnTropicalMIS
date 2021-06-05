module NoteOnTropicalMIS

using TropicalGEMM, TropicalNumbers
using OMEinsum

include("arithematics.jl")
include("tensors.jl")
include("graphs.jl")
include("independence_polynomial.jl")
include("configurations.jl")
include("random_regular_graph.jl")
include("ad.jl")
include("viz.jl")

end
