module NoteOnTropicalMIS

using TropicalGEMM, TropicalNumbers
using OMEinsum

include("tensors.jl")
include("graphs.jl")
include("independence_polynomial.jl")
include("random_regular_graph.jl")
include("viz.jl")

end
