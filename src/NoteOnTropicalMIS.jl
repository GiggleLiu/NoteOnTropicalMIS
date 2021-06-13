module NoteOnTropicalMIS

using TropicalGEMM, TropicalNumbers
using OMEinsum

project_relative_path(xs...) = normpath(joinpath(dirname(dirname(pathof(@__MODULE__))), xs...))

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
