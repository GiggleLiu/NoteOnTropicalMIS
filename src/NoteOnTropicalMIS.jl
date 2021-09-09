module NoteOnTropicalMIS

using GraphTensorNetworks

project_relative_path(xs...) = normpath(joinpath(dirname(dirname(pathof(@__MODULE__))), xs...))

include("widgets.jl")
include("demos.jl")

end
