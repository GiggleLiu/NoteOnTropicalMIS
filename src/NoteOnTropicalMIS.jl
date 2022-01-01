module NoteOnTropicalMIS

using GraphTensorNetworks

project_relative_path(xs...) = normpath(joinpath(dirname(dirname(pathof(@__MODULE__))), xs...))

end
