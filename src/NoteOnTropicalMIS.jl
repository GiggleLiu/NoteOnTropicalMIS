module NoteOnTropicalMIS

using GenericTensorNetworks

project_relative_path(xs...) = normpath(joinpath(dirname(dirname(pathof(@__MODULE__))), xs...))

end
