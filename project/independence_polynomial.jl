using GraphTensorNetworks, Random

using DelimitedFiles

function mis_maximal_counting(n, seed; writefile, sc_target=24, usecuda=false)
    folder = joinpath(homedir(), ".julia/dev/TropicalMIS", "project", "data")
    fname = joinpath(folder, "mis_degeneracy_L$n.dat")
    mask = Matrix{Bool}(reshape(readdlm(fname)[seed+1,4:end], n, n))
    g = diagonal_coupled_graph(mask)
    gp = MaximalIndependence(g; sc_target=sc_target, optmethod=:tree, niters=5)
    println("Graph $seed")
    res = graph_polynomial(gp, Val(:finitefield), usecuda=usecuda)[]
    folderout = joinpath(folder, "maximal_polynomial_L$(n)")
    if !isdir(folderout)
        mkdir(folderout)
    end
    ofname = joinpath(folderout, "$(seed).dat")
    writefile && writedlm(ofname, res.coeffs)
end

# patch
using GPUArrays, CUDA, LinearAlgebra
@static if !(VERSION > v"1.6")
function LinearAlgebra.permutedims!(dest::AbstractGPUArray, src::AbstractGPUArray,
                                    perm::NTuple)
    Base.checkdims_perm(dest, src, perm)
    function permutedims_kernel(ctx, dest, src, perm)
        I = @cartesianidx src
        @inbounds begin
            J = CartesianIndex(map(i->I[i], perm))
            dest[J] = src[I]
        end
        return
    end
    gpu_call(permutedims_kernel, dest, src, perm)
    return dest
end
end

const DEVICE = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : -1

if DEVICE >= 0
    CUDA.device!(DEVICE)
end

for i=0:999
    @time mis_maximal_counting(12, i; writefile=true, sc_target=0, usecuda=DEVICE>=0)
end
