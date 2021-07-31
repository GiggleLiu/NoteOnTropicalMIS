using NoteOnTropicalMIS, Polynomials, Random, OMEinsum

using DelimitedFiles
using NoteOnTropicalMIS.OMEinsum

function mis_maximal_counting(n, seed; writefile, sc_target=24, usecuda=false, imbalances=0.0:0.002:1.0)
    folder = joinpath(homedir(), ".julia/dev/TropicalMIS", "project", "data")
    fname = joinpath(folder, "mis_degeneracy_L$n.dat")
    mask = Matrix{Bool}(reshape(readdlm(fname)[seed+1,4:end], n, n))
    g = diagonal_coupled_graph(mask)
    println("Graph $seed")
    res = maximal_polynomial(Val(:finitefield), g; sc_target=sc_target, imbalances=imbalances, max_group_size=50, usecuda=usecuda)
    folderout = joinpath(folder, "maximal_polynomial_L$(n)")
    if !isdir(folderout)
        mkdir(folderout)
    end
    ofname = joinpath(folderout, "$(seed).dat")
    writefile && writedlm(ofname, res.coeffs)
end

# patch
using GPUArrays, CUDA, LinearAlgebra
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

const DEVICE = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : -1

if DEVICE >= 0
    CUDA.device!(DEVICE)
end

#for seed in [0, 2, 4, 9, 50]
#    @time mis_maximal_counting(10, seed; writefile=true, sc_target=20, usecuda=DEVICE>=0)
#end

for n = [11]
    for seed in 488:699
        try
            @time mis_maximal_counting(n, seed; writefile=true, sc_target=20, usecuda=DEVICE>=0, imbalances=0.03:0.001:1.0)
        catch e
        end
    end
end
