using NoteOnTropicalMIS, Polynomials, Random, OMEinsum

#Random.seed!(2)
#g = diagonal_coupled_graph(rand(13,13) .<= 0.8)
#@time maximal_polynomial(Val(:fft), g; sc_target=24, imbalances=0.0:0.00151:1.0, max_group_size=50)
#@time maximal_polynomial(Val(:polynomial), g; sc_target=19, imbalances=0.0:0.00151:1.0, max_group_size=50)

using DelimitedFiles
using NoteOnTropicalMIS.OMEinsum

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

function mis_maximal_counting(n, ρ, ; writefile=false, sc_target=24, usecuda=false, folder=".")
    mask = rand(n, n) .< ρ
    g = diagonal_coupled_graph(mask)
    println("Graph n=$n, ρ=$ρ")
    #res = maximal_polynomial(Val(:fft), g; sc_target=sc_target, imbalances=0.0:0.002:1.0, max_group_size=50, usecuda=usecuda, r=3.0)
    res = maximal_polynomial(Val(:finitefield), g; sc_target=sc_target, imbalances=0.0:0.002:1.0, max_group_size=50, usecuda=usecuda)
    ofname = joinpath(folder, "maximal_counting_$(n)x$(n)_rho$ρ.dat")
    writefile && writedlm(ofname, res.coeffs)
end

const DEVICE = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : -1

if DEVICE >= 0
    using CUDA
    CUDA.device!(DEVICE)
end

#mis_maximal_counting(15, 0.8, writefile=true, usecuda=DEVICE>=0, sc_target=26)
