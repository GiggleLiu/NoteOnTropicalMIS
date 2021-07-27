using NoteOnTropicalMIS, Polynomials, Random, OMEinsum

using DelimitedFiles
using NoteOnTropicalMIS.OMEinsum

function mis_maximal_counting(n, seed; writefile, rs=[1.0], sc_target=24, method=:finitefield, usecuda=false)
    folder = joinpath("/home/leo/.julia/dev/TropicalMIS", "project", "data")
    fname = joinpath(folder, "mis_degeneracy_L$n.dat")
    mask = Matrix{Bool}(reshape(readdlm(fname)[seed+1,4:end], n, n))
    g = diagonal_coupled_graph(mask)
    println("Graph $seed")
    if method == :fft
        for r in rs
            res = maximal_polynomial(Val(:fft), g; sc_target=sc_target, imbalances=0.0:0.002:1.0, max_group_size=50, r=r, usecuda=usecuda)
            ofname = joinpath(folder, "maximal_counting_$(n)x$(n)_$(seed)_r$(r).dat")
            writefile && writedlm(ofname, res.coeffs)
        end
    elseif method == :finitefield
        res = maximal_polynomial(Val(:finitefield), g; sc_target=sc_target, imbalances=0.0:0.002:1.0, max_group_size=50, usecuda=usecuda)
        ofname = joinpath(folder, "maximal_counting_$(n)x$(n)_$(seed).dat")
        writefile && writedlm(ofname, res.coeffs)
    end
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

for seed in [0, 2, 4, 9, 50]
    @time mis_maximal_counting(10, seed; writefile=true, sc_target=20, usecuda=DEVICE>=0)
end

for seed in [26, 40, 339, 542, 339, 8, 13, 25, 89, 108, 138]
    @time mis_maximal_counting(15, seed; writefile=true, sc_target=26, usecuda=DEVICE>=0)
end