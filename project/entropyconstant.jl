using GraphTensorNetworks, Random
using CUDA
CUDA.allowscalar(false)
using MKL

using DelimitedFiles
using LinearAlgebra

BLAS.set_num_threads(1)

function sequencing(n; writefile, sc_target, usecuda, nslices=1)
    g = square_lattice_graph(trues(n, n))
    gp = Independence(g; optimizer=TreeSA(sc_target=sc_target, sc_weight=1.0, nslices=nslices;
        ntrials=7, βs=0.01:0.05:22.0, niters=20, rw_weight=2.0), simplifier=MergeGreedy())
    println("Graph size $n, usecuda = $usecuda")
    @show timespace_complexity(gp)
    @time res = Array(GraphTensorNetworks.contractx(gp, 1.0; usecuda=usecuda))
    @show res
    ofname = joinpath(@__DIR__, "data", "entropyconstant-$n.dat")
    writefile && writedlm(ofname, res)
end

# patch
const DEVICE = parse(Int, ARGS[1])

if DEVICE >= 0
    CUDA.device!(DEVICE)
end

Random.seed!(2)
for L=34:37
    println("computing L = $L")
    @time sequencing(L; writefile=true, sc_target=23, usecuda=DEVICE>=0, nslices=L-23)
end
