using GraphTensorNetworks, Random
using CUDA
CUDA.allowscalar(false)

using DelimitedFiles

function sequencing(n; writefile, sc_target, usecuda, nslices=1)
    g = square_lattice_graph(trues(n, n))
    gp = Independence(g; optimizer=TreeSA(sc_target=sc_target, sc_weight=1.0, nslices=nslices;
        ntrials=5, βs=0.01:0.05:22.0, niters=20, rw_weight=0.2), simplifier=MergeGreedy())
    println("Graph size $n, usecuda = $usecuda")
    @show timespace_complexity(gp)
    return
    res = GraphTensorNetworks.big_integer_solve(Int32, 100) do T
        GraphTensorNetworks.contractx(gp, one(T); usecuda=usecuda)
    end
    #res = GraphTensorNetworks.contractx(gp, 1.0; usecuda=usecuda)[]
    @show res
    ofname = joinpath(@__DIR__, "data", "$n.dat")
    writefile && writedlm(ofname, res)
end

# patch
const DEVICE = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : -1

if DEVICE >= 0
    CUDA.device!(DEVICE)
end

for L = [11]
    println("computing L = $L")
    @time sequencing(L; writefile=false, sc_target=10, usecuda=DEVICE>=0)
end
