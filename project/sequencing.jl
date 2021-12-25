using GraphTensorNetworks, Random
using CUDA
CUDA.allowscalar(false)

using DelimitedFiles

function sequencing(n; writefile, sc_target, usecuda, nslices=1)
    g = square_lattice_graph(trues(n, n))
    gp = Independence(g; optimizer=TreeSA(sc_target=sc_target, sc_weight=1.0, nslices=nslices;
        ntrials=7, βs=0.01:0.05:22.0, niters=20, rw_weight=2.0), simplifier=MergeGreedy())
    println("Graph size $n, usecuda = $usecuda")
    @show timespace_complexity(gp)
    res = GraphTensorNetworks.big_integer_solve(Int32, 100) do T
        @time Array(GraphTensorNetworks.contractx(gp, one(T); usecuda=usecuda))
    end
    #res = GraphTensorNetworks.contractx(gp, 1.0; usecuda=usecuda)[]
    @show res
    ofname = joinpath(@__DIR__, "data", "$n.dat")
    writefile && writedlm(ofname, res)
end

# patch
const DEVICE = parse(Int, ARGS[1])

if DEVICE >= 0
    CUDA.device!(DEVICE)
end

for L=31:40
    println("computing L = $L")
    @time sequencing(L; writefile=true, sc_target=27, usecuda=DEVICE>=0, nslices=L-27)
end
