using GraphTensorNetworks, Random
using CUDA
CUDA.allowscalar(false)

using DelimitedFiles

function mis_counting(n, seed; writefile, sc_target, usecuda, maximal)
    folder = joinpath(homedir(), ".julia/dev/TropicalMIS", "project", "data")
    fname = joinpath(folder, "mis_degeneracy_L$n.dat")
    mask = Matrix{Bool}(reshape(readdlm(fname)[seed+1,4:end], n, n))
    if graph == "diag"
        g = diagonal_coupled_graph(mask)
    elseif graph == "square"
        g = square_lattice_graph(mask)
    end
    if maximal
        gp = MaximalIndependence(g; optimizer=TreeSA(sc_target=sc_target, niters=5))
    else
        gp = Independence(g; optimizer=TreeSA(sc_target=sc_target, sc_weight=1.0, ntrials=10, βs=0.01:0.05:15.0, niters=50, rw_weight=0.2), simplifier=MergeGreedy())
    end
    println("Graph $graph of size $n, seed= $seed, usecuda = $usecuda")
    res = graph_polynomial(gp, Val(:finitefield), usecuda=usecuda)[]
    if maximal
        folderout = joinpath(folder, "$graph", "maximal_polynomial_L$(n)")
    else
        folderout = joinpath(folder, "$graph", "independence_polynomial_L$(n)")
    end
    if !isdir(folderout)
        mkdir(folderout)
    end
    ofname = joinpath(folderout, "$(seed).dat")
    writefile && writedlm(ofname, res.coeffs)
end

# patch
const DEVICE = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : -1

if DEVICE >= 0
    CUDA.device!(DEVICE)
end

for L = [23]
    println("computing L = $L")
    for i=0:999
        @time mis_counting(L, i; writefile=true, sc_target=26, usecuda=DEVICE>=0, maximal=false)
    end
end
