using GraphTensorNetworks, Random

using DelimitedFiles

function mis_counting(n, seed; writefile, sc_target, usecuda, maximal)
    folder = joinpath(homedir(), ".julia/dev/TropicalMIS", "project", "data")
    fname = joinpath(folder, "mis_degeneracy_L$n.dat")
    mask = Matrix{Bool}(reshape(readdlm(fname)[seed+1,4:end], n, n))
    g = diagonal_coupled_graph(mask)
    if maximal
        gp = MaximalIndependence(g; sc_target=sc_target, optmethod=:tree, niters=5)
    else
        #gp = Independence(g; optmethod=:kahypar, sc_target=sc_target)
        gp = Independence(g; optmethod=:tree, sc_target=sc_target, sc_weight=1.0, ntrials=10, βs=0.01:0.05:15.0, niters=30, rw_weight=0.2)
    end
    println("Graph $seed")
    res = graph_polynomial(gp, Val(:finitefield), usecuda=usecuda)[]
    if maximal
        folderout = joinpath(folder, "maximal_polynomial_L$(n)")
    else
        folderout = joinpath(folder, "independence_polynomial_L$(n)")
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

for i=0:49
    @time mis_counting(24, i; writefile=true, sc_target=0, usecuda=DEVICE>=0, maximal=false)
end
