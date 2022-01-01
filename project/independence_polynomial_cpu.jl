using GraphTensorNetworks, Random
using Distributed, DelimitedFiles, GraphTensorNetworks.OMEinsumContractionOrders, GraphTensorNetworks.OMEinsum

@everywhere begin
using GraphTensorNetworks, Random, GraphTensorNetworks.OMEinsumContractionOrders, GraphTensorNetworks.OMEinsum
using DelimitedFiles

function do_work(f, jobs, results) # define work function everywhere
    while true
        job = take!(jobs)
        @info "running $job on device $(Distributed.myid())"
        res = f(job)
        put!(results, res)
    end
end

function mis_counting(graph, n, seed, mask; writefile, sc_target, usecuda, maximal, folderout)
    if graph == "diag"
        g = diagonal_coupled_graph(mask)
    elseif graph == "square"
        g = square_lattice_graph(mask)
    end
    if maximal
        gp = MaximalIndependence(g; optimizer=TreeSA(sc_target=sc_target, niters=5))
    else
        gp = Independence(g; optimizer=TreeSA(sc_target=sc_target, sc_weight=1.0, ntrials=3, βs=0.01:0.05:25.0, niters=20, rw_weight=1.0), simplifier=MergeGreedy())
        #gp = Independence(g)
    end
    println("Graph $graph of size $n, seed= $seed, usecuda = $usecuda")
    res = graph_polynomial(gp, Val(:finitefield), usecuda=usecuda)[]
    ofname = joinpath(folderout, "$(seed).dat")
    writefile && writedlm(ofname, res.coeffs)
    return nothing
end
end

function multiprocess_run(func, inputs::AbstractVector{T}) where T
    n = length(inputs)
    jobs = RemoteChannel(()->Channel{T}(n));
    results = RemoteChannel(()->Channel{Any}(n));
    for i in 1:n
        put!(jobs, inputs[i])
    end
    for p in workers() # start tasks on the workers to process requests in parallel
        remote_do(do_work, p, func, jobs, results)
    end
    return Any[take!(results) for i=1:n]
end

# patch
const graph = ARGS[1]

for L = 19:22
    println("computing L = $L")
    folder = joinpath(homedir(), ".julia/dev/TropicalMIS", "project", "data")
    fname = joinpath(folder, "mis_degeneracy_L$L.dat")
    maximal = false
    if maximal
        folderout = joinpath(folder, "$graph", "maximal_polynomial_L$(L)")
    else
        folderout = joinpath(folder, "$graph", "independence_polynomial_L$(L)")
    end
    if !isdir(dirname(folderout))
        mkdir(dirname(folderout))
    end
    if !isdir(folderout)
        mkdir(folderout)
    end
    masks = readdlm(fname)[:,4:end]
    multiprocess_run(collect(0:999)) do i
        mask = Matrix{Bool}(reshape(masks[i+1,:], L, L))
        @time mis_counting(graph, L, i, mask; writefile=true, sc_target=26, usecuda=false, maximal=maximal, folderout=folderout)
    end
end
