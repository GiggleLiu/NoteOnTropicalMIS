using GenericTensorNetworks, Random
using Distributed, DelimitedFiles, GenericTensorNetworks.OMEinsumContractionOrders, GenericTensorNetworks.OMEinsum

@everywhere begin
using GenericTensorNetworks, Random, GenericTensorNetworks.OMEinsumContractionOrders, GenericTensorNetworks.OMEinsum
using DelimitedFiles

function do_work(f, jobs, results) # define work function everywhere
    while true
        job = take!(jobs)
        @info "running $job on device $(Distributed.myid())"
        res = f(job)
        put!(results, res)
    end
end

function mis_counting(graph, n, seed, mask; writefile, sc_target, usecuda, maximal, folderout, overwrite)
    # avoid overwriting
    ofname = joinpath(folderout, "$(seed).dat")
    if !overwrite && isfile(ofname)
        return nothing
    end

    if graph == "diag"
        g = diagonal_coupled_graph(mask)
    elseif graph == "square"
        g = square_lattice_graph(mask)
    end
    if maximal
        gp = MaximalIndependence(g; optimizer=TreeSA(sc_target=sc_target, niters=5))
    else
        gp = Independence(g; optimizer=TreeSA(sc_target=sc_target, sc_weight=1.0, ntrials=1, βs=0.01:0.1:30.0, niters=10, rw_weight=1.0), simplifier=MergeGreedy())
    end
    println("Graph $graph of size $n, seed= $seed, usecuda = $usecuda")
    res = solve(gp, GraphPolynomial(); usecuda=usecuda)[]

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

const LS = 21:24
const folderout_list = String[]
const masks_list = []
const maximal = false

for L = LS
    println("computing L = $L")
    folder = joinpath(homedir(), ".julia/dev/TropicalMIS", "project", "data")
    fname = joinpath(folder, "mis_degeneracy2_L$L.dat")
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
    masks = readdlm(fname)[:,5:end]
    push!(folderout_list, folderout)
    push!(masks_list, masks)
end

const args = vcat([tuple.(L, Ref(folderout), collect(0:999), [Matrix{Bool}(reshape(masks[k,:], L, L)) for k=1:1000]) for (masks, folderout, L) in zip(masks_list, folderout_list, LS)]...)

multiprocess_run(args) do (L, folderout, i, mask)
    @time mis_counting(graph, L, i, mask; writefile=true, sc_target=26, usecuda=false, maximal=maximal, folderout=folderout, overwrite=false)
end
