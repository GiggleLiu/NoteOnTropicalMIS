using Distributed, CUDA, GraphTensorNetworks, Random, DelimitedFiles, GraphTensorNetworks.OMEinsumContractionOrders, GraphTensorNetworks.OMEinsum, StatsBase
USECUDA = parse(Bool, ARGS[1])
@show USECUDA
const gpus = collect(devices())
println("find $(length(gpus)) GPU devices")
const procs = addprocs(length(gpus))
const process_device_map = Dict(zip(procs, gpus))
@show process_device_map

@everywhere begin
using GraphTensorNetworks, Random, GraphTensorNetworks.OMEinsumContractionOrders, GraphTensorNetworks.OMEinsum
using CUDA, StatsBase
CUDA.allowscalar(false)
using DelimitedFiles

function do_work(f, jobs, results) # define work function everywhere
    while true
        job = take!(jobs)
        @info "running $job on device $(Distributed.myid())"
        res = f(job)
        put!(results, res)
    end
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

function generate_mask(Nx::Int, Ny::Int, natoms::Int)
    mask = zeros(Bool,Nx, Ny)
    mask[sample(1:Nx*Ny, natoms; replace=false)] .= true
    return mask
end

function multigpu_contract(se::SlicedEinsum{LT,ET}, xs::Tuple; size_info = nothing, process_device_map::Dict) where {LT, ET}
    length(se.slicing) == 0 && return se.eins(xs...; size_info=size_info)
    size_dict = size_info===nothing ? Dict{OMEinsum.labeltype(se),Int}() : copy(size_info)
    OMEinsum.get_size_dict!(se, xs, size_dict)

    it = OMEinsumContractionOrders.SliceIterator(se, size_dict)
    res = OMEinsum.get_output_array(xs, getindex.(Ref(size_dict), it.iyv))
    eins_sliced = OMEinsumContractionOrders.drop_slicedim(se.eins, se.slicing)
    inputs = collect(enumerate([copy(x) for x in it]))
    @info "start multiple process contraction!"
    results = multiprocess_run(inputs) do (k, slicemap)
        @info "computing slice $k/$(length(it))"
        device!(process_device_map[Distributed.myid()])
        xsi = ntuple(i->CuArray(OMEinsumContractionOrders.take_slice(xs[i], it.ixsv[i], slicemap)), length(xs))
        Array(einsum(eins_sliced, xsi, it.size_dict_sliced))
    end
    for (resi, (k, slicemap)) in zip(results, inputs)
        OMEinsumContractionOrders.fill_slice!(res, it.iyv, resi, slicemap)
    end
    @show res
    return res
end

function countall(grid, n; sc_target, usecuda, nslices, process_device_map, mask)
    if grid == "diag"
        g = diagonal_coupled_graph(mask)
    elseif grid == "square"
        g = square_lattice_graph(mask)
    end
    if n >= 28
        gp = Independence(g; optimizer=TreeSA(sc_target=sc_target, sc_weight=1.0, nslices=nslices;
            ntrials=3, βs=0.01:0.05:25.0, niters=6, rw_weight=2.0), simplifier=MergeGreedy())
    else
        gp = Independence(g)
    end
    println("Graph size $n, usecuda = $usecuda")
    @show timespace_complexity(gp)
    xs = GraphTensorNetworks.generate_tensors(x->1.0, gp)
    if n>= 25
        @time res = multigpu_contract(gp.code, (xs...,); process_device_map=process_device_map)
    else
        @time res = gp.code(xs...)
    end
    return res[]
end

function countall_multiple(graph, n; writefile, sc_target, usecuda, nslices, process_device_map)
    if graph == "diag"
        grid = "diag"
        ρ = 1.0
        nrepeat = 1
    elseif graph == "square"
        grid = "square"
        ρ = 1.0
        nrepeat = 1
    elseif graph == "diag-0.8"
        grid = "diag"
        ρ = 0.8
        nrepeat = 1000
    elseif graph == "square-0.8"
        grid = "square"
        ρ = 0.8
        nrepeat = 1000
    end
    res = zeros(nrepeat)
    for i=1:nrepeat
        mask = generate_mask(n, n, round(Int,n^2*ρ))
        println("$i: graph $grid, n=$n, natoms=$(sum(mask))")
        res[i] = countall(grid, n; sc_target=sc_target, usecuda=usecuda, nslices=nslices, process_device_map=process_device_map, mask=mask)
    end
    ofname = joinpath(@__DIR__, "data", "$graph", "entropyconstant-$n.dat")
    !ispath(dirname(ofname)) && mkdir(dirname(ofname))
    writefile && writedlm(ofname, res)
end

# current best for 38 = 204: 49.27
# current best for 39 = 7: 49.27
for L=2:32
    Random.seed!(2)
    graph = "square-0.8"
    println("computing L = $L")
    @time countall_multiple(graph, L; writefile=true, sc_target=28, usecuda=USECUDA, nslices=0, process_device_map=process_device_map)
end
