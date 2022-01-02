using Distributed, CUDA, GraphTensorNetworks, Random, DelimitedFiles, GraphTensorNetworks.OMEinsumContractionOrders, GraphTensorNetworks.OMEinsum
USECUDA = parse(Bool, ARGS[1])
@show USECUDA
const gpus = collect(devices())
println("find $(length(gpus)) GPU devices")
const procs = addprocs(length(gpus))
const process_device_map = Dict(zip(procs, gpus))
@show process_device_map

@everywhere begin
using GraphTensorNetworks, Random, GraphTensorNetworks.OMEinsumContractionOrders, GraphTensorNetworks.OMEinsum
using CUDA
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

function sequencing(graph, n; writefile, sc_target, usecuda, nslices=1, process_device_map)
    if graph == "diag"
        g = diagonal_coupled_graph(trues(n, n))
    elseif graph == "square"
        g = square_lattice_graph(trues(n, n))
    end
    gp = Independence(g; optimizer=TreeSA(sc_target=sc_target, sc_weight=1.0, nslices=nslices;
        ntrials=4, βs=0.01:0.05:25.0, niters=50, rw_weight=1.0), simplifier=MergeGreedy())
    println("Graph size $n, usecuda = $usecuda")
    @show timespace_complexity(gp)
    xs = GraphTensorNetworks.generate_tensors(x->1.0, gp)
    @time res = multigpu_contract(gp.code, (xs...,); process_device_map=process_device_map)
    @show res
    ofname = joinpath(@__DIR__, "data", "$graph", "entropyconstant-$n.dat")
    writefile && writedlm(ofname, res)
end

# current best for 38 = 204: 49.27
# current best for 39 = 7: 49.27
for L=27:38
    Random.seed!(5)
    graph = "diag"
    println("computing L = $L")
    @time sequencing(graph, L; writefile=true, sc_target=28, usecuda=USECUDA, nslices=L-26, process_device_map=process_device_map)
end
