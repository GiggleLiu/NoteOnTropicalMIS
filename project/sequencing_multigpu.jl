using Distributed, CUDA, GraphTensorNetworks, Random, DelimitedFiles, GraphTensorNetworks.OMEinsumContractionOrders, GraphTensorNetworks.OMEinsum
USECUDA = parse(Bool, ARGS[1])
@show USECUDA
if USECUDA
    println("find $(length(devices())) GPU devices")
    addprocs(length(devices()))
else
    addprocs(4)
end

@everywhere begin
using GraphTensorNetworks, Random, GraphTensorNetworks.OMEinsumContractionOrders, GraphTensorNetworks.OMEinsum
using CUDA
CUDA.allowscalar(false)
using DelimitedFiles
const process_device_map = Dict(zip(procs, gpus))
@show process_device_map

function do_work(f, jobs, results) # define work function everywhere
    while true
        job = take!(jobs)
        @info "running $job on device $(Distributed.myid())"
        res = f(job)
        put!(results, res)
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
end

function sequencing(n; writefile, sc_target, usecuda, nslices=1, process_device_map)
    g = square_lattice_graph(trues(n, n))
    gp = Independence(g; optimizer=TreeSA(sc_target=sc_target, sc_weight=1.0, nslices=nslices;
        ntrials=4, βs=0.01:0.05:25.0, niters=50, rw_weight=1.0), simplifier=MergeGreedy())
    println("Graph size $n, usecuda = $usecuda")
    @show timespace_complexity(gp)
    res = GraphTensorNetworks.big_integer_solve(Int32, 100) do T
        @info "T = $T"
        filename = joinpath(@__DIR__, "data", "$n-$(GraphTensorNetworks.modulus(one(T))).dat")
        if isfile(filename)
            fill(T(readdlm(filename)[]))
        else
            xs = GraphTensorNetworks.generate_tensors(one(T), gp)
            @time res = multigpu_contract(gp.code, xs; process_device_map=process_device_map)
            writedlm(filename, res[].val)
            res
        end
    end
    @show res
    ofname = joinpath(@__DIR__, "data", "$n.dat")
    writefile && writedlm(ofname, res)
end

# current best = 204: 49.27
Random.seed!(204) #parse(Int, ARGS[2]))
for L=38
    println("computing L = $L")
    @time sequencing(L; writefile=true, sc_target=28, usecuda=USECUDA, nslices=L-28, process_device_map=process_device_map)
end
