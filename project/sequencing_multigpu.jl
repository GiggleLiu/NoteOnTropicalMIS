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
const DEVICES = collect(devices())
using DelimitedFiles

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

function (se::SlicedEinsum{LT,ET})(@nospecialize(xs::AbstractArray...); size_info = nothing) where {LT, ET}
    xs = Array.(xs)  # convert to CPU first, then back to GPU
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
        device!(DEVICES[Distributed.myid()-1])
        flush(stdout)
        xsi = ntuple(i->CuArray(OMEinsumContractionOrders.take_slice(xs[i], it.ixsv[i], slicemap)), length(xs))
        Array(einsum(eins_sliced, xsi, it.size_dict_sliced))
    end
    @show results
    for (resi, (k, slicemap)) in zip(results, inputs)
        OMEinsumContractionOrders.fill_slice!(res, it.iyv, resi, slicemap)
    end
    return res
end
end

function sequencing(n; writefile, sc_target, usecuda, nslices=1)
    g = square_lattice_graph(trues(n, n))
    gp = Independence(g; optimizer=TreeSA(sc_target=sc_target, sc_weight=1.0, nslices=nslices;
        ntrials=10, βs=0.01:0.05:22.0, niters=20, rw_weight=2.0), simplifier=MergeGreedy())
    println("Graph size $n, usecuda = $usecuda")
    @show timespace_complexity(gp)
    res = GraphTensorNetworks.big_integer_solve(Int32, 100) do T
        @time Array(GraphTensorNetworks.contractx(gp, one(T); usecuda=usecuda))
    end
    @show res
    ofname = joinpath(@__DIR__, "data", "$n.dat")
    writefile && writedlm(ofname, res)
end

Random.seed!(2)
for L=31
    println("computing L = $L")
    @time sequencing(L; writefile=true, sc_target=29, usecuda=USECUDA, nslices=L-29)
end
