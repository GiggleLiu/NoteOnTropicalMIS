using GraphTensorNetworks, DelimitedFiles
using LightGraphs
using BenchmarkTools, Random
using LinearAlgebra
using CUDA
using OMEinsumContractionOrders

BLAS.set_num_threads(1)

function case_r3(n, k=3; sc_target, seed=2)
    # generate a random regular graph of size 100, degree 3
    graph = (Random.seed!(seed); LightGraphs.random_regular_graph(n, k))
    @assert length(connected_components(graph)) == 1  # connected graph
    # optimize the contraction order using KaHyPar + Greedy
    #optcode = MaximalIndependence(graph; optmethod=:kahypar, sc_target=sc_target, max_group_size=20, imbalances=0:0.001:1)
    optcode = MaximalIndependence(graph; optmethod=:tree, sc_target=sc_target, sc_weight=1.0, ntrials=20, βs=0.02:0.05:15.0, niters=50, initializer=:greedy)
    return optcode
end

# setup global arguments
const GRAPH = length(ARGS) >= 1 ? ARGS[1] : "r3"
const TASK = length(ARGS) >= 2 ? ARGS[2] : "counting_sum"
const DEVICE = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : -1

if DEVICE >= 0
    CUDA.device!(DEVICE)
    Base.ndims(::Base.Broadcast.Broadcasted{CUDA.CuArrayStyle{0}}) = 0
end

function run_benchmarks(cases; output_file)
    suite = BenchmarkGroup()
    for (case, f) in cases
        suite[case] = @benchmarkable $f()
    end

    tune!(suite)
    res = run(suite)

    times = zeros(length(cases))
    for (k, (case, f)) in enumerate(cases)
        times[k] = minimum(res[case].times)
    end

    println("Writing benchmark results to file: $output_file")
    mkpath(dirname(output_file))
    writedlm(output_file, times)
end

function runcase(;
        case_set=:r3,
        task,
        usecuda = false,
        ntruncate = 0,  # truncate benchmark cases
        seed=2
    )
    if case_set == :r3
        cases = [case_r3(n, 3; seed=seed, sc_target=s) for (n, s) in [
            (10, 6), (20, 8), (30, 9), (40, 11), (50, 16), (60, 17), (70, 16), (80, 20), (90, 26), (100, 26),
        ][1:end-ntruncate]]
        run_benchmarks([("n$(10*i)", ()->(usecuda ? (CUDA.@sync solve(case, replace(task, "_"=>" "); usecuda=true)) : solve((@show length(GraphTensorNetworks.labels(case.code)); case), replace(task, "_"=>" "); usecuda=false))) for (i, case) in enumerate(cases)],
                    output_file=joinpath(@__DIR__, "data", "maximal-$(task)-$(case_set)-$(usecuda ? "GPU" : "CPU").dat"))
    elseif case_set == :r3bk
        graphs = [(Random.seed!(seed); LightGraphs.random_regular_graph(i*10, 3)) for i=1:10-ntruncate]
        run_benchmarks([("n$(10*i)", ()->maximal_cliques(complement(g))) for (i, g) in enumerate(graphs)],
                    output_file=joinpath(@__DIR__, "data", "maximal-$(task)-$(case_set)-$(usecuda ? "GPU" : "CPU").dat"))
    end

end

const truncatedict = Dict(
    "r3"=>Dict([string(task)=>ntruncate for (task, ntruncate) in [
        ("counting_sum", 0), ("counting_all_(finitefield)", 0), ("configs_all", 3),
        ]]),
    "r3bk"=>Dict(["configs_all"=>4])
    )

runcase(case_set=Symbol(GRAPH), task=TASK, usecuda=DEVICE>=0, ntruncate=truncatedict[GRAPH][TASK])
