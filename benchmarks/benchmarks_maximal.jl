using GenericTensorNetworks, DelimitedFiles, TropicalGEMM
using Graphs
using BenchmarkTools, Random
using LinearAlgebra
using CUDA
using OMEinsumContractionOrders
using Comonicon

BLAS.set_num_threads(1)

function case_r3(n, k=3; maxsc, seed=2)
    # generate a random regular graph of size 100, degree 3
    graph = (Random.seed!(seed); Graphs.random_regular_graph(n, k))
    @assert length(connected_components(graph)) == 1  # connected graph
    # optimize the contraction order using KaHyPar + Greedy
    optcode = MaximalIndependence(graph; optimizer=TreeSA(sc_target=0, sc_weight=1.0, ntrials=10, βs=0.01:0.05:25.0, niters=20, rw_weight=2.0), simplifier=MergeGreedy())
    tw = timespace_complexity(optcode)[2]
    @info "n = $n, tw = $tw"
    if tw > maxsc
        optcode = MaximalIndependence(graph; optimizer=TreeSA(sc_target=maxsc, sc_weight=1.0, ntrials=10, βs=0.01:0.05:25.0, niters=20, rw_weight=2.0, nslices=Int(tw-maxsc)), simplifier=MergeGreedy())
        tw = timespace_complexity(optcode)[2]
        @info "sliced: n = $n, tw = $tw"
    end
    return optcode
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

const truncatedict = Dict([string(task)=>ntruncate for (task, ntruncate) in [
        ("counting_sum", 0), ("counting_all_(finitefield)", 0), ("configs_all", 3), ("bron", 4),
       ]])

@main function runcase(task::String, usecuda::Bool, maxsc::Int=27, seed::Int=8)
    if task ∈ ("tcsc", "treewidth")
        cases = [case_r3(n, 3; seed=seed, maxsc=task=="tcsc" ? maxsc : 100) for n=10:10:100]
        tcscs = zeros(length(cases), 2)
        for i=1:length(cases)
            tcscs[i,:] .= timespace_complexity(cases[i])
        end
        writedlm(joinpath(@__DIR__, "data", "maximal-$task.dat"), tcscs)
    elseif task == "bron"
        ntruncate=truncatedict[task]
        graphs = [(Random.seed!(seed); Graphs.random_regular_graph(i*10, 3)) for i=1:10-ntruncate]
        run_benchmarks([("n$(10*i)", ()->maximal_cliques(complement(g))) for (i, g) in enumerate(graphs)],
                    output_file=joinpath(@__DIR__, "data", "maximal-$(task)-r3-$(usecuda ? "GPU" : "CPU").dat"))
        cases = [case_r3(n, 3; seed=seed, maxsc=maxsc) for (n, s) in [
            (10, 6), (20, 8), (30, 9), (40, 11), (50, 16), (60, 17), (70, 16), (80, 20), (90, 26), (100, 26),
        ][1:end-ntruncate]]
    else
        ntruncate=truncatedict[task]
        run_benchmarks([("n$(10*i)", ()->(usecuda ? (CUDA.@sync solve(case, replace(task, "_"=>" "); usecuda=true)) : solve((@show length(GenericTensorNetworks.labels(case.code)); case), replace(task, "_"=>" "); usecuda=false))) for (i, case) in enumerate(cases)],
                    output_file=joinpath(@__DIR__, "data", "maximal-$(task)-r3-$(usecuda ? "GPU" : "CPU").dat"))
    end
end
