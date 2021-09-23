using GraphTensorNetworks, DelimitedFiles
using LightGraphs
using BenchmarkTools, Random
using LinearAlgebra
using CUDA

BLAS.set_num_threads(1)

function case_r3(n, k=3; sc_target, seed=2)
    # generate a random regular graph of size 100, degree 3
    graph = (Random.seed!(seed); LightGraphs.random_regular_graph(n, k))
    @assert length(connected_components(graph)) == 1  # connected graph
    # optimize the contraction order using KaHyPar + Greedy
    #optcode = Independence(graph; optmethod=:kahypar, sc_target=sc_target, max_group_size=40, imbalances=0:0.001:1)
    optcode = Independence(graph; optmethod=:tree, sc_target=sc_target, sc_weight=2.0, ntrials=20, βs=0.01:0.05:15.0, niters=50, rw_weight=0.2)
    return optcode
end

function case_dc(L::Int, ρ; sc_target, seed=2)
    # generate a random regular graph of size 100, degree 3
    Random.seed!(seed)
    graph = diagonal_coupled_graph(rand(L, L) .< ρ)
    # optimize the contraction order using KaHyPar + Greedy, target space complexity is 2^20
    optcode = Independence(graph; optmethod=:kahypar, sc_target=sc_target, max_group_size=40)
    return optcode
end

function case_sq(L::Int, ρ; sc_target, seed=2)
    # generate a random regular graph of size 100, degree 3
    Random.seed!(seed)
    graph = square_lattice_graph(rand(L, L) .< ρ)
    # optimize the contraction order using KaHyPar + Greedy, target space complexity is 2^20
    optcode = idp_code(graph; method=:kahypar, sc_target=sc_target, max_group_size=40)
    return optcode
end

# setup global arguments
const TASK = length(ARGS) >= 1 ? ARGS[1] : "size_max"
const DEVICE = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : -1

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

function runcase(cases; task = :maxsize, usecuda = false)

    run_benchmarks([("n$(10*i)", ()->(usecuda ? (CUDA.@sync solve(case, replace(task, "_"=>" "); usecuda=true)) : solve((@show length(GraphTensorNetworks.labels(case.code)); case), replace(task, "_"=>" "); usecuda=false))) for (i, case) in enumerate(cases)],
                   output_file=joinpath(@__DIR__, "data", "$(task)-r3-$(usecuda ? "GPU" : "CPU").dat"))
end

const truncatedict = Dict([string(task)=>ntruncate for (task, ntruncate) in [
        ("counting_sum", 0), ("size_max", 0), ("counting_max", 0), ("counting_max2", 0),
        ("counting_all", 3), ("counting_all_(fft)", 0), ("counting_all_(finitefield)", 0),
        ("config_max", 0), ("configs_max",5), ("configs_all", 16), ("configs_max2", 9), ("config_max_(bounded)", 0), ("configs_max_(bounded)", 0)
        ]])

cases = [case_r3(n, 3; seed=2, sc_target=s) for (n, s) in [
        (10, 3), (20, 4), (40, 5), (50, 8), (60, 8), (70, 8), (80, 10), (90, 13), (100, 13),
        (110, 15), (120, 15), (130, 13), (140, 17), (150, 18), (160, 20), (170, 19), (180, 24), (190, 24), (200, 25),
       ]]

function run_cpu(truncatedict)
    run = false
    for TASK in keys(truncatedict)
        println(TASK)
        run && runcase(cases[1:end-truncatedict[TASK]]; task=TASK, usecuda=false)
        if TASK == "counting_max2"
	    run = true
	end
    end
end

function run_gpu(truncatedict)
    for TASK in ["counting_sum", "size_max", "counting_max", "counting_max2",
        "counting_all_(fft)", "counting_all_(finitefield)",
        "config_max", "config_max_(bounded)"
        ]
        runcase(cases[1:end-truncatedict[TASK]]; task=TASK, usecuda=true)
    end
end


if DEVICE >= 0
    run_gpu(truncatedict)
else
    run_cpu(truncatedict)
end
