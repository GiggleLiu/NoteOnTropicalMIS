export case_r3, case_dc, run_task

using Random

function case_r3(n, k=3; sc_target, seed=2)
    # generate a random regular graph of size 100, degree 3
    graph = (Random.seed!(seed); LightGraphs.random_regular_graph(n, k))
    # generate einsum code, i.e. the labels of tensors
    code = EinCode(([minmax(e.src,e.dst) for e in LightGraphs.edges(graph)]..., # labels for edge tensors
                    [(i,) for i in LightGraphs.vertices(graph)]...), ())        # labels for vertex tensors
    size_dict = Dict([s=>2 for s in symbols(code)])
    # optimize the contraction order using KaHyPar + Greedy
    optimized_code = optimize_kahypar(code, size_dict; sc_target=sc_target, max_group_size=40)
    println("time/space complexity is $(OMEinsum.timespace_complexity(optimized_code, size_dict))")
    return optimized_code
end

function case_dc(L::Int, ρ; sc_target, seed=2)
    # generate a random regular graph of size 100, degree 3
    Random.seed!(seed)
    graph = diagonal_coupled_graph(rand(L, L) .< ρ)
    # generate einsum code, i.e. the labels of tensors
    code = EinCode(([minmax(e.src,e.dst) for e in LightGraphs.edges(graph)]..., # labels for edge tensors
                    [(i,) for i in LightGraphs.vertices(graph)]...), ())        # labels for vertex tensors
    size_dict = Dict([s=>2 for s in symbols(code)])
    # optimize the contraction order using KaHyPar + Greedy, target space complexity is 2^20
    optimized_code = optimize_kahypar(code, size_dict; sc_target=sc_target, max_group_size=40)
    println("time/space complexity is $(OMEinsum.timespace_complexity(optimized_code, size_dict))")
    return optimized_code
end

function run_task(code, task; usecuda=false)
    if task == :totalsize
        return mis_contract(1.0, code; usecuda=usecuda)
    elseif task == :maxsize
        return mis_size(code; usecuda=usecuda)
    elseif task == :counting
        return mis_count(code; usecuda=usecuda)
    elseif task == :idp_polynomial
        return independence_polynomial(Val(:polynomial), code; usecuda=usecuda)
    elseif task == :idp_fft
        return independence_polynomial(Val(:fft), code; usecuda=usecuda)
    elseif task == :idp_finitefield
        return independence_polynomial(Val(:finitefield), code; usecuda=usecuda)
    elseif task == :config_single
        return mis_config(code; all=false, bounding=false, usecuda=usecuda)
    elseif task == :config_single_bounded
        return mis_config(code; all=false, bounding=true, usecuda=usecuda)
    elseif task == :config_all
        return mis_config(code; all=true, bounding=false, usecuda=usecuda)
    elseif task == :config_all_bounded
        return mis_config(code; all=true, bounding=true, usecuda=usecuda)
    else
        error("unknown task $task.")
    end
end