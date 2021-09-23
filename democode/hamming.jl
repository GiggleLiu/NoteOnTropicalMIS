using GraphTensorNetworks, Random, LightGraphs
using DelimitedFiles

function case_r3(n, k=3; sc_target, seed=2)
    # generate a random regular graph of size 100, degree 3
    graph = (Random.seed!(seed); LightGraphs.random_regular_graph(n, k))
    @assert length(connected_components(graph)) == 1  # connected graph
    optcode = Independence(graph; optmethod=:tree, sc_target=sc_target, sc_weight=1.0, ntrials=20, βs=0.01:0.05:15.0, niters=50, rw_weight=0.2)
    return optcode
end

function compute_configs(filename; seed)
    se = best_solutions(case_r3(200, 3; sc_target=25, seed=seed); all=true)[].c
    save_configs(filename, se, format=:text)
end

function analyse_hamming(filename)
    hamming = zeros(201)
    configs = load_configs(filename; format=:text, nflavors=2)
    for j=1:length(configs)
        for i=j+1:length(configs)
            hamming[count_ones(configs[i] ⊻ configs[j])+1] += 1
        end
    end
    return hamming
end

for seed = 4:10
    filename = joinpath(@__DIR__, "data", "r3_n200_optconfigs_$seed.txt")
    compute_configs(filename, seed=seed)
    filename2 = joinpath(@__DIR__, "data", "r3_n200_hamming_$seed.txt")
    writedlm(filename2, analyse_hamming(filename))
end
