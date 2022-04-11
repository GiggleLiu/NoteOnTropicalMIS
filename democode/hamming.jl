using GenericTensorNetworks, Random, Graphs
using DelimitedFiles

function case_r3(n, k=3; sc_target, seed=2)
    # generate a random regular graph of size 100, degree 3
    graph = (Random.seed!(seed); Graphs.random_regular_graph(n, k))
    @assert length(connected_components(graph)) == 1  # connected graph
    optcode = Independence(graph; optimizer=TreeSA(sc_target=sc_target,
                sc_weight=1.0, ntrials=10, βs=0.01:0.05:25.0, niters=20, rw_weight=2.0))
    return optcode
end

function count_configs(n, seed)
    se = solve(case_r3(n, 3; sc_target=25, seed=seed), "counting max2")[].coeffs
    return se
end

function compute_configs(n; seed)
    se = solve(case_r3(n, 3; sc_target=25, seed=seed), "configs max2 (bounded)")[].coeffs
    filename = joinpath(@__DIR__, "data", "r3_n$(n)_optconfigs_$seed.dat")
    save_configs(filename, se[2], format=:binary)
    filename = joinpath(@__DIR__, "data", "r3_n$(n)_suboptconfigs_$seed.dat")
    save_configs(filename, se[1], format=:binary)
end

function analyse_hamming(n; seed)
    hamming = zeros(201)
    f1 = joinpath(@__DIR__, "data", "r3_n$(n)_optconfigs_$seed.dat")
    f2 = joinpath(@__DIR__, "data", "r3_n$(n)_suboptconfigs_$seed.dat")
    filename = joinpath(@__DIR__, "data", "r3_n$(n)_hamming_$seed.txt")
    C1 = load_configs(f1; format=:binary, nflavors=2, bitlength=n)
    C2 = load_configs(f2; format=:binary, nflavors=2, bitlength=n)
    configs = [C1.data..., C2.data...]
    for j=1:length(configs)
        for i=j+1:length(configs)
            hamming[count_ones(configs[i] ⊻ configs[j])+1] += 1
        end
    end
    writedlm(filename, hamming)
    return hamming
end

for seed = 2:3
    n = 100
    @show seed
    compute_configs(n, seed=seed)
    analyse_hamming(n; seed=seed)
end
