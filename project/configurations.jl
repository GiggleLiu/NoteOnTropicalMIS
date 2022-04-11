using DelimitedFiles, GenericTensorNetworks

function mis_configurations(n, seed; order, writefile, sc_target=12)
    folder = joinpath(homedir(), ".julia/dev/TropicalMIS", "project", "data", "configurations_L$(n)")
    if !isdir(folder)
        mkdir(folder)
    end
    fname = joinpath(dirname(folder), "mis_degeneracy_L$n.dat")
    mask = Matrix{Bool}(reshape(readdlm(fname)[seed+1,4:end], n, n))
    g = diagonal_coupled_graph(mask)
    gp = Independence(g; optimizer=GreedyMethod(), simplifier=MergeGreedy())
    res_c = solve(gp, "counting max$order")[]
    println("Graph $seed, n = $n, max2 count = $res_c")
    if sum(res_c.coeffs) > 10000000
        @warn "degeneracy too high!"
        return res_c, nothing
    end
    gp = Independence(g; optimizer=TreeSA(sc_target=sc_target), simplifier=MergeGreedy())
    s2, configs = (res = solve(gp, "configs max$order (bounded)")[]; (res.maxorder, res.coeffs))
    @assert all(res_c.coeffs .== length.(configs))
    for i=1:order
        ofname = joinpath(folder, "$seed-$(i-1).dat")
        writefile && write(ofname, toMatrix(configs[order-i+1]))
    end
    return res_c, configs
end

function toMatrix(x::ConfigEnumerator{N,K,C}) where {N,K,C}
    m = zeros(UInt64, C, length(x))
    for i=1:length(x)
        m[:,i] .= x.data[i].data
    end
    return m
end

for (n, seeds) in [
    #(5, [410, 407, 396]),
    #(6, [667, 557, 78]),
    #(7, [189, 623, 354]),
    #(10, [983, 828, 61]),
    #(11, [571, 808, 438]),
    #(15, [612, 907, 758])
    #(10, [982, 456, 843]),  # easy, ~25
    #(10, [347, 793, 82]),  # intermediate, ~85
    #(8, [188, 970, 91]),
    #(8, [100, 72]),
    #(8, [188, 216]),
    #(9, [805, 144, 560, 651, 234]),
    #(8, [188, 970, 91, 100, 72, 316, 747, 216, 168, 852,  7, 743, 32,
    #   573, 991, 957, 555, 936, 342, 950]),
    #(9, [805, 144, 560, 651, 234, 442, 408, 866, 263, 873, 91, 99, 566,
       #849, 854, 57, 64, 953, 43, 210])
    #(8, [100, 72, 316, 747, 216, 168, 852, 7, 743, 32, 573, 991, 957, 555, 936, 342, 950]),
    #(9, [805, 144, 560, 651, 234, 442, 408, 866, 263, 873, 91, 99, 566, 849, 854, 57, 64, 953, 43, 210])
    #(7, [189,623,354,40,323,173,661,345,813,35,162,965,336,667,870,1,156,901,576,346])
    #(7, [40, 323, 173, 661, 345, 813, 35, 162, 965, 336, 667, 870, 1, 156, 901, 576, 346])
    #(6, 807)
    #(6, [667, 557, 78, 312, 807, 776, 485, 980, 71, 50, 521, 773, 549, 523, 374, 515, 669, 344, 21, 107, 201,
    #    851, 736, 508, 286, 526, 385, 116, 20, 999, 357, 149, 872, 233, 528, 603, 912, 820])
    #(10, [828, 61, 156, 773, 105, 222, 417, 269, 309, 350, 786, 590, 109, 83, 243, 699, 425, 174, 925])
    #(11, [571, 808, 438, 748, 802, 454, 13, 401, 596, 126, 412, 977, 645, 263, 208, 622, 971, 725, 328, 895])
    #(10, 81:99)
    #(11, 0:99)
    #(12, 51:99)
    (10, [156, 773, 105, 222, 417, 269, 309, 350, 786, 590, 109, 83, 243, 699, 425, 174, 925])
    ]
    for seed in seeds
        @time mis_configurations(n, seed; order=2, writefile=true, sc_target=14)
    end
end