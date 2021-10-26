using DelimitedFiles, GraphTensorNetworks
#using GraphTensorNetworks.OMEinsumContractionOrders
#using GraphTensorNetworks.OMEinsum
#using GraphTensorNetworks.TropicalNumbers

function mis_configurations(n, seed; writefile, sc_target=12)
    folder = joinpath(homedir(), ".julia/dev/TropicalMIS", "project", "data")
    fname = joinpath(folder, "mis_degeneracy_L$n.dat")
    mask = Matrix{Bool}(reshape(readdlm(fname)[seed+1,4:end], n, n))
    g = diagonal_coupled_graph(mask)
    code = Independence(g; optimizer=TreeSA(sc_target=sc_target))
    s, c = max_size_count(code)
    println("Graph $seed, n = $n, MIS size = $s, degeneracy = $c")
    if c > 1000000
        @warn "degeneracy too high, got: $c"
    end
    #s2, config1, config0 = (res = best2_solutions(code; all=true)[]; (res.maxorder, res.coeffs...))
    s2, config2, config1, config0 = (res = solve(code, "configs max3")[]; (res.maxorder, res.coeffs...))
    @assert s2 == s2
    @assert length(config0) == c
    ofname = joinpath(folder, "configurations_$(n)x$(n)_$seed.dat")
    ofname1 = joinpath(folder, "configurations_suboptimal_$(n)x$(n)_$seed.dat")
    ofname2 = joinpath(folder, "configurations_subsuboptimal_$(n)x$(n)_$seed.dat")
    writefile && write(ofname, toMatrix(config0))
    writefile && write(ofname1, toMatrix(config1))
    writefile && write(ofname2, toMatrix(config2))
    return s, config0, config1, config2
end

function toMatrix(x::ConfigEnumerator{N,C}) where {N,C}
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
    (8, [188])
    ]
    for seed in seeds
        @time mis_configurations(n, seed; writefile=true, sc_target=14)
    end
end