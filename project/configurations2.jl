using DelimitedFiles, NoteOnTropicalMIS
using NoteOnTropicalMIS.OMEinsumContractionOrders
using NoteOnTropicalMIS.OMEinsum
using NoteOnTropicalMIS.TropicalNumbers

function mis_configurations(n, seed; writefile, sc_target=12)
    folder = joinpath("/home/leo/.julia/dev/TropicalMIS", "project", "data")
    fname = joinpath(folder, "mis_degeneracy_L$n.dat")
    mask = Matrix{Bool}(reshape(readdlm(fname)[seed+1,4:end], n, n))
    g = diagonal_coupled_graph(mask)
    code = NoteOnTropicalMIS.idp_code(g; method=:kahypar, sc_target=sc_target, max_group_size=50, imbalances=0.0:0.003:0.8)
    s, c = (res = mis_contract(CountingTropical{Int,Int}(1,1), code)[]; (res.n, res.c))
    println("Graph $seed, n = $n, MIS size = $s, degeneracy = $c")
    if c > 1000000
        @warn "degeneracy too high, got: $c"
    end
    s2, config = (res = mis_config(code; all=true)[]; (res.n, res.c))
    @assert s == s2
    @assert length(config) == c
    ofname = joinpath(folder, "configurations_$(n)x$(n)_$seed.dat")
    writefile && write(ofname, toMatrix(config))
    return s, config
end

function toMatrix(x::ConfigEnumerator{N,C}) where {N,C}
    m = zeros(UInt64, C, length(x))
    for i=1:length(x)
        m[:,i] .= x.data[i].data
    end
    return m
end

for seed in [0, 2, 4, 9, 50]
    @time mis_configurations(10, seed; writefile=true)
end

for seed in  [4, 6, 7, 17] #[26, 40, 339, 339, 8, 13, 25, 89, 108, 138]
    @time mis_configurations(15, seed; writefile=true, sc_target=16)
end