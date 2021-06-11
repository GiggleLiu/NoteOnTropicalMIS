using DelimitedFiles, NoteOnTropicalMIS
using NoteOnTropicalMIS.OMEinsumContractionOrders
using NoteOnTropicalMIS.OMEinsum

function mis_configurations(n, seed, round; writefile=false)
    folder = joinpath("/home/leo/.julia/dev/TropicalMIS/", "project", "data")
    if round == 'A'
        fname = joinpath(folder, "atom_configuration_$(n)x$(n)_$seed.txt")
    elseif round == 'B'
        fname = joinpath(folder, "MISGraph_$(n)x$(n)_5.3x5.3_$seed.txt")
    end
    mask = Matrix{Bool}(readdlm(fname))
    rawcode = diagonal_coupled_eincode(mask)
    code = optimize_kahypar(rawcode, uniformsize(rawcode, 2); sc_target=17, max_group_size=40)
    nc = mis_count(code)[]
    println("Graph $seed, n = $n, MIS size = $(nc.n), degeneracy = $(nc.c)")
    if nc.c < 10000000
        res = mis_config(code; all=true)[]
        @show res.n, nc.n
        @assert res.n == nc.n
        config = res.c
        @assert length(config) == nc.c

        ofname = joinpath(folder, "configurations_$(n)x$(n)_5.3x5.3_$seed.dat")
        writefile && write(ofname, Matrix(config))
        return nc.n, config
    end
end

if false
    mis_configurations(15, 1, 'A')
    mis_configurations(15, 3, 'A')
    mis_configurations(13, 1, 'A')
    mis_configurations(13, 3, 'A')
else
    for n in [7, 10], seed in 1:6
        mis_configurations(n, seed, 'B'; writefile=true)
    end
    #for n in [13, 15], seed in [1, 2, 3, "Easy_1", "Medium_1", "Hard_1"]
    #    mis_configurations(n, seed, 'B', writefile=true)
    #end
end