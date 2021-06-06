using DelimitedFiles, NoteOnTropicalMIS
using NoteOnTropicalMIS.OMEinsumContractionOrders
using NoteOnTropicalMIS.OMEinsum

function mis_configurations(n, seed; writefile=false)
    folder = joinpath("project", "data")
    mask = Matrix{Bool}(readdlm(fname))
    rawcode = diagonal_coupled_eincode(mask)
    code = optimize_kahypar(rawcode, uniformsize(rawcode, 2); sc_target=17, max_group_size=40)
    if task == :counting
        nc = mis_count(code)[]
        println("n = $n, MIS size = $(nc.n), degeneracy = $(nc.c)")
    elseif task == :config
    elseif task == :allconfigs
        res = mis_config(code; all=true)[]
        @assert res.n == nc.n
        config = res.c
        @assert length(config) == nc.c

        ofname = joinpath(folder, "configurations_$(n)x$(n)_5.3x5.3_$seed.dat")
        writefile && write(ofname, Matrix(config))
        return nc.n, config
    end
end

if true
    #for n in [7, 10], seed in 1:6
    #    mis_configurations(n, seed; writefile=true)
    #end
    for n in [13, 15], seed in [1, 2, 3, "Easy_1", "Medium_1", "Hard_1"]
        mis_configurations(n, seed, writefile=true)
    end
end