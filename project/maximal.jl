using NoteOnTropicalMIS, Polynomials, Random, OMEinsum

#Random.seed!(2)
#g = diagonal_coupled_graph(rand(13,13) .<= 0.8)
#@time maximal_polynomial(Val(:fft), g; sc_target=24, imbalances=0.0:0.00151:1.0, max_group_size=50)
#@time maximal_polynomial(Val(:polynomial), g; sc_target=19, imbalances=0.0:0.00151:1.0, max_group_size=50)

using DelimitedFiles
using NoteOnTropicalMIS.OMEinsum

function mis_maximal_counting(n, seed, round; writefile=false, rs, sc_target=24)
    folder = joinpath("/home/leo/.julia/dev/TropicalMIS/", "project", "data")
    if round == 'A'
        fname = joinpath(folder, "atom_configuration_$(n)x$(n)_$seed.txt")
    elseif round == 'B'
        fname = joinpath(folder, "MISGraph_$(n)x$(n)_5.3x5.3_$seed.txt")
    end
    mask = Matrix{Bool}(readdlm(fname))
    g = diagonal_coupled_graph(mask)
    println("Graph $seed")
    for r in rs
        res = maximal_polynomial(Val(:fft), g; sc_target=sc_target, imbalances=0.0:0.002:1.0, max_group_size=50, r=r)
        ofname = joinpath(folder, "maximal_counting_$(n)x$(n)_5.3x5.3_$(seed)_r$(r).dat")
        writefile && writedlm(ofname, res.coeffs)
    end
end

if false
    mis_maximal_counting(15, 1, 'A')
    mis_maximal_counting(15, 3, 'A')
    mis_maximal_counting(13, 1, 'A')
    mis_maximal_counting(13, 3, 'A')
else
    #for n in [7, 10], seed in 1:6
    #    mis_maximal_counting(n, seed, 'B'; writefile=true, rs=[1.0])
    #end
    for n in [13], seed in [1, 2, 3, "Easy_1", "Medium_1", "Hard_1"]
        mis_maximal_counting(n, seed, 'B', writefile=true)
    end
    #for n in [13, 15], seed in [1, 2, 3, "Easy_1", "Medium_1", "Hard_1"]
    #    mis_maximal_counting(n, seed, 'B', writefile=true)
    #end
    #for n in [15], seed in ["Wang"]
    #    mis_maximal_counting(n, seed, 'B', writefile=true)
    #end
end