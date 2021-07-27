using NoteOnTropicalMIS, Polynomials, Random, OMEinsum

using DelimitedFiles
using NoteOnTropicalMIS.OMEinsum

function mis_maximal_counting(n, seed; writefile, rs=[1.0], sc_target=24, method=:finitefield)
    folder = joinpath("/home/leo/.julia/dev/TropicalMIS", "project", "data")
    fname = joinpath(folder, "mis_degeneracy_L$n.dat")
    mask = Matrix{Bool}(reshape(readdlm(fname)[seed+1,4:end], n, n))
    g = diagonal_coupled_graph(mask)
    println("Graph $seed")
    if method == :fft
        for r in rs
            res = maximal_polynomial(Val(:fft), g; sc_target=sc_target, imbalances=0.0:0.002:1.0, max_group_size=50, r=r)
            ofname = joinpath(folder, "maximal_counting_$(n)x$(n)_$(seed)_r$(r).dat")
            writefile && writedlm(ofname, res.coeffs)
        end
    elseif method == :finitefield
        res = maximal_polynomial(Val(:finitefield), g; sc_target=sc_target, imbalances=0.0:0.002:1.0, max_group_size=50)
        ofname = joinpath(folder, "maximal_counting_$(n)x$(n)_$(seed).dat")
        writefile && writedlm(ofname, res.coeffs)
    end
end

for seed in [0, 2, 4, 9, 50]
    @time mis_maximal_counting(10, seed; writefile=true, sc_target=20)
end