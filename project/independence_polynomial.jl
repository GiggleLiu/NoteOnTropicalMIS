using NoteOnTropicalMIS, Polynomials, Random, OMEinsum
using DelimitedFiles
using NoteOnTropicalMIS.OMEinsum

function idp(which, n::Int; writefile=false, rs, sc_target=24)
    if which == :diagonal_coupled
        g = diagonal_coupled_graph(trues(n, n))
    elseif which == :square_lattice
        g = square_lattice_graph(trues(n, n))
    end
    println("Graph n=$n")
    optcode = idp_code(g; method=:kahypar, sc_target=sc_target, imbalances=0.0:0.002:1.0, max_group_size=50)
    for r in rs
        res = independence_polynomial(Val(:fft), optcode; r=r)
        println(res)
        #ofname = joinpath(folder, "maximal_counting_$(n)x$(n)_5.3x5.3_$(seed)_r$(r).dat")
        #writefile && writedlm(ofname, res.coeffs)
    end
end

idp(:diagonal_coupled, 10; rs=[1.0])