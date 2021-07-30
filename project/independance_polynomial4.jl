using DelimitedFiles, NoteOnTropicalMIS
using NoteOnTropicalMIS: unitdisk_graph

const DEVICE = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : -1

if DEVICE >= 0
    using CUDA
    CUDA.device!(DEVICE)
end

function load_graph()
    file = joinpath(homedir(), ".julia", "dev", "TropicalMIS", "project", "data", "atoms-leo.dat")
    atoms = readdlm(file)
    unitdisk_graph([atoms[i,:] for i=1:size(atoms, 1)], 1.1)
end

function count_degeneracy(; write, usecuda=false, r=8.0)
    code = NoteOnTropicalMIS.idp_code(load_graph(); method=:kahypar, imbalances=0.0:0.002:1.0, sc_target=19)
    #n, c = (res = run_task(code, :counting)[]; (res.n, res.c))
    idp = independence_polynomial(Val(:fft), code; usecuda=usecuda, r=r)
    #println("Graph MIS size = $n, degeneracy = $(idp.coeffs[end]) (should be same as $c)")
    if write
        writedlm(joinpath(homedir(), ".julia", "dev", "TropicalMIS", "project", "data", "independance_polynomial_leo_r$r.dat"), idp)
    end
    return idp
end

#count_degeneracy(write=true, usecuda=DEVICE>=0, r=8.0)
count_degeneracy(write=true, usecuda=DEVICE>=0, r=0.1)
