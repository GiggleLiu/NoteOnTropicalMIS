using DelimitedFiles, NoteOnTropicalMIS
using NoteOnTropicalMIS: unitdisk_graph

const DEVICE = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : -1

if DEVICE >= 0
    using CUDA, CuYao
    CUDA.device!(DEVICE)
end

function load_graph()
    file = joinpath(homedir(), ".julia", "dev", "TropicalMIS", "project", "data", "atoms-leo.dat")
    atoms = readdlm(file)
    unitdisk_graph([atoms[i,:] for i=1:size(atoms, 1)], 1.1)
end

function count_degeneracy(; write, usecuda=false)
    code = NoteOnTropicalMIS.idp_code(load_graph(); method=:kahypar, imbalances=0.0:0.002:1.0, sc_target=19)
    n, c = (res = run_task(code, :counting)[]; (res.n, res.c))
    idp = run_task(code, :idp_finitefield; usecuda=usecuda)
    println("Graph MIS size = $n, degeneracy = $(idp.coeffs[end]) (should be same as $c)")
    if write
        writedlm(joinpath(homedir(), ".julia", "dev", "TropicalMIS", "project", "data", "independance_polynomial_leo.dat"), idp)
    end
    return idp
end

count_degeneracy(write=false)