using DelimitedFiles, NoteOnTropicalMIS
using NoteOnTropicalMIS.OMEinsumContractionOrders
using NoteOnTropicalMIS.OMEinsum
using TensorOperations, LinearAlgebra

function LinearAlgebra.permutedims!(C::Array{T,N}, A::StridedArray{T,N}, perm) where {T,N}
    if isbitstype(T)
        TensorOperations.tensorcopy!(A, ntuple(identity,N), C, perm)
    else
        invoke(permutedims!, Tuple{Any,AbstractArray,Any}, C, A, perm)
    end
end

function mis_configurations(task, code; writefile=false)
    folder = joinpath("project", "data")
    if task == :missize
        nc = mis_size(code)
    elseif task == :miscounting
        nc = mis_count(code)
        println("MIS size = $(nc.n), degeneracy = $(nc.c)")
    elseif task == :misconfig
        nc = mis_config(code; all=false, usemask=true)[]
    elseif task == :allconfigs
        res = mis_config(code; all=true, usemask=true)[]
        return res
    elseif task == :independencepolynomial
        independence_polynomial(Val(:polynomial), code)
    elseif task == :numofis
        NoteOnTropicalMIS.mis_contract(1.0, code)[]
    else
        error("unknown task: $task")
    end
end

n = 200
rawcode = random_regular_eincode(n, 3);
code = optimize_kahypar(rawcode, uniformsize(rawcode, 2); sc_target=27, max_group_size=40);
tc, sc = OMEinsum.timespace_complexity(code, uniformsize(code,2))
println("time complexity = $tc, space complexity = $sc")
@time mis_configurations(:numofis, code)
@time mis_configurations(:missize, code)
@time mis_configurations(:miscounting, code)
@time mis_configurations(:misconfig, code)
@time mis_configurations(:allconfigs, code)
@time mis_configurations(:independencepolynomial, code)