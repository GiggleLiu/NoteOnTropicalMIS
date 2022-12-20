# To run this script, you need to install required packages in a Julia REPL:
# ```julia
# using Pkg; Pkg.add(["GenericTensorNetworks", "CUDA", "Comonicon", "DelimitedFiles"])
# ```
#
# Then open a terminal, go to the file folder and type
# ```bash
# $ julia fib2d.jl 20
# ```
#
# Report an issue here: https://github.com/QuEraComputing/GenericTensorNetworks.jl/issues
# If you want the multi-GPU version of this code, please let me know: cacate0129@gmail.com

using GenericTensorNetworks, Random, Comonicon, DelimitedFiles
using CUDA; CUDA.allowscalar(false)

"""
Computing the two dimensional generalizatioon of Fibonacci numbers.

# Arguments
* `L`: the lattice size L x L.

# Options
* `--usecuda`: the CUDA device, -1 for not using GPU. Default is -1.
* `--nslices`: the number of slices, which is used to reduce the space complexity. Default is `L-28`.
* `--filename`: the filename to store the output. Default is `""` (not writting to a file).
* `--seed`: the seed for finding contraction order. Default is `2`.
"""
@main function fib2d(L::Int; nslices::Int=max(0, L-28), usecuda::Int=-1, filename::String="", seed::Int=2)
    usecuda >= 0 && CUDA.device!(usecuda)
    Random.seed!(seed)
    # the 2D Fibonacci number is defined by the total number of independent sets of a square lattice
    g = square_lattice_graph(trues(L, L))
    # both `optimizer` and `simplifier` are for tensor network contraction order optimization.
    gp = Independence(g; optimizer=TreeSA(sc_target=28, sc_weight=1.0, nslices=nslices;
        ntrials=7, βs=0.01:0.05:22.0, niters=20, rw_weight=2.0), simplifier=MergeGreedy())
    @info "The graph size $L x $L, usecuda = $usecuda"
    tc, sc, rw = timespacereadwrite_complexity(gp)
    # complexities are defiened by the log2 number of arithmetic operations/elements/read write operations.
    @info "Time complexity = $tc, space complexity = $sc, read-write complexity = $rw"
   
    # `big_integer_solve` uses the Chinese remainder theorem to combine many small remainders into a big integer.
    # Note the number of independent sets is a very big number that can not be stored with fixed width integers.
    res = GenericTensorNetworks.big_integer_solve(Int32, 100) do T
        @info "Computing on finite field algebra: $T"
        # `contractx` is for contracting the tensor network and evaluate the independence polynomial at `x = 1`,
        # which corresponds to the number of independent sets.
        # Here, `x` has the same value 1, but different finite field algebra, i.e. the modulus is different.
        @time Array(GenericTensorNetworks.contractx(gp, one(T); usecuda=usecuda>=0))
    end
    @info "result = $res"
    if !isempty(filename)
        @info "writing result to $filename"
        writedlm(filename, res)
    end
end
