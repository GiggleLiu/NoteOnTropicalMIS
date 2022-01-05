using Distributed, CUDA, GraphTensorNetworks, Random, DelimitedFiles, GraphTensorNetworks.OMEinsumContractionOrders, GraphTensorNetworks.OMEinsum
USECUDA = parse(Bool, ARGS[1])
@show USECUDA

function countmax(graph, n; writefile, sc_target, usecuda, nslices)
    mask = trues(n, n)
    if graph == "diag"
        g = diagonal_coupled_graph(mask)
    elseif graph == "square"
        g = square_lattice_graph(mask)
    end
    gp = Independence(g; optimizer=TreeSA(sc_target=sc_target, sc_weight=1.0, nslices=nslices;
        ntrials=4, βs=0.01:0.05:25.0, niters=50, rw_weight=1.0), simplifier=MergeGreedy())
    println("Graph size $n, usecuda = $usecuda")
    @show timespace_complexity(gp)
    res = GraphTensorNetworks.big_integer_solve(Int32, 100) do T
        @info "T = $T"
        filename = joinpath(@__DIR__, "data", graph, "maxcount-$n-$(GraphTensorNetworks.modulus(one(T))).dat")
        if isfile(filename)
            fill(T(readdlm(filename)[]))
        else
            xs = GraphTensorNetworks.generate_tensors(x->CountingTropical(Int32(1), one(T)), gp)
            if usecuda
                xs = CUDA.CuArray.(xs)
            end
            @time res = Array(gp.code(xs...))
            writedlm(filename, res[].c.val)
            fill(res[].c)
        end
    end
    @show res
    ofname = joinpath(@__DIR__, "data", graph, "maxcount-$n.dat")
    writefile && writedlm(ofname, res)
end

# current best for 38 = 204: 49.27
# current best for 39 = 7: 49.27
for L=28
    Random.seed!(2)
    graph = "diag"
    println("computing graph $graph, L = $L")
    @time countmax(graph, L; writefile=true, sc_target=28, usecuda=USECUDA, nslices=(L-26)*2-1)
end
