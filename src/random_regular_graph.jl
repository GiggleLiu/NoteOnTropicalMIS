using LightGraphs
export random_regular_eincode

function random_regular_eincode(n, k)
	g = LightGraphs.random_regular_graph(n, k)
	ixs = [minmax(e.src,e.dst) for e in LightGraphs.edges(g)]
	return EinCode((ixs..., [(i,) for i in LightGraphs.vertices(g)]...), ())
end

function solve_regular_graph(n::Int, k::Int; sc_target=20)
    code = random_regular_eincode(n, k)
    optcode = optimize_kahypar(code, uniformsize(code, 2); sc_target=sc_target, max_group_size=40)
    mis_size = mis_solve(optcode)[].n |> Int

	xs = (0:mis_size)
	ys = [mis_polysolve(optcode, x)[] for x in xs]
	fit(xs, ys, mis_size)
end