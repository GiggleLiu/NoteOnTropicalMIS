using LightGraphs, OMEinsumContractionOrders
export random_regular_eincode

function random_regular_eincode(n, k; optimize=nothing)
	g = LightGraphs.random_regular_graph(n, k)
	ixs = [minmax(e.src,e.dst) for e in LightGraphs.edges(g)]
	return EinCode((ixs..., [(i,) for i in LightGraphs.vertices(g)]...), ())
end