content.(res2) .- 4.0
res3 = content.(compress!(copy(res2))) .- 4
config = let
	vals = ones(16)
	ForwardDiff.gradient(x->contract_a(Tropical{eltype(x)}, x)[2].n, vals)
end
vizconfig(nodes, edges, config)

edges_simple = find_edges(getindex.(nodes_simple, 2), 0.23)
@test length(edges_simple) == 28
vizconfig(nodes_simple, edges_simple)
nodes_2 = let
	a = 0.12
	ymid = xmid = 0.5
	X = 0.33
	Y = 0.17
	D = 0.15
	y = [ymid-Y, ymid-Y+D, ymid-a/2, ymid+a/2, ymid+Y-D, ymid+Y]
	x = [xmid-X, xmid-X+D, xmid-1.5a, xmid-a/2, xmid+a/2, xmid+1.5a, xmid+X-D, xmid+X]
	xmin, xmax, ymin, ymax = x[1], x[end], y[1], y[end]
	offsetx = 0.5
	set1 = [""=>(xmid*0.5, y[1]*0.5), ""=>(xmin*0.5, ymid*0.5), ""=>(xmid*0.5, ymax*0.5)]
	set2 = [""=>(xmid*0.5+offsetx, y[1]*0.5), ""=>(xmid*0.5+offsetx, ymax*0.5), ""=>(xmax*0.5+offsetx, ymid*0.5)]
	[set1..., set2...]
	
	offsety = 0.5
	set3 = [""=>(xmid*0.5, y[1]*0.5+offsety), ""=>(xmin*0.5, ymid*0.5+offsety), ""=>(xmid*0.5, ymax*0.5+offsety), ""=>(xmax*0.5, ymid*0.5+offsety)]
	set4 = [""=>(xmid*0.5+offsetx, y[1]*0.5+offsety), ""=>(xmin*0.5+offsetx, ymid*0.5+offsety), ""=>(xmid*0.5+offsetx, ymax*0.5+offsety), ""=>(xmax*0.5+offsetx, ymid*0.5+offsety)]
	set5 = [""=>(0.5, 0.7), ""=>(0.5, 0.8)]
	vcat(set1, set2, set3, set4, set5)
end

vizconfig(nodes_2, [1=>3, 2=>6, 4=>5, 7=>9, 8=>10, 11=>13, 12=>14, 10=>15, 12=>16, 15=>16])

nodes_simple2 = let
	offsetx = 0.4
	set1 = map(x->(""=>x.second .* 0.5), nodes_simple)
	set2 = map(x->(""=>x.second .* 0.5 .+ (offsetx, 0.0)), nodes_simple)
	[set1..., set2...]
end

vizconfig(nodes_simple2, find_edges(getindex.(nodes_simple2, 2), 0.115))

1-0-0-1 => 1-0-1-0