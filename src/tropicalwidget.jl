using Revise, Viznet, Compose, TropicalNumbers, OMEinsum
using ForwardDiff

# ╔═╡ 4053c327-43d2-4203-aa60-48b02a8c4dac
nodes = let
	a = 0.12
	ymid = xmid = 0.5
	X = 0.45
	Y = 0.33
	D = 0.15
	y = [ymid-Y, ymid-Y+D, ymid-a/2, ymid+a/2, ymid+Y-D, ymid+Y]
	x = [xmid-X, xmid-X+D, xmid-1.5a, xmid-a/2, xmid+a/2, xmid+1.5a, xmid+X-D, xmid+X]
	xmin, xmax, ymin, ymax = x[1], x[end], y[1], y[end]
	["a"=>(xmid, y[1]), "b"=>(xmin, ymid), "c"=>(xmid, ymax), "d"=>(xmax, ymid),
		"e"=>(xmid, y[2]), "f"=>(x[2], ymid), "g"=>(xmid, y[end-1]),
		"h"=>(x[end-1], ymid), "i"=>(x[3], y[3]), "j"=>(x[4], y[3]),
		"k"=>(x[5], y[3]), "l"=>(x[6], y[3]), "m"=>(x[3], y[4]),
		"n"=>(x[4], y[4]), "o"=>(x[5], y[4]), "p"=>(x[6], y[4])]
end

# ╔═╡ 6c9ed9a9-4579-4260-adba-6af2df243ec1
function find_edges(nodes, distance)
	edges = Tuple{Int,Int}[]
	for (i, p) in enumerate(nodes)
		for (j,p2) in enumerate(nodes)
			if i<j && sqrt(sum(abs2, p2 .- p)) < distance
				push!(edges, (i,j))
			end
		end
	end
	edges
end

# ╔═╡ ef28826a-24a5-48c7-b819-e42be90d05ab
edges = find_edges(getindex.(nodes, 2), sqrt(0.05))

# ╔═╡ e35ec714-70be-47bf-982f-1ef20d7ed470
@assert length(edges) == 32

# ╔═╡ 30b768c5-0ce0-4f4a-9981-a23f169b4a96
function vizconfig(nodes, edges, config=zeros(Int, length(nodes)))
	Compose.set_default_graphic_size(12cm, 12cm)
	tb = textstyle(:default, fill("white"))
	nb = nodestyle(:default)
	nb2 = nodestyle(:default, fill("red"))
	eb = bondstyle(:default)
	canvas() do
		for (i, (t, p)) in enumerate(nodes)
			(config[i]==1 ? nb2 : nb) >> p
			tb >> (p, t)
		end
		for (i,j) in edges
			eb >> (nodes[i].second, nodes[j].second)
		end
	end
end

# ╔═╡ b57d89c7-b46a-474b-a701-b8632338a01b
vizconfig(nodes, edges)

# ╔═╡ efed2d39-25e3-4b15-a184-c33789949aa6
res2 = contract_a(TropicalF64, ones(16))

# ╔═╡ 7da26aae-83a6-11eb-35ea-b3b0a276855a
md"## Prooving the equivalence"

# ╔═╡ 201acfaa-8358-11eb-094b-e7241fd58324
content.(res2) .- 4.0

# ╔═╡ b98adcc4-7b23-4d1e-bcf2-626ea5f4cc47
res3 = content.(compress!(copy(res2))) .- 4

# ╔═╡ e278731d-3a77-460a-8ea2-d94ddb0851ac
md"**They must be different by a constant!**"

# ╔═╡ 113aaf3f-7471-442a-a5f6-8eca78b44030
config = let
	vals = ones(16)
	ForwardDiff.gradient(x->contract_a(Tropical{eltype(x)}, x)[2].n, vals)
end

# ╔═╡ 3a568b30-6407-44c7-8a19-1b2d03a44179
vizconfig(nodes, edges, config)

# ╔═╡ df6060bd-b4b9-4209-b610-f3d27945a6db
md"# Simpler graph"

# ╔═╡ efdc5d3e-58e0-410d-b15a-2363e9c9f156
nodes_simple = let
	a = 0.12
	ymid = xmid = 0.5
	X = 0.33
	Y = 0.17
	D = 0.15
	y = [ymid-Y, ymid-Y+D, ymid-a/2, ymid+a/2, ymid+Y-D, ymid+Y]
	x = [xmid-X, xmid-X+D, xmid-1.5a, xmid-a/2, xmid+a/2, xmid+1.5a, xmid+X-D, xmid+X]
	xmin, xmax, ymin, ymax = x[1], x[end], y[1], y[end]
	["a"=>(xmid, y[1]), "b"=>(xmin, ymid), "c"=>(xmid, ymax), "d"=>(xmax, ymid),
		"i"=>(x[3], y[3]), "j"=>(x[4], y[3]),
		"k"=>(x[5], y[3]), "l"=>(x[6], y[3]), "m"=>(x[3], y[4]),
		"n"=>(x[4], y[4]), "o"=>(x[5], y[4]), "p"=>(x[6], y[4])]
end

# ╔═╡ cefb1c82-7ef3-4f43-b5b0-c42eed297d8a
edges_simple = find_edges(getindex.(nodes_simple, 2), 0.23)

# ╔═╡ d6a305f3-ddf7-445c-81f1-d9b9a3ec40bb
@assert length(edges_simple) == 28

# ╔═╡ 78ee8772-83e4-40fb-8151-0123370481d9
vizconfig(nodes_simple, edges_simple)

# ╔═╡ 2f1ac010-83ac-11eb-2cde-13653ed10cbc
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

# ╔═╡ 39cdf52c-83ac-11eb-26f5-157131b9d365
vizconfig(nodes_2, [1=>3, 2=>6, 4=>5, 7=>9, 8=>10, 11=>13, 12=>14, 10=>15, 12=>16, 15=>16])

# ╔═╡ 2d8bdb56-83ac-11eb-2e38-833260805680
nodes_simple2 = let
	offsetx = 0.4
	set1 = map(x->(""=>x.second .* 0.5), nodes_simple)
	set2 = map(x->(""=>x.second .* 0.5 .+ (offsetx, 0.0)), nodes_simple)
	[set1..., set2...]
end

# ╔═╡ 44ba9b02-83ac-11eb-3098-f7fa5f338177
vizconfig(nodes_simple2, find_edges(getindex.(nodes_simple2, 2), 0.115))

# ╔═╡ 17ea4dbc-88ce-11eb-3339-917291f3e932
1-0-0-1 => 1-0-1-0

# ╔═╡ 3b89fdf8-83a6-11eb-3761-2d9f346fb157
md"## connecting multiple crossing widgets"

# ╔═╡ 5a91ae70-83a9-11eb-3385-772978273c2a
# all widgets compressed locally -> no violation
# there are violations -> not all widgets are compressed.

# ╔═╡ 87f380eb-90b5-4846-b1b1-f4b75c61a1a0
md"# Search"

# ╔═╡ 91e8b86c-9040-482a-b4f1-1419e720799c
function is_diff_by_const(t1::AbstractArray{T}, t2::AbstractArray{T}) where T
	x = NaN
	for (a, b) in zip(t1, t2)
		if isinf(a) && isinf(b)
			continue
		end
		if isinf(a) || isinf(b)
			return false, 0
		end
		if isnan(x)
			x = (a - b)
		elseif x != a - b
			return false, 0
		end
	end
	return true, x
end

# ╔═╡ e8b1ffa5-ea35-4a5e-8787-a8aa45574db5
@assert is_diff_by_const(content.(res1), res3)[1]

# ╔═╡ b00671b1-916f-49b4-a8a5-ef6f0398ce05
function contract(edges, x)
	n = length(x)
	code = EinCode(([(i,) for i=1:n]..., edges...), (1,2,3,4))
	code([v(Tropical{eltype(x)}, 1, xi) for xi in x]..., [b(Tropical{eltype(x)}) for i=1:length(edges)]...)
end

# ╔═╡ 6e9b602c-847a-11eb-3faf-57f38c071f63
function contract_b(::Type{T}, x) where T
	edges = find_edges(getindex.(nodes_simple, 2), 0.23)
	contract(edges, x)
end

# ╔═╡ ba980320-847a-11eb-3da5-8d655b20c14f
contract_b(TropicalF64, ones(12))

# ╔═╡ 3a5acd0e-83a6-11eb-3f29-67de4e769f1e
function contract2(::Type{T}, x) where T
	w1 = contract_b(T, x[1:12])
	w2 = contract_b(T, x[13:24])
	fun = ein"abcd,df,efgh->abcegh"
	case0 = fun(w1, b(T), w2)
	# either a == 0 or c == 0, either g == 0 or e == 0
	blr = b(T)[1:1,1:1]
	# a = 0, g = 0
	case1 = fun(w1[1:1,:,:,1:1], blr, w2[:,1:1,1:1,:])
	# a = 0, e = 0
	case2 = fun(w1[1:1,:,:,1:1], blr, w2[1:1,1:1,:,:])
	# c = 0, e = 0
	case3 = fun(w1[:,:,1:1,1:1], blr, w2[1:1,1:1,:,:])
	# c = 0, g = 0
	case4 = fun(w1[:,:,1:1,1:1], blr, w2[:,1:1,1:1,:])
	case1# .+ case2 .+ case3 .+ case4
end

# ╔═╡ 032bdae8-83a7-11eb-0e8c-b782cf7c45a3
maximum(content.(contract2(Tropical{Float64}, ones(24))))

# ╔═╡ 8cf330e6-83a7-11eb-011c-ab5566255ca9
contract2(Tropical{Float64}, ones(24))[:,2:2,:,:,:,2:2]

# ╔═╡ 858340cc-8470-11eb-3ebf-15e5c5597929
configs2 = ForwardDiff.gradient(x->contract2(Tropical{eltype(x)}, x)[:,2:2,:,:,:,2:2][end,end].n, ones(32))

# ╔═╡ e6eef592-8476-11eb-2fdd-5b6d52730e66
vizconfig(nodes_simple2, find_edges(getindex.(nodes_simple2, 2), 0.115), configs2)

# ╔═╡ c05885f2-3465-4cf3-b8f2-e2269044c63f
function trytry(locs, distance)
	edges = find_edges(locs, distance)
	@show length(edges)
	y = contract(edges, ones(length(locs)))
	b, diff = is_diff_by_const(content.(res1), content.(compress!(y)))
	if b
		return true, diff
	else
		return false, 0
	end
end

# ╔═╡ ece0e9b8-2d33-4a24-8090-ea49eabdfcb8
trytry(getindex.(nodes_simple, 2), 0.23)

# ╔═╡ e831bc1f-0d92-44ce-9a0c-6b7c5a4d0f68
function run(nancs, distance)
	ymid = xmid = 0.5
	X = 0.45
	xmin, xmax = xmid-X, xmid+X
	Y = 0.45
	ymin, ymax = ymid-Y, ymid+Y
	local locs
	var = 0.6
	for i=1:1
		@show i
		loc_ancs = zip(rand(nancs) .* var .+ (0.5-var/2), rand(nancs) .* var .+ (0.5-var/2))
		locs = [(xmid, ymin), (xmin, ymid), (xmid, ymax), (xmax, ymid), loc_ancs...]
		res, diff = trytry(locs, distance)
		if res
			return true, locs, diff
		end
	end
	return false, locs, 0
end

# ╔═╡ 4da55afd-cdd8-49e1-bdf4-662b212f9153
@bind clock Clock(0.1)

# ╔═╡ 1ddc529f-d1cb-4bcc-876f-5501e4b8ef77
distance = 0.6

# ╔═╡ 5301304a-8203-11eb-246e-6799c32b6f5b
begin
	Profile.clear()
	@profile for i=1:3 run(6, distance) end
	Profile.print(mincount=100, format=:flat)
end

# ╔═╡ 86589fd8-19a6-4236-af20-c0ae395c2418
success, mylocs = let
	clock
	run(6, distance)
end

# ╔═╡ 57115702-ab8c-45d4-afd5-18d4fdcf18bd
if success
	sleep(100000)
end

# ╔═╡ 87adb73e-4514-49a9-9196-0d02f06702f1
vizconfig([""=>loc for loc in mylocs], find_edges(mylocs, distance))

# ╔═╡ Cell order:
# ╠═53ce616e-8055-11eb-0caf-e968c041ff4f
# ╠═69180ced-3729-4cf7-9fb4-db89de76469e
# ╠═5df8b62e-8357-11eb-2322-f54c02f1dc6c
# ╠═f005035b-b982-47f8-b292-28aac5b5c00d
# ╠═6c50f6de-8357-11eb-1a29-4f742afaeaa2
# ╠═bb36b716-8357-11eb-1007-536e94c427c9
# ╠═bfda24e2-8357-11eb-1fb9-7390cc2b0988
# ╠═c18e8111-16a8-42a8-8fa8-67ac76f4a13e
# ╠═0086090b-a61f-42c5-b210-c4d7d53e8e25
# ╠═db52bcdc-1a89-4e24-898f-d3fdbbed27ea
# ╠═4053c327-43d2-4203-aa60-48b02a8c4dac
# ╠═6c9ed9a9-4579-4260-adba-6af2df243ec1
# ╠═ef28826a-24a5-48c7-b819-e42be90d05ab
# ╠═e35ec714-70be-47bf-982f-1ef20d7ed470
# ╠═30b768c5-0ce0-4f4a-9981-a23f169b4a96
# ╠═b57d89c7-b46a-474b-a701-b8632338a01b
# ╠═d583d1ef-5f72-4331-95f5-2aa4693cdca9
# ╠═efed2d39-25e3-4b15-a184-c33789949aa6
# ╟─7da26aae-83a6-11eb-35ea-b3b0a276855a
# ╠═201acfaa-8358-11eb-094b-e7241fd58324
# ╠═b98adcc4-7b23-4d1e-bcf2-626ea5f4cc47
# ╟─e278731d-3a77-460a-8ea2-d94ddb0851ac
# ╠═e8b1ffa5-ea35-4a5e-8787-a8aa45574db5
# ╠═fadb61d5-0aac-496e-83e5-da89be0589ba
# ╠═113aaf3f-7471-442a-a5f6-8eca78b44030
# ╠═3a568b30-6407-44c7-8a19-1b2d03a44179
# ╟─df6060bd-b4b9-4209-b610-f3d27945a6db
# ╠═efdc5d3e-58e0-410d-b15a-2363e9c9f156
# ╠═cefb1c82-7ef3-4f43-b5b0-c42eed297d8a
# ╠═d6a305f3-ddf7-445c-81f1-d9b9a3ec40bb
# ╠═78ee8772-83e4-40fb-8151-0123370481d9
# ╠═ece0e9b8-2d33-4a24-8090-ea49eabdfcb8
# ╠═2f1ac010-83ac-11eb-2cde-13653ed10cbc
# ╠═39cdf52c-83ac-11eb-26f5-157131b9d365
# ╠═2d8bdb56-83ac-11eb-2e38-833260805680
# ╠═44ba9b02-83ac-11eb-3098-f7fa5f338177
# ╠═17ea4dbc-88ce-11eb-3339-917291f3e932
# ╟─3b89fdf8-83a6-11eb-3761-2d9f346fb157
# ╠═5a91ae70-83a9-11eb-3385-772978273c2a
# ╠═6e9b602c-847a-11eb-3faf-57f38c071f63
# ╠═ba980320-847a-11eb-3da5-8d655b20c14f
# ╠═3a5acd0e-83a6-11eb-3f29-67de4e769f1e
# ╠═032bdae8-83a7-11eb-0e8c-b782cf7c45a3
# ╠═8cf330e6-83a7-11eb-011c-ab5566255ca9
# ╠═858340cc-8470-11eb-3ebf-15e5c5597929
# ╠═e6eef592-8476-11eb-2fdd-5b6d52730e66
# ╟─87f380eb-90b5-4846-b1b1-f4b75c61a1a0
# ╠═91e8b86c-9040-482a-b4f1-1419e720799c
# ╠═b00671b1-916f-49b4-a8a5-ef6f0398ce05
# ╠═c05885f2-3465-4cf3-b8f2-e2269044c63f
# ╠═e831bc1f-0d92-44ce-9a0c-6b7c5a4d0f68
# ╠═a4faf523-31f2-4446-82d2-059e39793e62
# ╠═4da55afd-cdd8-49e1-bdf4-662b212f9153
# ╠═1ddc529f-d1cb-4bcc-876f-5501e4b8ef77
# ╠═448c3e6a-8203-11eb-2c9e-fd7d2a3db9d2
# ╠═5301304a-8203-11eb-246e-6799c32b6f5b
# ╠═86589fd8-19a6-4236-af20-c0ae395c2418
# ╠═57115702-ab8c-45d4-afd5-18d4fdcf18bd
# ╠═87adb73e-4514-49a9-9196-0d02f06702f1
