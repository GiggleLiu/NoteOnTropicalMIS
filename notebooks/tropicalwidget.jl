### A Pluto.jl notebook ###
# v0.15.0

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 53ce616e-8055-11eb-0caf-e968c041ff4f
using Revise, Viznet, Compose, TropicalNumbers, OMEinsum

# ╔═╡ a4faf523-31f2-4446-82d2-059e39793e62
using PlutoUI

# ╔═╡ 448c3e6a-8203-11eb-2c9e-fd7d2a3db9d2
using Profile

# ╔═╡ 30b768c5-0ce0-4f4a-9981-a23f169b4a96
function vizconfig(nodes, edges, config=zeros(Int, length(nodes)); xyratio=1.0)
	X = 12cm
	Y = xyratio * 12cm
	Compose.set_default_graphic_size(X, Y)
	tb = textstyle(:default, fill("white"))
	nb = nodestyle(:default)
	nb2 = nodestyle(:default, fill("red"))
	eb = bondstyle(:default)
	img = canvas() do
		for (i, (t, p)) in enumerate(nodes)
			(config[i]==1 ? nb2 : nb) >> p
			tb >> (p, t)
		end
		for (i,j) in edges
			eb >> (nodes[i].second, nodes[j].second)
		end
	end
	Compose.compose(context(0, (1-X/Y)/2, 1.0, X/Y), img)
end

# ╔═╡ 100e3eb1-9d5c-4ae0-9a1b-d279927443fc
md"## We have a unit disk gadget"

# ╔═╡ 6c9ed9a9-4579-4260-adba-6af2df243ec1
# find edges in a unit disk graph, returns a vector of tuples
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

# ╔═╡ efdc5d3e-58e0-410d-b15a-2363e9c9f156
nodes = let
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

# ╔═╡ a7dec0bc-f9a2-43e1-9ec1-3fe0e38de8e6
md"""
*Note*
* `=>` defines a [pair](https://docs.julialang.org/en/v1/base/collections/)
"""

# ╔═╡ cefb1c82-7ef3-4f43-b5b0-c42eed297d8a
edges = find_edges(getindex.(nodes, 2), 0.23)

# ╔═╡ 6a766790-d4ff-4d4d-9a7d-4810963a52a1
md"""
*Note*
* dot means [broadcasting](https://docs.julialang.org/en/v1/manual/arrays/#Broadcasting),
* `getindex(pair, 2)` means getting the second field of a pair, here the location.
"""

# ╔═╡ d6a305f3-ddf7-445c-81f1-d9b9a3ec40bb
@assert length(edges) == 28

# ╔═╡ 78ee8772-83e4-40fb-8151-0123370481d9
vizconfig(nodes, edges, rand(Bool, 12); xyratio=0.5)

# ╔═╡ 56f3dcd1-11f1-4de1-bd73-356cc332d38d
md"## We convert it to a generic einsum network and contract it"

# ╔═╡ 07665e9c-6d2c-42e7-a577-2de1c3216883
md"""
We use Tropical numbers and its variants for computing.
Tropical number map `+` to `max`, `*` to `+`, its zero and one elements are changed correspondingly.

Check [this paper](https://arxiv.org/abs/2008.06888) and this [note](https://github.com/Happy-Diode/NoteOnTropicalMIS/blob/master/paper/paper.pdf) for details.
"""

# ╔═╡ bb36b716-8357-11eb-1007-536e94c427c9
zero(TropicalF64)

# ╔═╡ bfda24e2-8357-11eb-1fb9-7390cc2b0988
one(TropicalF64)

# ╔═╡ 82b1c65b-244d-4e38-b17c-7887cd07e19f
md"""
We construct a einsum network (a generalization to tensor network to hypergraph) by
* placing a bond tensor of rank 2 at an edge,
* placing a vertex tensor of rank one at a vertex.
"""

# ╔═╡ 69180ced-3729-4cf7-9fb4-db89de76469e
function bondtensor(::Type{T}) where T
	res = ones(T, 2, 2)
	res[2, 2] = zero(T)
	return res
end

# ╔═╡ 56ec7baa-57bc-4eae-b965-b997d959d24e
md"*Note*
* The input argument is a type, `where` statement is used to get the input argument type. Note the type of a type `T` is `Type{T}`. Here, the input argument is a type, and the `where` statement infers the type. We write a function like this rather than using `if-else` is because Julia compiles you code in a just in time way and type is a **static** information.
"

# ╔═╡ 5df8b62e-8357-11eb-2322-f54c02f1dc6c
bondtensor(TropicalF64)

# ╔═╡ f005035b-b982-47f8-b292-28aac5b5c00d
function vertextensor(val::T) where T
	[one(T), val]
end

# ╔═╡ 6c50f6de-8357-11eb-1a29-4f742afaeaa2
vertextensor(TropicalF64(1))

# ╔═╡ efed2d39-25e3-4b15-a184-c33789949aa6
res_gadget = ein"a,b,c,d,i,j,k,l,m,n,o,p,ai,aj,ak,al,bi,bm,bm,cn,co,cp,dl,dp,ij,
im,in,jk,jm,jn,jo,kl,kn,ko,kp,lo,lp,mn,no,op->abcd"(
	[vertextensor(TropicalF64(1)) for i=1:12]...,
	[bondtensor(TropicalF64) for i=1:28]...)

# ╔═╡ 0586bb1a-70ea-4923-bbfe-5373e1f45f9d
md"""
*Notes*

* `ein"a,b,c,d,i,j,k,l,m,n,o,p,ai,aj,ak,al,bi,bm,bm,cn,co,cp,dl,dp,ij,im,in,jk,jm,jn,jo,kl,kn,ko,kp,lo,lp,mn,no,op->abcd"` defines an einsum contraction pattern, meaning enumerating all symbols and accumulate the results to the output tensor. It is a callable, with tensors as inputs.
"""

# ╔═╡ 90dfe54c-349c-4302-a8f1-cf69e66fdbbc
md"Do it in a programming style"

# ╔═╡ 74661ef2-e74b-4b90-86b9-2d946a4e097f
function contract(edges, x::Vector{T}) where T
	n = length(x)
	code = EinCode(([(i,) for i=1:n]..., edges...), (1,2,3,4))
	code([vertextensor(xi) for xi in x]..., [bondtensor(T) for i=1:length(edges)]...)
end

# ╔═╡ 5ce13497-f034-461f-8c60-de66a2989722
md"""
*Note*

The `EinCode` constructor takes two arguments. The first is a tuple of input argument labels (tuple of tuple). The second argument is the output labels (tuple).
"""

# ╔═╡ 90c88c7f-e0c7-4845-a8b9-7f7562bd256f
contract(edges, [Tropical(1.0) for i=1:12])

# ╔═╡ f657f321-255e-44f2-a1d5-9093fa8eca28
md"## Obtaining configurations"

# ╔═╡ 0b184beb-33cf-488f-b40c-ad284c8cf993
md"Obtaining just one optimal configuration"

# ╔═╡ 8081596a-fa8a-47d3-b239-5eca34d34063
res_config = ein"a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,
ae,bf,cg,dh,ei,ej,ek,el,fi,fm,gm,gn,go,gp,hl,hp,ij,
im,in,jk,jm,jn,jo,kl,kn,ko,kp,lo,lp,mn,no,op->abcd"(
	[vertextensor(ConfigTropical(1.0, TropicalNumbers.onehot(StaticBitVector{16,1},i))) for i=1:16]...,
	[bondtensor(ConfigTropical{Float64,16,1}) for i=1:32]...)

# ╔═╡ 859125d9-d331-433a-9f3e-c72a21191c3c
md"""
`StaticBitVector{N,C}` defines a bit string. It has two type parameters. The first is the number of vertices. The second is the number of integers (UInt64) to store the configuration, it is computed as `C = (N-1) ÷ 64 + 1`. `TropicalNumbers.onehot` creates a one hot bit string. Its first and second argument are type and location.
"""

# ╔═╡ 8b300219-eca2-4db3-90ba-dc108ca46f04
bs1 = TropicalNumbers.onehot(StaticBitVector{16,1},2)

# ╔═╡ c70f064b-31bc-4c73-8d49-81f0bb801424
bs2 = TropicalNumbers.onehot(StaticBitVector{16,1},4)

# ╔═╡ 8bb43564-c87f-4a32-92b1-8f3029b46069
md"""
when we multiply two tropical numbers, we do the following operation to bitstring.
"""

# ╔═╡ def9abff-7e12-4022-8582-b361fa878961
bs1 | bs2

# ╔═╡ eab19833-5767-48ed-b120-9f2b85e0ad6e
md"""
To obtaining all optimal configurations, we redefine the algebra of the counting field of a counting tropical number.
"""

# ╔═╡ 3857ee49-5877-4cfa-b137-847cdeca0614
begin
struct ConfigEnumerator{N,C}
    data::Vector{StaticBitVector{N,C}}
end

Base.length(x::ConfigEnumerator{N}) where N = length(x.data)
Base.:(==)(x::ConfigEnumerator{N,C}, y::ConfigEnumerator{N,C}) where {N,C} = x.data == y.data

function Base.:+(x::ConfigEnumerator{N,C}, y::ConfigEnumerator{N,C}) where {N,C}
    length(x) == 0 && return y
    length(y) == 0 && return x
    return ConfigEnumerator{N,C}(vcat(x.data, y.data))
end

function Base.:*(x::ConfigEnumerator{L,C}, y::ConfigEnumerator{L,C}) where {L,C}
    M, N = length(x), length(y)
    M == 0 && return x
    N == 0 && return y
    z = Vector{StaticBitVector{L,C}}(undef, M*N)
    @inbounds for j=1:N, i=1:M
        z[(j-1)*M+i] = x.data[i] | y.data[j]
    end
    return ConfigEnumerator{L,C}(z)
end

Base.zero(::Type{ConfigEnumerator{N,C}}) where {N,C} = ConfigEnumerator{N,C}(StaticBitVector{N,C}[])
Base.one(::Type{ConfigEnumerator{N,C}}) where {N,C} = ConfigEnumerator{N,C}([TropicalNumbers.staticfalses(StaticBitVector{N,C})])
end

# ╔═╡ 4d93ee8e-38e4-4556-b8a1-21990ca12925
res_configs = contract(edges, [CountingTropical(1.0, ConfigEnumerator([TropicalNumbers.onehot(StaticBitVector{16,1},i)])) for i=1:16])

# ╔═╡ be48b33b-90d6-4fa5-a1ba-9b2e52e52d29
res_configs[1,1,1,1].c

# ╔═╡ caa3bc03-f6ff-4ebe-9608-518b8924adae
vizconfig(nodes,edges,res_configs[1,1,1,1].c.data[1]; xyratio=0.5)

# ╔═╡ 87f32eb9-9829-4d98-a666-de91eb4d21d1
vizconfig(nodes,edges,res_configs[1,1,1,1].c.data[2]; xyratio=0.5)

# ╔═╡ 7da26aae-83a6-11eb-35ea-b3b0a276855a
md"## Prooving the equivalence"

# ╔═╡ db52bcdc-1a89-4e24-898f-d3fdbbed27ea
res1 = let
	T = TropicalF64
	ein"a,b,c,d,ac,bd->abcd"([vertextensor(T(1)) for i=1:4]..., [bondtensor(T) for i=1:2]...)
end

# ╔═╡ c18e8111-16a8-42a8-8fa8-67ac76f4a13e
begin
	⪰(a, b) = b ⪯ a
	⪯(a, b) = (a & b) == a
end

# ╔═╡ 0086090b-a61f-42c5-b210-c4d7d53e8e25
function compress!(a::AbstractArray{T}) where T
	for (ind_a, val_a) in enumerate(a)
		for (ind_b, val_b) in enumerate(a)
			bs_a = ind_a - 1
			bs_b = ind_b - 1
			if bs_a != bs_b && val_a <= val_b && bs_a ⪰ bs_b
				a[ind_a] = zero(T)
			end
		end
	end
	return a
end

# ╔═╡ 201acfaa-8358-11eb-094b-e7241fd58324
content.(res_gadget) .- 2.0

# ╔═╡ b98adcc4-7b23-4d1e-bcf2-626ea5f4cc47
res3 = content.(compress!(copy(res_gadget))) .- 2

# ╔═╡ e278731d-3a77-460a-8ea2-d94ddb0851ac
md"**They must be different by a constant!**"

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

# ╔═╡ 87f380eb-90b5-4846-b1b1-f4b75c61a1a0
md"# Search"

# ╔═╡ c05885f2-3465-4cf3-b8f2-e2269044c63f
function trytry(locs, distance)
	edges = find_edges(locs, distance)
	y = contract(edges, Tropical.(ones(length(locs))))
	b, diff = is_diff_by_const(content.(res1), content.(compress!(y)))
	if b
		return true, diff
	else
		return false, 0
	end
end

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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Compose = "a81c6b42-2e10-5240-aca2-a61377ecd94b"
OMEinsum = "ebe7aa44-baf0-506c-a96f-8464559b3922"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Profile = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"
Revise = "295af30f-e4ad-537b-8983-00126c2a3abe"
TropicalNumbers = "b3a74e9c-7526-4576-a4eb-79c0d4c32334"
Viznet = "52a3aca4-6234-47fd-b74a-806bdf78ede9"

[compat]
Compose = "~0.9.2"
OMEinsum = "~0.4.5"
PlutoUI = "~0.7.9"
Revise = "~3.1.17"
TropicalNumbers = "~0.4.1"
Viznet = "~0.3.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "485ee0867925449198280d4af84bdb46a2a404d0"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.0.1"

[[AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[BFloat16s]]
deps = ["LinearAlgebra", "Test"]
git-tree-sha1 = "4af69e205efc343068dc8722b8dfec1ade89254a"
uuid = "ab4f0b2a-ad5b-11e8-123f-65d77653426b"
version = "0.1.0"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BatchedRoutines]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "8ee75390ba4bbfaf9aa48c121857b0da9a914265"
uuid = "a9ab73d0-e05c-5df1-8fde-d6a4645b8d8e"
version = "0.2.1"

[[BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Statistics", "UUIDs"]
git-tree-sha1 = "ffabdf5297c9038973a0a3724132aa269f38c448"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.1.0"

[[CEnum]]
git-tree-sha1 = "215a9aa4a1f23fbd05b92769fdd62559488d70e9"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.1"

[[CUDA]]
deps = ["AbstractFFTs", "Adapt", "BFloat16s", "CEnum", "CompilerSupportLibraries_jll", "DataStructures", "ExprTools", "GPUArrays", "GPUCompiler", "LLVM", "LazyArtifacts", "Libdl", "LinearAlgebra", "Logging", "Printf", "Random", "Random123", "RandomNumbers", "Reexport", "Requires", "SparseArrays", "SpecialFunctions", "TimerOutputs"]
git-tree-sha1 = "8ef71bf6d6602cf227196b43650924bf9ef7babc"
uuid = "052768ef-5323-5732-b1bb-66c8b64840ba"
version = "3.3.3"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "4b28f88cecf5d9a07c85b9ce5209a361ecaff34a"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "0.9.45"

[[CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "8ad457cfeb0bca98732c97958ef81000a543e73e"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.0.5"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "dc7dedc2c2aa9faf59a55c622760a25cbefbe941"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.31.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[Compose]]
deps = ["Base64", "Colors", "DataStructures", "Dates", "IterTools", "JSON", "LinearAlgebra", "Measures", "Printf", "Random", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "c6461fc7c35a4bb8d00905df7adafcff1fe3a6bc"
uuid = "a81c6b42-2e10-5240-aca2-a61377ecd94b"
version = "0.9.2"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4437b64df1e0adccc3e5d1adbc3ac741095e4677"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.9"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Dierckx]]
deps = ["Dierckx_jll"]
git-tree-sha1 = "5fefbe52e9a6e55b8f87cb89352d469bd3a3a090"
uuid = "39dd38d3-220a-591b-8e3c-4c3a8c710a94"
version = "0.5.1"

[[Dierckx_jll]]
deps = ["CompilerSupportLibraries_jll", "Libdl", "Pkg"]
git-tree-sha1 = "a580560f526f6fc6973e8bad2b036514a4e3b013"
uuid = "cd4c43a9-7502-52ba-aa6d-59fb2a88580b"
version = "0.0.1+0"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "a32185f5428d3986f47c2ab78b1f216d5e6cc96f"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.5"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[ExprTools]]
git-tree-sha1 = "b7e3d17636b348f005f11040025ae8c6f645fe92"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.6"

[[FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[GPUArrays]]
deps = ["AbstractFFTs", "Adapt", "LinearAlgebra", "Printf", "Random", "Serialization", "Statistics"]
git-tree-sha1 = "ececbf05f8904c92814bdbd0aafd5540b0bf2e9a"
uuid = "0c68f7d7-f131-5f86-a1c3-88cf8149b2d7"
version = "7.0.1"

[[GPUCompiler]]
deps = ["DataStructures", "ExprTools", "InteractiveUtils", "LLVM", "Libdl", "Logging", "TimerOutputs", "UUIDs"]
git-tree-sha1 = "e8a09182a4440489e2e3dedff5ad3f6bbe555396"
uuid = "61eb1bfa-7361-4325-ad38-22787b887f55"
version = "0.12.5"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "81690084b6198a2e1da36fcfda16eeca9f9f24e4"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.1"

[[JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "31c2eee64c1eee6e8e3f30d5a03d4b5b7086ab29"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.8.18"

[[LLVM]]
deps = ["CEnum", "LLVMExtra_jll", "Libdl", "Printf", "Unicode"]
git-tree-sha1 = "1b7ba36ea7aa6fa2278118951bad114fbb8359f2"
uuid = "929cbde3-209d-540e-8aea-75f648917ca0"
version = "4.1.0"

[[LLVMExtra_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b36c0677a0549c7d1dc8719899a4133abbfacf7d"
uuid = "dad2f222-ce93-54a1-a47d-0025e8a3acab"
version = "0.0.6+0"

[[LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[LeftChildRightSiblingTrees]]
deps = ["AbstractTrees"]
git-tree-sha1 = "71be1eb5ad19cb4f61fa8c73395c0338fd092ae0"
uuid = "1d6d02ad-be62-4b6b-8a6d-2f90e265016e"
version = "0.1.2"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["DocStringExtensions", "LinearAlgebra"]
git-tree-sha1 = "7bd5f6565d80b6bf753738d2bc40a5dfea072070"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.2.5"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "4bfb8b57df913f3b28a6bd3bdbebe9a50538e689"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.1.0"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "6a8a2a625ab0dea913aba95c11370589e0239ff0"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.6"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[OMEinsum]]
deps = ["AbstractTrees", "BatchedRoutines", "CUDA", "ChainRulesCore", "Combinatorics", "LinearAlgebra", "MacroTools", "PkgBenchmark", "Requires", "Test", "TupleTools"]
git-tree-sha1 = "135a0e57deed5fb80122375488194876c842ecc1"
uuid = "ebe7aa44-baf0-506c-a96f-8464559b3922"
version = "0.4.5"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "c8abc88faa3f7a3950832ac5d6e690881590d6dc"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "1.1.0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PkgBenchmark]]
deps = ["BenchmarkTools", "Dates", "InteractiveUtils", "JSON", "LibGit2", "Logging", "Pkg", "Printf", "TerminalLoggers", "UUIDs"]
git-tree-sha1 = "e4a10b7cdb7ec836850e43a4cee196f4e7b02756"
uuid = "32113eaa-f34f-5b0d-bd6c-c81e245fc73d"
version = "0.2.12"

[[PlutoUI]]
deps = ["Base64", "Dates", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "Suppressor"]
git-tree-sha1 = "44e225d5837e2a2345e69a1d1e01ac2443ff9fcb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.9"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[ProgressLogging]]
deps = ["Logging", "SHA", "UUIDs"]
git-tree-sha1 = "80d919dee55b9c50e8d9e2da5eeafff3fe58b539"
uuid = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
version = "0.1.4"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Random123]]
deps = ["Libdl", "Random", "RandomNumbers"]
git-tree-sha1 = "0e8b146557ad1c6deb1367655e052276690e71a3"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.4.2"

[[RandomNumbers]]
deps = ["Random", "Requires"]
git-tree-sha1 = "441e6fc35597524ada7f85e13df1f4e10137d16f"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.4.0"

[[Reexport]]
git-tree-sha1 = "5f6c21241f0f655da3952fd60aa18477cf96c220"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.1.0"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "410bbe13d9a7816e862ed72ac119bda7fb988c08"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.1.17"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "LogExpFunctions", "OpenSpecFun_jll"]
git-tree-sha1 = "a50550fa3164a8c46747e62063b4d774ac1bcf49"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.5.1"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[TerminalLoggers]]
deps = ["LeftChildRightSiblingTrees", "Logging", "Markdown", "Printf", "ProgressLogging", "UUIDs"]
git-tree-sha1 = "d620a061cb2a56930b52bdf5cf908a5c4fa8e76a"
uuid = "5d786b92-1e48-4d6f-9151-6b4477ca9bed"
version = "0.1.4"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "209a8326c4f955e2442c07b56029e88bb48299c7"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.12"

[[TropicalNumbers]]
git-tree-sha1 = "ed96a2230559fcb69418d524efc821e524148695"
uuid = "b3a74e9c-7526-4576-a4eb-79c0d4c32334"
version = "0.4.1"

[[TupleTools]]
git-tree-sha1 = "3c712976c47707ff893cf6ba4354aa14db1d8938"
uuid = "9d95972d-f1c8-5527-a6e0-b4b365fa01f6"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Viznet]]
deps = ["Compose", "Dierckx"]
git-tree-sha1 = "7a022ae6ac8b153d47617ed8c196ce60645689f1"
uuid = "52a3aca4-6234-47fd-b74a-806bdf78ede9"
version = "0.3.3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═53ce616e-8055-11eb-0caf-e968c041ff4f
# ╠═30b768c5-0ce0-4f4a-9981-a23f169b4a96
# ╟─100e3eb1-9d5c-4ae0-9a1b-d279927443fc
# ╠═6c9ed9a9-4579-4260-adba-6af2df243ec1
# ╠═efdc5d3e-58e0-410d-b15a-2363e9c9f156
# ╟─a7dec0bc-f9a2-43e1-9ec1-3fe0e38de8e6
# ╠═cefb1c82-7ef3-4f43-b5b0-c42eed297d8a
# ╟─6a766790-d4ff-4d4d-9a7d-4810963a52a1
# ╠═d6a305f3-ddf7-445c-81f1-d9b9a3ec40bb
# ╠═78ee8772-83e4-40fb-8151-0123370481d9
# ╟─56f3dcd1-11f1-4de1-bd73-356cc332d38d
# ╟─07665e9c-6d2c-42e7-a577-2de1c3216883
# ╠═bb36b716-8357-11eb-1007-536e94c427c9
# ╠═bfda24e2-8357-11eb-1fb9-7390cc2b0988
# ╟─82b1c65b-244d-4e38-b17c-7887cd07e19f
# ╠═69180ced-3729-4cf7-9fb4-db89de76469e
# ╟─56ec7baa-57bc-4eae-b965-b997d959d24e
# ╠═5df8b62e-8357-11eb-2322-f54c02f1dc6c
# ╠═f005035b-b982-47f8-b292-28aac5b5c00d
# ╠═6c50f6de-8357-11eb-1a29-4f742afaeaa2
# ╠═efed2d39-25e3-4b15-a184-c33789949aa6
# ╟─0586bb1a-70ea-4923-bbfe-5373e1f45f9d
# ╟─90dfe54c-349c-4302-a8f1-cf69e66fdbbc
# ╠═74661ef2-e74b-4b90-86b9-2d946a4e097f
# ╟─5ce13497-f034-461f-8c60-de66a2989722
# ╠═90c88c7f-e0c7-4845-a8b9-7f7562bd256f
# ╟─f657f321-255e-44f2-a1d5-9093fa8eca28
# ╟─0b184beb-33cf-488f-b40c-ad284c8cf993
# ╠═8081596a-fa8a-47d3-b239-5eca34d34063
# ╟─859125d9-d331-433a-9f3e-c72a21191c3c
# ╠═8b300219-eca2-4db3-90ba-dc108ca46f04
# ╠═c70f064b-31bc-4c73-8d49-81f0bb801424
# ╟─8bb43564-c87f-4a32-92b1-8f3029b46069
# ╠═def9abff-7e12-4022-8582-b361fa878961
# ╟─eab19833-5767-48ed-b120-9f2b85e0ad6e
# ╠═3857ee49-5877-4cfa-b137-847cdeca0614
# ╠═4d93ee8e-38e4-4556-b8a1-21990ca12925
# ╠═be48b33b-90d6-4fa5-a1ba-9b2e52e52d29
# ╠═caa3bc03-f6ff-4ebe-9608-518b8924adae
# ╠═87f32eb9-9829-4d98-a666-de91eb4d21d1
# ╟─7da26aae-83a6-11eb-35ea-b3b0a276855a
# ╠═db52bcdc-1a89-4e24-898f-d3fdbbed27ea
# ╠═c18e8111-16a8-42a8-8fa8-67ac76f4a13e
# ╠═0086090b-a61f-42c5-b210-c4d7d53e8e25
# ╠═201acfaa-8358-11eb-094b-e7241fd58324
# ╠═b98adcc4-7b23-4d1e-bcf2-626ea5f4cc47
# ╟─e278731d-3a77-460a-8ea2-d94ddb0851ac
# ╠═e8b1ffa5-ea35-4a5e-8787-a8aa45574db5
# ╠═91e8b86c-9040-482a-b4f1-1419e720799c
# ╟─87f380eb-90b5-4846-b1b1-f4b75c61a1a0
# ╠═c05885f2-3465-4cf3-b8f2-e2269044c63f
# ╠═e831bc1f-0d92-44ce-9a0c-6b7c5a4d0f68
# ╠═a4faf523-31f2-4446-82d2-059e39793e62
# ╠═4da55afd-cdd8-49e1-bdf4-662b212f9153
# ╠═1ddc529f-d1cb-4bcc-876f-5501e4b8ef77
# ╠═448c3e6a-8203-11eb-2c9e-fd7d2a3db9d2
# ╠═86589fd8-19a6-4236-af20-c0ae395c2418
# ╠═57115702-ab8c-45d4-afd5-18d4fdcf18bd
# ╠═87adb73e-4514-49a9-9196-0d02f06702f1
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
