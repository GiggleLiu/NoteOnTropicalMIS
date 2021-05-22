### A Pluto.jl notebook ###
# v0.14.5

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

# ╔═╡ c8c14900-b7f5-11eb-3dbe-4390dd0fcc0a
using Revise, OMEinsum, TropicalNumbers, NoteOnTropicalMIS, Compose, Viznet

# ╔═╡ d1a8f940-c8b8-4ce6-af68-dee3a286e632
using LightGraphs

# ╔═╡ 8e2263d2-ffbf-4aa5-bf7e-29f16b2dbade
using PlutoUI

# ╔═╡ e56f2ba7-6a50-4400-b761-9f587076778f
using Polynomials

# ╔═╡ f2c1633a-0e2b-4ecc-81bf-f54fa0efcb44
using FFTW

# ╔═╡ 0a5712bd-4f8f-41a4-8fa3-7f947a078a4a
begin
	using Mods
	function chinese_remainder(remainders, primes)
		M = prod(primes)
		Ms = M .÷ primes
		ti = invmod.(Ms, primes)
		mod(sum(remainders .* ti .* Ms), M)
	end
	devisor(::Mod{N}) where N = N
	value(r::Mod{N}) where N = r.val
	function chinese_remainder(rs::Mod...)
		chinese_remainder(value.(rs), devisor.(rs))
	end
end

# ╔═╡ e872ca19-4635-4601-85c4-311330bc28d5
using Test

# ╔═╡ a97d1c41-fb19-42ce-a412-68491f02a3fa
using Primes

# ╔═╡ b91528fe-6aa1-4c58-968e-a8997e895af3
md"## Rule 0"

# ╔═╡ 8d10a46a-99ab-4847-8fcb-51318c3f15a0
md"A vertex $v$ is in $I$ or one of $N(v)$ is in $I$, otherwise it violates the tropical tensor network definition: computing the local maximum."

# ╔═╡ 3d08e981-b041-4cbf-a808-76532ccc9359
code1 = ein"v,va,vb,vc->abc"

# ╔═╡ 83fb52ab-634e-48e5-8231-6196c240fb1e
mis_solve(code1)

# ╔═╡ 10737ca0-5419-432b-ba10-6e63d6ddddf0
vizeinsum(code1, ['v'=>(0.0, 0.0), 'a'=>(0.0, 0.5), 'b'=>(0.5, 0.0), 'c'=>(0.5, 0.5)], graphsize=5cm, unit=5.0, textoffset=(0.07, 0.0), rescale=0.9)

# ╔═╡ 57033470-7163-4c10-bf2b-5126f059a1eb
let
	img = canvas() do
		unit = 4.0
		textcolor = "black"
		nb = nodestyle(:circle, fill("black"), linewidth(0mm); r=0.005*unit)
		bt = nodestyle(:square, fill("black"), linewidth(0mm); r=0.015*unit)
		bt1 = nodestyle(:circle, fill("transparent"), stroke("black"), linewidth(0.08mm*unit); r=0.009*unit)
		nb2 = nodestyle(:circle, fill("red"), linewidth(0mm); r=0.005*unit)
		eb = bondstyle(:default, linewidth(0.08mm*unit))
		tb = textstyle(:default, fill(textcolor), fontsize(3pt*unit))
		y = 0.2
		nb >> (0.25, y)
		eb >> ((0.1, y), (0.4, y))
		tb >> ((1.0, y), "vertex")
		y = 0.4
		nb2 >> (0.25, y)
		eb >> ((0.1, y), (0.4, y))
		tb >> ((1.0, y), "boundary vertex")
		y = 0.6
		bt1 >> (0.25, y)
		#eb >> ((0.1, y), (0.4, y))
		tb >> ((1.0, y), "vertex tensor")
		y = 0.8
		bt >> (0.25, y)
		eb >> ((0.1, y), (0.4, y))
		tb >> ((1.0, y), "edge tensor")
	end
	Compose.set_default_graphic_size(10cm, 5cm)
	Compose.compose(context(0, 0, 0.5, 1.0), img)
end

# ╔═╡ 327414f7-7e8d-4d74-883c-25e35adf1a9a
code2 = ein"v,va,vb,vc,a,b,c->abcv"

# ╔═╡ 1891feb3-4e3a-4fdd-aca9-5b14fef13c30
res2 = mis_solve(code2)

# ╔═╡ ecd0f2af-338c-432a-9ff8-5c7bdfb410d2
vizeinsum(code2, ['v'=>(0.0, 0.0), 'a'=>(0.0, 0.5), 'b'=>(0.5, 0.0), 'c'=>(0.5, 0.5)], graphsize=5cm, unit=5.0, textoffset=(0.08, 0.0), rescale=0.9)

# ╔═╡ 545961ab-8145-470e-bc62-e7ea7ddacbea
compress!(res2)

# ╔═╡ 5dbb8308-911d-4140-89d9-108c2c59eea0
md"## Mirrors

A vertex $w ∈ N^2(v)$ is called a mirror of $v$ if $N(v) \backslash N(w)$ is a clique.
"

# ╔═╡ 7d368640-ea6d-4943-8557-9bf2b9187209
clique = ein"i,j,k,ij,jk,ik->ijk"

# ╔═╡ 38c7499d-05f2-4789-bf63-e3b904a5602e
mis_solve(clique)

# ╔═╡ 6d3847b6-25b6-4182-a1c3-624f81e541c3
mirror = ein"v,w,i,j,k,l,m,vi,vj,vk,ij,jk,ik,wl,vl,wm,vm->ijklmw"

# ╔═╡ 979df4a5-ab35-44fc-b0ff-6133f46c4427
vizeinsum(mirror, ['v'=>(-0.5, 0.0), 'i'=>(-0.3, -0.2), 'j'=>(-0.3, -0.6), 'k'=>(0.0, -0.1), 'l'=>(0.0, 0.3), 'm'=>(0.0, 0.1), 'w'=>(0.5, 0.0)]; graphsize=10cm, unit=3.0, textoffset=(0.07, 0.0), rescale=0.9)

# ╔═╡ ebbf31e2-52d0-45c1-bbe2-eea3c2a2356e
md"Mirror lemma states that
```math
\alpha(G) = \max(1+\alpha(G\backslash N[v]),\alpha(G\backslash (M(v)\cup \{v\}))
```"

# ╔═╡ 5927c5bc-5e10-4749-a965-0964523aaabe
md"By rule 0, either $v$ or at least one of $N(v)$ is in $I$. The mirror lemma states that if $v$ is not in $M$, there exists an MIS $I$ that $M(v)\notin I$. If $v$ is not in the MIS, then there must be one of $N(v)$ (ijkml) in the MIS (rule 0). If $w$ is in $I$, then none of $N(v) \cap N(w)$ (ml) is in $I$, then there must be one of node in the clique $N(v)\backslash N(w)$ (ijk) in $I$ (rule 0), since clique has at most one node in the MIS, the tensor compression will eliminate this solution by moving the occupied node to the interia. Hence, the tensor compression rule do captures the mirror rule."

# ╔═╡ d7c60181-419b-40c4-b386-aa9d571ac3db
compress!(mis_solve(mirror))

# ╔═╡ 81180c12-8d51-4ee1-95f4-ab0fa81dc1b4
md"## Rule 2"

# ╔═╡ ca891829-e93d-4456-a29c-a5afa5ea3f4d
code3 = ein"v,a,b,w,va,vb,vw,wa,wb->bav"

# ╔═╡ 79ee475d-6d66-4368-8bf9-26a340b6d369
vizeinsum(code3, ['v'=>(-0.5, 0.0), 'a'=>(0.0, 0.5), 'b'=>(0.0, -0.5), 'w'=>(0.5, 0.0)], graphsize=5cm, unit=5.0, textoffset=(0.15, 0.0), rescale=0.9)

# ╔═╡ 43874b55-d3e5-41fa-bdad-487092191770
res3 = mis_solve(code3)

# ╔═╡ 485b0764-6e2f-49ed-981f-cadd61175f76
mis_solve(ein"v,a,b,va,vb->bav")

# ╔═╡ 4ff5c3cc-d2f4-45e8-b4f2-7212ce09a4a5
compress!(res3)

# ╔═╡ b3e073f5-9faa-4f0a-99ed-0df045302593
md"""
```fortran
if |V| = 0 then
	return 0
```
"""

# ╔═╡ 301b8224-4285-4eb3-9b1b-dee0af55729b
md"""
```fortran
if ∃v ∈ V with d(v) ≤ 1 then
	return 1 + mis2(G \ N[v])
```
"""

# ╔═╡ cb2e4113-5d10-4113-a17a-417617688aed
compress!(mis_solve(ein"v,j,vj->j"))

# ╔═╡ 47af75c7-1c23-436c-8fde-8a3027e08824
md"""
```fortran
if ∃v ∈ V with d(v) = 2 then
    (let u₁ and u₂ be the neighbors of v)
    if {u₁, u₂} ∈ E then
		return 1 + mis2(G \ N[v])
	if {u₁, u₂} ∉ E then
		if |N²(v)| = 1 then
			(let N²(v) = {w})
			return max(2 + mis2(G \ (N²[v] ∪ N[w])), 2 + mis2(G \ N²[v]))
		if |N²(v)| > 1 then
			return max(mis2(G \ N[v]), mis2(G \ (M(v) ∪ {v}))
```
"""

# ╔═╡ 743aad4f-8a03-449c-a809-188b7481dbb5
md"case: `{u₁, u₂} ∈ E`"

# ╔═╡ cf8a67c1-8a11-4d5e-b49c-8c94f7698250
compress!(mis_solve(ein"v,i,j,vi,vj,ij->ij"))

# ╔═╡ 4aca79a5-2bad-4dd1-8453-ecc342cb5e1d
md"case: `{u₁, u₂} ∉ E and |N₂(v)=1|`"

# ╔═╡ 036c0a7e-4b0f-4e91-aa76-abb4a25b31d7
compress!(mis_solve(ein"v,i,j,vi,vj,iw->w"))

# ╔═╡ 4d9a5bfc-54c8-4e8a-b0cb-d6f505748395
compress!(mis_solve(ein"v,i,j,vi,vj,jw->w"))

# ╔═╡ 9f803d43-fdd0-4e38-a80c-f976da616abf
compress!(mis_solve(ein"v,i,j,vi,vj,jw,iw->w"))

# ╔═╡ 0b5186d6-749f-4a43-bca9-26580ebceb95
md"Note: the first term is not nessesary"

# ╔═╡ 15a090b6-ac7e-4824-bbe6-b656e3515a04
md"case: `{u₁, u₂} ∉ E and |N₂(v)>1|` is the mirror rule."

# ╔═╡ 56be5571-bb89-49bf-822b-458575c76e61
md"""
```fortran
if ∃v ∈V with d(v) = 3 then
	(let u₁u₂ and u₃ be the neighbors of v)
	if G[N(v)] has no edge then
		if v has a mirror then
			return max(1+mis2(G\N[v]),mis2(G\(M(v)∪{v}))
		if v has no mirror then
			return max(1+mis2(G\N[v]),2+mis2(G\N[{u₁,u₂}]),
                       2+mis2(G\(N[{u₁,u₃}]∪{u₂})),
                       2+mis2(G\(N[{u₂,u₃}]∪{u₁})))
	if G[N(v)] has one or two edges then
		return max(1+mis2(G\N[v]),mis2(G\(M(v)∪{v}))
	if G[N(v)] has three edges then
		return 1+mis2(G\N[v])
```
"""

# ╔═╡ a82a478e-6e46-445a-af4a-3da4a5c27f46
md"case: G[N(v)] has no edge and v has a mirror is a mirror rule."

# ╔═╡ 7b444fb1-7dc9-42ab-80a4-fc63005b4ead
md"case: G[N(v)] has no edge and v has no mirror."

# ╔═╡ b79dddc2-d4a3-434a-b093-5de8c54e7baa
code311 = ein"v,i,j,k,vi,vj,vk->ijk"

# ╔═╡ c114b32e-6710-4a8e-bc8e-d190f0a3a1b4
vizeinsum(code311, ['v'=>(-0.5, 0.0), 'i'=>(-0.2, -0.3), 'j'=>(-0.2, 0.3), 'k'=>(-0.1, 0.0)]; graphsize=10cm, unit=3.0, textoffset=(0.07, 0.0), rescale=0.9)

# ╔═╡ 63255dae-8a58-42f0-91b9-95bbb030cd24
compress!(mis_solve(code311))

# ╔═╡ 1dc38a7c-d874-4b6c-b534-4294cc673c40
md"The entries $T_{111} = 3.0$ and $T_{110} = 2.0$ corresponds to `2+mis2(G\N[{u₁,u₂}])`, $T_{011} = 2.0$ and $T_{101} = 2.0$ correspond to `2+mis2(G\(N[{u₂,u₃}]∪{u₁})` and `2+mis2(G\(N[{u₁,u₃}]∪{u₂}))`"

# ╔═╡ 761205b6-ec01-4472-93c1-e7a63b566351
md"case: G[N(v)] has one or two edges is the mirror rule."

# ╔═╡ 9d0823c9-488c-4363-972b-3b7b1f7535a8
md"case: G[N(v)] has three edges"

# ╔═╡ d6fd7c0d-03dc-451e-ba92-58934100b7b4
code331 = ein"v,i,j,k,vi,vj,vk,ij,jk,ik->ijk"

# ╔═╡ b9e9d629-5325-4c97-a65a-13297f9342b0
compress!(mis_solve(code331))

# ╔═╡ 468ae3cb-dc4e-4726-b4ab-a088827e47bb
md"""
```fortran
if ∆(G) ≥ 6 then
	choose a vertex v of maximum degree in G
	return max(1+mis2(G\N[v]),mis2(G\v))
```
"""

# ╔═╡ 9f4e7d1b-2d2f-4382-9de8-852bbff1e55b
md"This is rule 0"

# ╔═╡ c7df0876-1a65-44a3-b372-4f32af50672f
md"""
```fortran
if G is disconnected then
	(let C ⊆V be a component of G)
	return mis2(G[C])+mis2(G\C)
```
"""

# ╔═╡ 41ca79c5-0d35-480b-86f5-37c0d33410db
md"This is obvious"

# ╔═╡ c71a62e8-ee77-4b41-b321-9e47db50cbb4
md"""
```fortran
if G is 4 or 5-regular then
	choose any vertex v of G
	return max(1+mis2(G\N[v]),mis2(G\ (M(v)∪{v}))
```
"""

# ╔═╡ c419fcce-2a09-4871-90d6-3af8673520ee
md"This is the mirror rule."

# ╔═╡ e1f474e6-4dcb-4934-808b-3719de3f906a
md"""
```fortran
if ∆(G) = 5 and δ(G) = 4 then
	choose adjacent vertices v and w with d(v) = 5 and d(w) = 4 in G
	return max(1+mis2(G\N[v]),1+mis2(G \ ({v}∪M(v)∪N[w])), mis2(G \ (M(v)∪{v,w})))
```
"""

# ╔═╡ 90d676d8-4139-433c-afcb-dc3c16cc0be2
md"This is a combination of mirror rule 0 and mirror rule by deviding the problem into 3 categories: $v$ is the set, $v$ is not in the set but $w$ is in the set and both $v$ and $w$ are not in the set. After contracting and compress $v$ and $w$ as well as their environment, ."

# ╔═╡ c8dbb330-efdd-4067-acd7-ed0fc1facf48
md"## Independance polynomial"

# ╔═╡ 08b92a4a-3db4-4548-9516-e5a513f78237
function random_regular_eincode(n, k)
	g = LightGraphs.random_regular_graph(n, k)
	ixs = [minmax(e.src,e.dst) for e in LightGraphs.edges(g)]
	code = EinCode((ixs..., [(i,) for i in LightGraphs.vertices(g)]...), ())
end

# ╔═╡ 1b8bb513-7e2a-4077-b0e9-ead9764a93fd
@bind n NumberField(2:80; default=10)

# ╔═╡ a3b8e20a-4214-4f12-ad60-de699718b727
code = random_regular_eincode(n, 3)

# ╔═╡ db7e73f2-d1e8-4c77-b4c4-3f1f936915f2
optcode = optimize_greedy(code, NoteOnTropicalMIS.uniformsize(code, 2))

# ╔═╡ 6781f8f5-97b5-43ae-9afc-f914f53e170e
OMEinsum.timespace_complexity(optcode, uniformsize(optcode, 2))

# ╔═╡ 1a63c2ea-0b51-4900-8fbb-2aaec508fcd3
mis_result = mis_solve(optcode)[].n |> Int

# ╔═╡ 916ed7f0-c98f-4f4e-b5c9-32e4eeef2ccd
mis_count(optcode)

# ╔═╡ 99609748-c999-4903-b31b-d883a2ce3250
begin
_auto_mispolytensor(x::T, ix::NTuple{2}) where T = T[1 1; 1 0]
_auto_mispolytensor(x::T, ix::NTuple{1}) where T = T[1, x]
function generate_polyxs!(x::T, code::OMEinsum.NestedEinsum, xs) where {T}
    for (ix, arg) in zip(OMEinsum.getixs(code.eins), code.args)
		if arg isa Integer
			xs[arg] = _auto_mispolytensor(x, ix)
		else
        	generate_polyxs!(x, arg, xs)
		end
    end
	return xs
end

function generate_polyxs!(x::T, code::EinCode, xs) where {T}
    for (i,ix) in enumerate(OMEinsum.getixs(code))
		xs[i] = _auto_mispolytensor(x, ix)
    end
	return xs
end
end

# ╔═╡ f48bf3da-acfd-4d07-ac72-daaacbe987b7
function mis_polysolve(code, x::T) where {T}
	xs = generate_polyxs!(x, code, Vector{Any}(undef, NoteOnTropicalMIS.ninput(code)))
	code(xs...)
end

# ╔═╡ 52077825-a293-4bd7-91bc-8bed0c3a447a
md"## Polynomial fit"

# ╔═╡ 75e06a22-7107-45fb-adb7-7f3888b2f6f1
let
	xs = (0:mis_result)
	ys = [mis_polysolve(optcode, x)[] for x in xs]
	fit(xs, ys, mis_result)
	#@show xs, ys
end

# ╔═╡ daa47c15-207d-49c1-aa1c-fd5e148d5a55
md"## Fourier"

# ╔═╡ 826e6fd6-3a7c-4aef-896a-667582b744b7
let
	ω = exp(-2im*π/(mis_result+1))
	xs = 5.0 .* collect(ω .^ (0:mis_result))
	ys = [mis_polysolve(optcode, x)[] for x in xs]
	ifft(ys) ./ (5.0 .^ (0:mis_result))
end

# ╔═╡ b9f15328-a0f2-4538-9758-fe7e2a2190d7
md"## Symbolic"

# ╔═╡ ec4ef314-d42b-4a52-ae75-d109dab02426
OMEinsum.asarray(x::Polynomial, y::AbstractArray{<:Polynomial}) = fill(x)

# ╔═╡ 2f20b087-e990-44d5-abaf-c6605d1515ed
mis_polysolve(optcode, Polynomial([0, 1]))

# ╔═╡ 5f422a16-c709-472a-b380-5b49fa0085e7
md"## Chinese remainder theorem"

# ╔═╡ 7fa19fa3-9f88-47c2-b74a-c98981c75697
@test chinese_remainder([2,3,2], [3,5,7]) == 23

# ╔═╡ 675af106-c92c-45df-bfd5-3e3fb0778842
@test CRT(Mod{3}(2), Mod{5}(3), Mod{7}(2)) == 23

# ╔═╡ 4efa563f-8fd4-48bd-8d92-057209538559
CRT(Mod{3}(2), Mod{5}(3), Mod{13}(2))

# ╔═╡ 850c0c70-7256-4b6e-ac13-6b28fb5b71da
# a set of numbers prime to each other
function bigest_prime_table(::Type{T}, n::Int) where T
	res = T[]
	for k = typemax(T):-1:T(2)
		if isempty(res) || all(x->gcd(k, x)==1, res)
			push!(res, k)
		end
		if length(res) == n
			return res
		end
	end
	return res
end

# ╔═╡ 513616a5-4edf-45d6-87fc-f3d1df222b8b
bigest_prime_table(Int, 3)

# ╔═╡ 830384d5-d729-46f1-bbc7-2cc8a02cd1e1
@test bigest_prime_table(Int, 4) == [typemax(Int):-1:typemax(Int)-2..., typemax(Int)-6]

# ╔═╡ 33d82884-0585-4a48-b734-f97186454058
YS = let
	xs = (0:mis_result)
	ys = [mis_polysolve(optcode, x)[] for x in xs]
	
	YS = []
	for P in bigest_prime_table(Int, 5)
		ys = [mis_polysolve(optcode, Mod{P,Int}(x))[] for x in xs]
		push!(YS, ys)
	end
	YS
	#fit(xs, ys, mis_result)
end

# ╔═╡ 9a9e18fd-cfcf-4321-934a-c681d770e2ea
r5 = map(yi->chinese_remainder(BigInt.(Mods.value.(yi)), BigInt.(Mods.modulus.(yi))), zip(YS[1:4]...))

# ╔═╡ c6ce1a10-81ca-47d7-b5b2-8bff52e55b40
r4 = map(yi->chinese_remainder(BigInt.(Mods.value.(yi)), BigInt.(Mods.modulus.(yi))), zip(YS[1:4]...))

# ╔═╡ 3080838b-22ec-44c3-a60e-2b2531bef037
r3 = map(yi->chinese_remainder(BigInt.(Mods.value.(yi)), BigInt.(Mods.modulus.(yi))), zip(YS[1:3]...))

# ╔═╡ e25f34bf-1611-49e1-8933-16f245a0c460
r2 = map(yi->chinese_remainder(BigInt.(Mods.value.(yi)), BigInt.(Mods.modulus.(yi))), zip(YS[1:2]...))

# ╔═╡ fc6c553b-7b6f-4abf-9136-ee428462516f
r1 = map(yi->chinese_remainder(BigInt.(Mods.value.(yi)), BigInt.(Mods.modulus.(yi))), YS[1])

# ╔═╡ bb42928b-0589-425a-bcbd-346ed9a5ef76
r4 == r3

# ╔═╡ c13f5fe7-3a2c-4369-a18b-a7388190bd4a
md"## Ring arithemetics"

# ╔═╡ 02918811-f120-4291-809b-caeef49468f2
function gaussian_eliminate(A, b)
	A \ b
end

# ╔═╡ f8b0dffb-a6af-4f7d-a608-013c50ffda75
A, b = let
	xs = 0:mis_result
	N = prevprime(typemax(Int))
	T = Mod{N,Int}
	ys = [mis_polysolve(optcode, T(x))[] for x in xs]
	A = zeros(T, mis_result+1, mis_result+1)
	for j=1:mis_result+1, i=1:mis_result+1
		A[j,i] = T(xs[j])^(i-1)
	end
	A, T.(ys)
end

# ╔═╡ a458b905-9434-43ab-b62d-c653c276b9f5
Base.abs(x::Mod) = x

# ╔═╡ d30c76c8-2bde-405a-aed2-11b46f919394
Base.isless(x::Mod{N}, y::Mod{N}) where N = mod(x.val, N) < mod(y.val, N)

# ╔═╡ 0319a259-4bcc-4ad2-a2e7-bb07882a5670
A \ b

# ╔═╡ Cell order:
# ╠═c8c14900-b7f5-11eb-3dbe-4390dd0fcc0a
# ╟─b91528fe-6aa1-4c58-968e-a8997e895af3
# ╟─8d10a46a-99ab-4847-8fcb-51318c3f15a0
# ╠═3d08e981-b041-4cbf-a808-76532ccc9359
# ╠═83fb52ab-634e-48e5-8231-6196c240fb1e
# ╠═10737ca0-5419-432b-ba10-6e63d6ddddf0
# ╟─57033470-7163-4c10-bf2b-5126f059a1eb
# ╠═327414f7-7e8d-4d74-883c-25e35adf1a9a
# ╠═1891feb3-4e3a-4fdd-aca9-5b14fef13c30
# ╠═ecd0f2af-338c-432a-9ff8-5c7bdfb410d2
# ╠═545961ab-8145-470e-bc62-e7ea7ddacbea
# ╟─5dbb8308-911d-4140-89d9-108c2c59eea0
# ╠═7d368640-ea6d-4943-8557-9bf2b9187209
# ╠═38c7499d-05f2-4789-bf63-e3b904a5602e
# ╠═6d3847b6-25b6-4182-a1c3-624f81e541c3
# ╠═979df4a5-ab35-44fc-b0ff-6133f46c4427
# ╟─ebbf31e2-52d0-45c1-bbe2-eea3c2a2356e
# ╟─5927c5bc-5e10-4749-a965-0964523aaabe
# ╠═d7c60181-419b-40c4-b386-aa9d571ac3db
# ╟─81180c12-8d51-4ee1-95f4-ab0fa81dc1b4
# ╠═ca891829-e93d-4456-a29c-a5afa5ea3f4d
# ╠═79ee475d-6d66-4368-8bf9-26a340b6d369
# ╠═43874b55-d3e5-41fa-bdad-487092191770
# ╠═485b0764-6e2f-49ed-981f-cadd61175f76
# ╠═4ff5c3cc-d2f4-45e8-b4f2-7212ce09a4a5
# ╟─b3e073f5-9faa-4f0a-99ed-0df045302593
# ╟─301b8224-4285-4eb3-9b1b-dee0af55729b
# ╠═cb2e4113-5d10-4113-a17a-417617688aed
# ╟─47af75c7-1c23-436c-8fde-8a3027e08824
# ╟─743aad4f-8a03-449c-a809-188b7481dbb5
# ╠═cf8a67c1-8a11-4d5e-b49c-8c94f7698250
# ╟─4aca79a5-2bad-4dd1-8453-ecc342cb5e1d
# ╠═036c0a7e-4b0f-4e91-aa76-abb4a25b31d7
# ╠═4d9a5bfc-54c8-4e8a-b0cb-d6f505748395
# ╠═9f803d43-fdd0-4e38-a80c-f976da616abf
# ╟─0b5186d6-749f-4a43-bca9-26580ebceb95
# ╟─15a090b6-ac7e-4824-bbe6-b656e3515a04
# ╟─56be5571-bb89-49bf-822b-458575c76e61
# ╟─a82a478e-6e46-445a-af4a-3da4a5c27f46
# ╟─7b444fb1-7dc9-42ab-80a4-fc63005b4ead
# ╠═b79dddc2-d4a3-434a-b093-5de8c54e7baa
# ╠═c114b32e-6710-4a8e-bc8e-d190f0a3a1b4
# ╠═63255dae-8a58-42f0-91b9-95bbb030cd24
# ╟─1dc38a7c-d874-4b6c-b534-4294cc673c40
# ╟─761205b6-ec01-4472-93c1-e7a63b566351
# ╟─9d0823c9-488c-4363-972b-3b7b1f7535a8
# ╠═d6fd7c0d-03dc-451e-ba92-58934100b7b4
# ╠═b9e9d629-5325-4c97-a65a-13297f9342b0
# ╟─468ae3cb-dc4e-4726-b4ab-a088827e47bb
# ╟─9f4e7d1b-2d2f-4382-9de8-852bbff1e55b
# ╟─c7df0876-1a65-44a3-b372-4f32af50672f
# ╟─41ca79c5-0d35-480b-86f5-37c0d33410db
# ╟─c71a62e8-ee77-4b41-b321-9e47db50cbb4
# ╟─c419fcce-2a09-4871-90d6-3af8673520ee
# ╟─e1f474e6-4dcb-4934-808b-3719de3f906a
# ╟─90d676d8-4139-433c-afcb-dc3c16cc0be2
# ╟─c8dbb330-efdd-4067-acd7-ed0fc1facf48
# ╠═d1a8f940-c8b8-4ce6-af68-dee3a286e632
# ╠═08b92a4a-3db4-4548-9516-e5a513f78237
# ╠═8e2263d2-ffbf-4aa5-bf7e-29f16b2dbade
# ╠═1b8bb513-7e2a-4077-b0e9-ead9764a93fd
# ╟─a3b8e20a-4214-4f12-ad60-de699718b727
# ╠═db7e73f2-d1e8-4c77-b4c4-3f1f936915f2
# ╠═6781f8f5-97b5-43ae-9afc-f914f53e170e
# ╠═1a63c2ea-0b51-4900-8fbb-2aaec508fcd3
# ╠═916ed7f0-c98f-4f4e-b5c9-32e4eeef2ccd
# ╠═f48bf3da-acfd-4d07-ac72-daaacbe987b7
# ╠═99609748-c999-4903-b31b-d883a2ce3250
# ╟─52077825-a293-4bd7-91bc-8bed0c3a447a
# ╠═e56f2ba7-6a50-4400-b761-9f587076778f
# ╠═75e06a22-7107-45fb-adb7-7f3888b2f6f1
# ╟─daa47c15-207d-49c1-aa1c-fd5e148d5a55
# ╠═f2c1633a-0e2b-4ecc-81bf-f54fa0efcb44
# ╠═826e6fd6-3a7c-4aef-896a-667582b744b7
# ╟─b9f15328-a0f2-4538-9758-fe7e2a2190d7
# ╠═ec4ef314-d42b-4a52-ae75-d109dab02426
# ╠═2f20b087-e990-44d5-abaf-c6605d1515ed
# ╟─5f422a16-c709-472a-b380-5b49fa0085e7
# ╠═0a5712bd-4f8f-41a4-8fa3-7f947a078a4a
# ╠═e872ca19-4635-4601-85c4-311330bc28d5
# ╠═7fa19fa3-9f88-47c2-b74a-c98981c75697
# ╠═675af106-c92c-45df-bfd5-3e3fb0778842
# ╠═4efa563f-8fd4-48bd-8d92-057209538559
# ╠═850c0c70-7256-4b6e-ac13-6b28fb5b71da
# ╠═513616a5-4edf-45d6-87fc-f3d1df222b8b
# ╠═830384d5-d729-46f1-bbc7-2cc8a02cd1e1
# ╠═33d82884-0585-4a48-b734-f97186454058
# ╠═9a9e18fd-cfcf-4321-934a-c681d770e2ea
# ╠═c6ce1a10-81ca-47d7-b5b2-8bff52e55b40
# ╠═3080838b-22ec-44c3-a60e-2b2531bef037
# ╠═e25f34bf-1611-49e1-8933-16f245a0c460
# ╠═fc6c553b-7b6f-4abf-9136-ee428462516f
# ╠═bb42928b-0589-425a-bcbd-346ed9a5ef76
# ╟─c13f5fe7-3a2c-4369-a18b-a7388190bd4a
# ╠═02918811-f120-4291-809b-caeef49468f2
# ╠═a97d1c41-fb19-42ce-a412-68491f02a3fa
# ╠═f8b0dffb-a6af-4f7d-a608-013c50ffda75
# ╠═a458b905-9434-43ab-b62d-c653c276b9f5
# ╠═d30c76c8-2bde-405a-aed2-11b46f919394
# ╠═0319a259-4bcc-4ad2-a2e7-bb07882a5670
