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
begin
	using Revise, OMEinsum, OMEinsumContractionOrders, TropicalNumbers, NoteOnTropicalMIS, Compose, Viznet
	using OMEinsumContractionOrders: uniformsize
	function savetodisk(obj, filename)
		obj |> SVG(filename*".svg")
	end
end;

# ╔═╡ e56f2ba7-6a50-4400-b761-9f587076778f
using Polynomials

# ╔═╡ d1a8f940-c8b8-4ce6-af68-dee3a286e632
using LightGraphs

# ╔═╡ 8e2263d2-ffbf-4aa5-bf7e-29f16b2dbade
using PlutoUI

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

# ╔═╡ 2a37780e-53b6-4a40-9738-5daba20c1001
begin
	node_black = nodestyle(:circle, r=0.03)
	node_white = nodestyle(:circle, fill("white"), r=0.08)
	edge_black = bondstyle(:default)
	text_black = textstyle(:default)
	text_white = textstyle(:default, fill("white"))
	
	function save_pdf(img, fname, save)
		if save
			img |> SVG("$fname.svg")
		end
		run(`rsvg-convert -f pdf -o $fname.pdf $fname.svg`)
		return img
	end
end;

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

# ╔═╡ 06f4c475-685d-4127-bb1b-e2ed1b9a5cd4
md" $(@bind x CheckBox()) save to disk"

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

# ╔═╡ 29e224b6-0b4c-4a83-811e-6d6cde561278
md" $(@bind save_fig12 CheckBox()) save fig1"

# ╔═╡ e2082ba8-bc1a-4114-a43a-c8b430d303fd
let
	img = canvas() do
		Compose.set_default_graphic_size(8cm, 8cm)
		locs = [(0.2, 0.2), (0.2, 0.6), (0.6,0.6), (0.6, 0.2), (0.9, 0.9)]
		for loc in locs
			node_black >> loc
		end
		for (loc, t) in zip(locs, "abcde")
			text_white >> (loc, string(t))
		end
		for (i,j) in [(1, 2), (2, 3), (3,4), (4,1), (3,5), (2,4)]
			edge_black >> (locs[i], locs[j])
		end
	end
	save_pdf(img, "fig1", save_fig12)
end

# ╔═╡ b91a8b8e-8bc5-4ebc-addc-0893144f2b1a
let
	img = canvas() do
		Compose.set_default_graphic_size(8cm, 8cm)
		locs = [(0.2, 0.2), (0.2, 0.6), (0.6,0.6), (0.6, 0.2), (0.9, 0.9)]
		for loc in locs
			text_black >> (loc, "[1, x]")
		end
		for (i,j) in [(1, 2), (2, 3), (3,4), (4,1), (3,5), (2,4)]
			text_black >> ((locs[i] .+ locs[j]) ./ 2, "[1 1\n1 0]")
		end
		for loc in locs
			node_white >> loc
		end
		for (i,j) in [(1, 2), (2, 3), (3,4), (4,1), (3,5), (2,4)]
			node_white >> ((locs[i] .+ locs[j]) ./ 2)
		end
		for (i,j) in [(1, 2), (2, 3), (3,4), (4,1), (3,5), (2,4)]
			edge_black >> (locs[i], locs[j])
		end
	end
	save_pdf(img, "fig2", save_fig12)
end

# ╔═╡ c1910146-a9bb-472a-814c-419aa75b077e
md"""
The result is
```math
I(G, x) = \sum_{k=1}^{\alpha(G)} c_k x^k,
```
where ``c_k`` is the the degeneracy for independant set size $k$.

[Computing the Independence Polynomial: from the Tree Threshold down to the Roots](https://arxiv.org/abs/1608.02282)


[[Ferrin, Gregory Matthew. "Independence polynomials." (2014)]](https://scholarcommons.sc.edu/cgi/viewcontent.cgi?article=3632&context=etd).
"""

# ╔═╡ c5d0d8c7-0fc2-4f32-803b-a15d836e80bb
democode = ein"a,b,c,d,ab,ad,bc,bd,cd,ce->"

# ╔═╡ 8b4b532e-f149-49e0-be60-7b3d674ba4d8
demo_xs = NoteOnTropicalMIS.generate_polyxs!(2.0, democode)

# ╔═╡ 1e1c2e60-8c8f-4b2c-b615-1ad4497c0e5a
f_2 = democode(demo_xs...)

# ╔═╡ 598c0ff3-786e-4911-a6f7-1f1d4956e154
let
	xs = (0:2)
	ys = map(xs) do x
		demo_xs = NoteOnTropicalMIS.generate_polyxs!(x, democode)
		democode(demo_xs...)[]
	end
	fit(xs, ys, 2)
end

# ╔═╡ 0daab8e4-f099-414f-aca3-51b845ef88bd
md"# Random 3 regular graph"

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
optcode = optimize_kahypar(code, uniformsize(code, 2); sc_target=20, max_group_size=30)

# ╔═╡ 6781f8f5-97b5-43ae-9afc-f914f53e170e
OMEinsum.timespace_complexity(optcode, uniformsize(optcode, 2))

# ╔═╡ 1a63c2ea-0b51-4900-8fbb-2aaec508fcd3
mis_result = mis_solve(optcode)[].n |> Int

# ╔═╡ 916ed7f0-c98f-4f4e-b5c9-32e4eeef2ccd
mis_count(optcode)

# ╔═╡ 52077825-a293-4bd7-91bc-8bed0c3a447a
md"## Polynomial fit"

# ╔═╡ 25544e9e-7b20-41d2-932d-7557389defbe
md"""
```math
\begin{align}
a_0 + a_1 x_1 + a_1 x_1^2 + \ldots + a_m x_1^m &= y_0\\
a_0 + a_1 x_2 + a_2 x_2^2 + \ldots + a_m x_2^m &= y_1\\
\ldots&\\
a_0 + a_1 x_m + a_2 x_m^2 + \ldots + a_m x_m^m& = y_m
\end{align}
```
"""

# ╔═╡ 75e06a22-7107-45fb-adb7-7f3888b2f6f1
let
	xs = (0:mis_result)
	ys = [mis_polysolve(optcode, x)[] for x in xs]
	fit(xs, ys, mis_result)
	#@show xs, ys
end

# ╔═╡ daa47c15-207d-49c1-aa1c-fd5e148d5a55
md"## Fourier"

# ╔═╡ 003bd44a-b0be-4f84-94b6-7e1ae9490bcd
md"""

```math
\begin{align}
\left(\begin{matrix}
1 & x_1 & x_1^2 & \ldots & x_1^m \\
1 & x_2 & x_2^2 & \ldots & x_2^m \\
\vdots & \vdots & \vdots &\ddots & \vdots \\
1 & x_m & x_m^2 & \ldots & x_m^m
\end{matrix}\right)
\left(\begin{matrix}
a_0 \\ a_1 \\ \vdots \\ a_m
\end{matrix}\right)
= \left(\begin{matrix}
y_0 \\ y_1 \\ \vdots \\ y_m
\end{matrix}\right)\\
\end{align}
```
"""

# ╔═╡ 00ba67fe-7f53-4b47-ad98-cfd6171c8bfe
md"""
Let $x_k = r\omega^k$
"""

# ╔═╡ 0537b8ea-8007-4428-b46a-1a7b91ba30c7
md"""
```math
\begin{align}
\left(\begin{matrix}
1 & r\omega & r^2\omega^2 & \ldots & r^m\omega^m \\
1 & r\omega^2 & r^2\omega^4 & \ldots & r^m\omega^{2m} \\
\vdots & \vdots & \vdots &\ddots & \vdots \\
1 & r\omega^m & r^2\omega^{2m} & \ldots & r^m\omega^{m^2}
\end{matrix}\right)
\left(\begin{matrix}
a_0 \\ a_1 \\ \vdots \\ a_m
\end{matrix}\right)
= \left(\begin{matrix}
y_0 \\ y_1 \\ \vdots \\ y_m
\end{matrix}\right)\\
\end{align}
```
"""

# ╔═╡ 8502b87f-3e59-4b58-ab4b-2814537d4979
md"When $r=1$, the left side is a DFT matrix. We can obtain the factors using the relation $\vec a = {\rm FFT^{-1}}(\omega) \cdot \vec y$. In the special case that $\omega = e^{-2πi/(m+1)}$, it is directly solvable with inverse fast fourier transformation algorithm in package `FFTW`."

# ╔═╡ 2c5c9e61-2e1b-43f2-a33e-9ea621f8c092
md"""
```math
{\rm FFT}(\omega) \cdot \vec a_r = \vec y
```
where $(\vec a_r)_k = a_k r ^k$, by choosing diferent $r$, we can obtain better precision in low independant set size region  (``\omega<1``) and high independant set size region (``\omega>1``).
"""

# ╔═╡ 826e6fd6-3a7c-4aef-896a-667582b744b7
let
	r = 1.2
	ω = exp(-2im*π/(mis_result+1))
	xs = r .* collect(ω .^ (0:mis_result))
	ys = [mis_polysolve(optcode, x)[] for x in xs]
	ifft(ys) ./ (r .^ (0:mis_result))
end

# ╔═╡ b9f15328-a0f2-4538-9758-fe7e2a2190d7
md"## Symbolic"

# ╔═╡ ec4ef314-d42b-4a52-ae75-d109dab02426
OMEinsum.asarray(x::Polynomial, y::AbstractArray) = fill(x)

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
		ys = [mis_polysolve(optcode, Mod{P}(x))[] for x in xs]
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

# ╔═╡ 8a7d7bf7-3816-4619-a75e-5ee03777add7
let
	xs = (0:mis_result)
	A = zeros(BigInt, mis_result+1, mis_result+1)
	for j=1:mis_result+1, i=1:mis_result+1
		A[j,i] = BigInt(xs[j])^(i-1)
	end
	A \ r3
end

# ╔═╡ c13f5fe7-3a2c-4369-a18b-a7388190bd4a
md"## Finite field arithemetics"

# ╔═╡ a458b905-9434-43ab-b62d-c653c276b9f5
Base.abs(x::Mod) = x

# ╔═╡ d30c76c8-2bde-405a-aed2-11b46f919394
Base.isless(x::Mod{N}, y::Mod{N}) where N = mod(x.val, N) < mod(y.val, N)

# ╔═╡ f8b0dffb-a6af-4f7d-a608-013c50ffda75
let
	xs = 0:mis_result
	N = prevprime(typemax(Int))
	T = Mod{N}
	ys = [mis_polysolve(optcode, T(x))[] for x in xs]
	A = zeros(T, mis_result+1, mis_result+1)
	for j=1:mis_result+1, i=1:mis_result+1
		A[j,i] = T(xs[j])^(i-1)
	end
	A \ T.(ys)
end

# ╔═╡ 04563702-ddb0-437a-8ae5-4132a13d1f49
prevprime(1000)

# ╔═╡ c69f3eb5-aec6-496b-8872-3f92d295f93b
prevprime(typemax(Int))

# ╔═╡ bd980a2e-f60a-4bde-9275-0737ab20d3c2
Mod{997}(40) + Mod{997}(2341234)

# ╔═╡ 01177d08-50ba-4b6a-8a46-bd446ba8c56e
Mod{997}(40) * Mod{997}(2341234)

# ╔═╡ 5eee76e3-d8f4-4ab6-bf8f-8493d91b4a50
Mod{997}(1) == Mod{997}(998)

# ╔═╡ 0ec2b495-07f3-40c8-833a-01223514521f
Mod{997}(40) / Mod{997}(2341234)

# ╔═╡ 13cb14b6-50e6-4f48-8243-b8793bdd4cbf
inv(Mod{997}(37))

# ╔═╡ f14b15c2-f066-4a2e-a4fe-5e5fbe9cf987
md"# Contraction order comparison"

# ╔═╡ ba163c37-011c-455b-9899-38e9a108b6a6
md"number of nodes = $(@bind number_of_nodes Slider(4:2:500; show_value=true))"

# ╔═╡ 33c23f29-6e7f-4a8e-b164-69b5605a3f3b
let
	code = random_regular_eincode(number_of_nodes, 3)
	optcode = optimize_kahypar(code, uniformsize(code, 2); sc_target=number_of_nodes^0.7, max_group_size=50)
	OMEinsum.timespace_complexity(optcode, uniformsize(optcode, 2))
end

# ╔═╡ b1fcf86d-6e2c-4f1c-8ad5-910310b08165
let
	code = random_regular_eincode(number_of_nodes, 3)
	optcode = optimize_greedy(code, uniformsize(code, 2); nrepeat=10, method=OMEinsum.MinSpaceOut())
	OMEinsum.timespace_complexity(optcode, uniformsize(optcode, 2))
end

# ╔═╡ e31cd630-e972-4a1e-9c40-06f0bd63731d
mis_solve(ein"i,j,v,vi,vj,ij->ijv")

# ╔═╡ Cell order:
# ╠═c8c14900-b7f5-11eb-3dbe-4390dd0fcc0a
# ╠═2a37780e-53b6-4a40-9738-5daba20c1001
# ╟─b91528fe-6aa1-4c58-968e-a8997e895af3
# ╟─8d10a46a-99ab-4847-8fcb-51318c3f15a0
# ╠═3d08e981-b041-4cbf-a808-76532ccc9359
# ╠═83fb52ab-634e-48e5-8231-6196c240fb1e
# ╠═10737ca0-5419-432b-ba10-6e63d6ddddf0
# ╟─06f4c475-685d-4127-bb1b-e2ed1b9a5cd4
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
# ╟─29e224b6-0b4c-4a83-811e-6d6cde561278
# ╟─e2082ba8-bc1a-4114-a43a-c8b430d303fd
# ╟─b91a8b8e-8bc5-4ebc-addc-0893144f2b1a
# ╟─c1910146-a9bb-472a-814c-419aa75b077e
# ╠═c5d0d8c7-0fc2-4f32-803b-a15d836e80bb
# ╠═8b4b532e-f149-49e0-be60-7b3d674ba4d8
# ╠═1e1c2e60-8c8f-4b2c-b615-1ad4497c0e5a
# ╠═e56f2ba7-6a50-4400-b761-9f587076778f
# ╠═598c0ff3-786e-4911-a6f7-1f1d4956e154
# ╟─0daab8e4-f099-414f-aca3-51b845ef88bd
# ╠═d1a8f940-c8b8-4ce6-af68-dee3a286e632
# ╠═08b92a4a-3db4-4548-9516-e5a513f78237
# ╠═8e2263d2-ffbf-4aa5-bf7e-29f16b2dbade
# ╟─1b8bb513-7e2a-4077-b0e9-ead9764a93fd
# ╟─a3b8e20a-4214-4f12-ad60-de699718b727
# ╠═db7e73f2-d1e8-4c77-b4c4-3f1f936915f2
# ╠═6781f8f5-97b5-43ae-9afc-f914f53e170e
# ╠═1a63c2ea-0b51-4900-8fbb-2aaec508fcd3
# ╠═916ed7f0-c98f-4f4e-b5c9-32e4eeef2ccd
# ╟─52077825-a293-4bd7-91bc-8bed0c3a447a
# ╟─25544e9e-7b20-41d2-932d-7557389defbe
# ╠═75e06a22-7107-45fb-adb7-7f3888b2f6f1
# ╟─daa47c15-207d-49c1-aa1c-fd5e148d5a55
# ╟─003bd44a-b0be-4f84-94b6-7e1ae9490bcd
# ╟─00ba67fe-7f53-4b47-ad98-cfd6171c8bfe
# ╟─0537b8ea-8007-4428-b46a-1a7b91ba30c7
# ╟─8502b87f-3e59-4b58-ab4b-2814537d4979
# ╟─2c5c9e61-2e1b-43f2-a33e-9ea621f8c092
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
# ╠═8a7d7bf7-3816-4619-a75e-5ee03777add7
# ╟─c13f5fe7-3a2c-4369-a18b-a7388190bd4a
# ╠═a97d1c41-fb19-42ce-a412-68491f02a3fa
# ╠═a458b905-9434-43ab-b62d-c653c276b9f5
# ╠═d30c76c8-2bde-405a-aed2-11b46f919394
# ╠═f8b0dffb-a6af-4f7d-a608-013c50ffda75
# ╠═04563702-ddb0-437a-8ae5-4132a13d1f49
# ╠═c69f3eb5-aec6-496b-8872-3f92d295f93b
# ╠═bd980a2e-f60a-4bde-9275-0737ab20d3c2
# ╠═01177d08-50ba-4b6a-8a46-bd446ba8c56e
# ╠═5eee76e3-d8f4-4ab6-bf8f-8493d91b4a50
# ╠═0ec2b495-07f3-40c8-833a-01223514521f
# ╠═13cb14b6-50e6-4f48-8243-b8793bdd4cbf
# ╟─f14b15c2-f066-4a2e-a4fe-5e5fbe9cf987
# ╟─ba163c37-011c-455b-9899-38e9a108b6a6
# ╠═33c23f29-6e7f-4a8e-b164-69b5605a3f3b
# ╠═b1fcf86d-6e2c-4f1c-8ad5-910310b08165
# ╠═e31cd630-e972-4a1e-9c40-06f0bd63731d
