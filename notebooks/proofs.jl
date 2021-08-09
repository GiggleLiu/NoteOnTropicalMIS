### A Pluto.jl notebook ###
# v0.14.2

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
	using Pkg
	Pkg.activate(pwd())
	using Revise, OMEinsum, OMEinsumContractionOrders, TropicalNumbers, NoteOnTropicalMIS, Compose, Viznet, PlutoUI
	using OMEinsumContractionOrders: uniformsize
	function savetodisk(obj, filename)
		obj |> SVG(filename*".svg")
	end
end;

# ╔═╡ 2a37780e-53b6-4a40-9738-5daba20c1001
begin
	node_black = nodestyle(:circle, r=0.03)
	node_red = nodestyle(:circle, fill("red"), r=0.03)
	node_white = nodestyle(:circle, fill("white"), r=0.08)
	edge_black = bondstyle(:default)
	text_black = textstyle(:default)
	text_white = textstyle(:default, fill("white"))
	
	function save_pdf(img, fname, save)
		if save
			img |> SVG("$fname.svg")
			run(`rsvg-convert -f pdf -o $fname.pdf $fname.svg`)
		end
		return img
	end
	
	mis_solve(code) = mis_contract(Tropical(1.0), OMEinsum.optimize_greedy(code, uniformsize(code, 2)))
	mis_configs(code) = all_config(OMEinsum.optimize_greedy(code, uniformsize(code, 2)))
end;

# ╔═╡ a8ac2b58-86c7-4964-ac10-46c842cbe731
canvas() do
	Compose.set_default_graphic_size(8cm, 8cm)
	xs = [(0.1, 0.5), (0.1, 0.7), (0.3, 0.5), (0.4, 0.8), (0.5, 0.6), (0.3, 0.3), (0.5, 0.2), (0.6, 0.7), (0.7, 0.4), (0.9, 0.5)]
	for (i, loc) in enumerate(xs)
		(i ∈ (3,4,5,6,7) ? node_red : node_black) >> loc
		text_white >> (loc, string('a'+i-1))
	end
	for (i,j) in [(1, 2), (2,3), (1,3), (3, 4), (4, 5), (5, 6), (6,7), (7,8), (8, 9), (9, 10), (3, 5), (3, 6), (7,9), (6,9), (10,8), (3,9)]
		edge_black >> (xs[i], xs[j])
	end
end

# ╔═╡ e57abcc2-7b48-4e38-a90a-fe6a6bc1eeb7
let
	img1 = canvas() do
		xs = [(0.1, 0.5), (0.1, 0.7), (0.3, 0.5), (0.4, 0.8), (0.5, 0.6), (0.3, 0.3), (0.5, 0.2), (0.6, 0.7), (0.7, 0.4), (0.9, 0.5)]
		for (i, loc) in enumerate(xs)
			if i ∈ (3,4,5,6,7)
				(node_red >> loc)
				text_white >> (loc, string('a'+i-1))
			end
		end
		for (i,j) in [(3, 4), (4, 5), (5, 6), (6,7), (3, 5), (3, 6)]
			edge_black >> (xs[i], xs[j])
		end
	end
	img2 = vizeinsum(ein"d,c,e,f,g,dc,de,ce,ef,fg,cf->cfg", ['c'=>(0.3, 0.5), 'd'=>(0.4, 0.8), 'e'=>(0.5, 0.6), 'f'=>(0.3, 0.3), 'g'=>(0.5, 0.2)], unit=2.5, rescale=0.8, textcolor="transparent")
	img3 = canvas() do
		unit = 4.0
		textcolor = "black"
		nb = nodestyle(:circle, fill("black"), linewidth(0mm); r=0.01*unit)
		bt = nodestyle(:square, fill("black"), linewidth(0mm); r=0.015*unit)
		bt1 = nodestyle(:circle, fill("transparent"), stroke("black"), linewidth(0.08mm*unit); r=0.015*unit)
		nb2 = nodestyle(:circle, fill("red"), linewidth(0mm); r=0.01*unit)
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
	img4 = canvas() do
		text_black >> ((0.2, 0.9), "(a)")
		text_black >> ((0.6, 0.9), "(b)")
	end
	Compose.set_default_graphic_size(14cm, 7cm)
	Compose.compose(context(), (context(0.0,0.0,0.5,1.0), img1), (context(0.4,0.0,0.17,1.0), img2), (context(0.6,0.3,0.2,0.4), img3), (context(), img4))
end

# ╔═╡ 2989b2da-6951-470d-8ef1-cd08628f50f4
let
	img1 = canvas() do
		xs = [(0.1, 0.5), (0.1, 0.7), (0.3, 0.5), (0.4, 0.8), (0.5, 0.6), (0.3, 0.3), (0.5, 0.2), (0.6, 0.7), (0.7, 0.4), (0.9, 0.5)]
		for (i, loc) in enumerate(xs)
			if i ∈ (3,4,5,6,7)
				(node_red >> loc)
				text_white >> (loc, string('a'+i-1))
			end
		end
		for (i,j) in [(3, 4), (4, 5), (5, 6), (6,7), (3, 5), (3, 6)]
			edge_black >> (xs[i], xs[j])
		end
	end

	img2 = canvas() do
		xs = [(0.1, 0.5), (0.1, 0.7), (0.3, 0.5), (0.4, 0.8), (0.5, 0.6), (0.3, 0.3), (0.5, 0.2), (0.6, 0.7), (0.7, 0.4), (0.9, 0.5)]
		for (i, loc) in enumerate(xs)
			if i ∈ (3,4,5,6,7)
				(node_red >> loc)
				text_white >> (loc, string('a'+i-1))
			end
		end
		for (i,j) in [(3, 4), (4, 5), (5, 6), (6,7), (3, 5), (3, 6)]
			edge_black >> (xs[i], xs[j])
		end
	end
	Compose.set_default_graphic_size(14cm, 7cm)
	Compose.compose(context(), (context(0.0,0.0,0.5,1.0), img1), (context(0.5,0.0,0.5,1.0), img2))
end

# ╔═╡ d4921f23-459a-4436-a526-71951b2d438e
mis_solve(ein"d,c,e,f,g,dc,de,ce,ef,fg,cf->cfg")

# ╔═╡ 80fb581d-c2ab-4d81-97cb-b65bd038a7ce
compress!(mis_solve(ein"d,c,e,f,g,dc,de,ce,ef,fg,cf->cfg"))

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

# ╔═╡ 3504e107-dcf3-4527-853b-b6a9fe3ced4b
md"## Basic Rule"

# ╔═╡ 300114b5-d8f7-47bb-a569-2835321f9695
md" $(@bind save_basic CheckBox()) save to disk"

# ╔═╡ 6b590957-d846-407a-99fe-b570d5e92107
basic_res = let
	code = ein"u,v,w,j,k,vw,vj,vk,vu->juwk"
	mis_solve(code)
end;

# ╔═╡ 91fb9192-22c6-4135-81c7-efa7ea681722
compress!(basic_res)

# ╔═╡ 1cb4f4fe-2b5c-4ecf-88b5-3c7fec5d9b87
basic_res2 = let
	code = ein"v,u,w,j,k,vw,vj,vk,vu->juwk"
	mis_configs(code)
end

# ╔═╡ cb1d664e-4efe-4d65-a2be-239c8cd3cfcb
md"## Satelitte Rule"

# ╔═╡ 30b2968f-413d-4e10-bb80-e71a5596dc5c
function vizgraph(nodes, edges)
	node_black = nodestyle(:circle, r=0.05)
	node_red = nodestyle(:circle, fill("red"), r=0.05)
	edge_black = bondstyle(:default)
	text_white = textstyle(:default, fill("white"))
	canvas() do
		Compose.set_default_graphic_size(8cm, 8cm)
		for (label, loc, color) in nodes
			(color=="red" ? node_red : node_black) >> loc
			text_white >> (loc, label)
		end
		for (i,j) in edges
			edge_black >> (nodes[i][2], nodes[j][2])
		end
	end
end;

# ╔═╡ e291fac5-ef7a-4308-a9d8-53adcc0bbf11
let
	img = vizgraph([("j", (0.1, 0.2), "red"), ("w", (0.3, 0.5), "red"), ("v", (0.5, 0.2), "black"), ("k", (0.7, 0.5), "red"), ("u", (0.9, 0.2), "red")], [(3,2), (4,3), (3,1), (3,5)])
	save_pdf(img, "basic", save_basic)
end

# ╔═╡ 6766e2fa-795d-47fb-9e95-47471f9cd8e2
md"``\alpha(G) = \max\{\alpha(G \backslash \{v\}), \alpha(G \backslash N[S[v]]) + |S(v)| + 1\}``"

# ╔═╡ 8c97d833-0c37-40cc-9460-58ea5f1d2e0b
md" $(@bind save_satellite CheckBox()) save to disk"

# ╔═╡ 9bdd69ac-352c-4b88-acbb-5564e57322e5
let
	img = vizgraph([("j", (0.1, 0.2), "red"), ("w", (0.3, 0.5), "black"), ("v", (0.5, 0.2), "black"), ("k", (0.7, 0.5), "red"), ("u", (0.9, 0.2), "red")], [(2,3), (2,4), (1,3), (3, 4), (4, 5), (2,5)])
	save_pdf(img, "satellite", save_satellite)
end

# ╔═╡ 687a88e1-76ad-4563-82aa-8487f65c441f
satelite_res = let
	code = ein"u,v,w,j,k,vw,vj,vk,ku,wu,wk->juk"
	mis_solve(code)
end

# ╔═╡ 23fcfe1a-e7ce-4cdc-abcb-651b8f069606
compress!(copy(satelite_res))

# ╔═╡ c9c09456-767a-4634-bca0-6897b7c5e923
md"## Mirror rule"

# ╔═╡ 70b13ee9-7282-4677-8d27-e7f3f83dc79e
md"``\alpha(G) = \max(1 + \alpha(G \backslash N[v]), \alpha(G \backslash (M(v) \cup \{v\}))``"

# ╔═╡ 792d4aaf-b247-48d7-a321-2a3ae73cb988
md" $(@bind save_mirror CheckBox()) save to disk"

# ╔═╡ a7d0bcd6-a303-4133-9bb1-ba0bb30d294b
let
	img = vizgraph([("j", (0.1, 0.2), "red"), ("w", (0.3, 0.5), "red"), ("v", (0.5, 0.2), "black"), ("k", (0.7, 0.5), "red"), ("u", (0.9, 0.2), "red")], [(2,3), (2,4), (1,3), (3, 4), (4, 5), (2,1)])
	save_pdf(img, "mirror", save_mirror)
end

# ╔═╡ 185b0604-3cc6-42d5-b67f-68c1cb985649
mirror_res = let
	code = ein"u,v,w,j,k,vw,vj,vk,ku,wj,wk->juwk"
	mis_solve(code)
end

# ╔═╡ 5b3d225f-a168-4bdb-827f-9bdf85b87e46
compress!(copy(mirror_res))

# ╔═╡ 7a70783a-b9d6-409f-93e7-73443e55d4ca
mirror_res2 = let
	code = ein"u,v,w,j,k,vw,vj,vk,ku,wj,wk->juwk"
	mis_configs(code)
end

# ╔═╡ Cell order:
# ╠═c8c14900-b7f5-11eb-3dbe-4390dd0fcc0a
# ╠═2a37780e-53b6-4a40-9738-5daba20c1001
# ╠═a8ac2b58-86c7-4964-ac10-46c842cbe731
# ╟─e57abcc2-7b48-4e38-a90a-fe6a6bc1eeb7
# ╠═2989b2da-6951-470d-8ef1-cd08628f50f4
# ╠═d4921f23-459a-4436-a526-71951b2d438e
# ╠═80fb581d-c2ab-4d81-97cb-b65bd038a7ce
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
# ╟─3504e107-dcf3-4527-853b-b6a9fe3ced4b
# ╠═e291fac5-ef7a-4308-a9d8-53adcc0bbf11
# ╟─300114b5-d8f7-47bb-a569-2835321f9695
# ╠═6b590957-d846-407a-99fe-b570d5e92107
# ╠═91fb9192-22c6-4135-81c7-efa7ea681722
# ╠═1cb4f4fe-2b5c-4ecf-88b5-3c7fec5d9b87
# ╟─cb1d664e-4efe-4d65-a2be-239c8cd3cfcb
# ╟─30b2968f-413d-4e10-bb80-e71a5596dc5c
# ╟─6766e2fa-795d-47fb-9e95-47471f9cd8e2
# ╠═9bdd69ac-352c-4b88-acbb-5564e57322e5
# ╟─8c97d833-0c37-40cc-9460-58ea5f1d2e0b
# ╠═687a88e1-76ad-4563-82aa-8487f65c441f
# ╠═23fcfe1a-e7ce-4cdc-abcb-651b8f069606
# ╟─c9c09456-767a-4634-bca0-6897b7c5e923
# ╟─70b13ee9-7282-4677-8d27-e7f3f83dc79e
# ╠═a7d0bcd6-a303-4133-9bb1-ba0bb30d294b
# ╟─792d4aaf-b247-48d7-a321-2a3ae73cb988
# ╠═185b0604-3cc6-42d5-b67f-68c1cb985649
# ╠═5b3d225f-a168-4bdb-827f-9bdf85b87e46
# ╠═7a70783a-b9d6-409f-93e7-73443e55d4ca
