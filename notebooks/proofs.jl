### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ c8c14900-b7f5-11eb-3dbe-4390dd0fcc0a
using Revise, OMEinsum, TropicalNumbers, NoteOnTropicalMIS, Compose, Viznet

# ╔═╡ 3d08e981-b041-4cbf-a808-76532ccc9359
code1 = ein"v,va,vb,vc->abc"

# ╔═╡ 83fb52ab-634e-48e5-8231-6196c240fb1e
mis_solve(code1)

# ╔═╡ 10737ca0-5419-432b-ba10-6e63d6ddddf0
vizeinsum(code1, ['v'=>(0.0, 0.0), 'a'=>(0.0, 1.0), 'b'=>(1.0, 0.0), 'c'=>(1.0, 1.0)], graphsize=5cm, unit=5.0, textoffset=(0.15, 0.0), rescale=0.9)

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
code2 = ein"v,va,vb,vc,a,b,c->abc"

# ╔═╡ 1891feb3-4e3a-4fdd-aca9-5b14fef13c30
res2 = mis_solve(code2)

# ╔═╡ ecd0f2af-338c-432a-9ff8-5c7bdfb410d2
vizeinsum(code2, ['v'=>(0.0, 0.0), 'a'=>(0.0, 1.0), 'b'=>(1.0, 0.0), 'c'=>(1.0, 1.0)], graphsize=5cm, unit=5.0, textoffset=(0.15, 0.0), rescale=0.9)

# ╔═╡ 545961ab-8145-470e-bc62-e7ea7ddacbea
compress!(res2)

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
md"case: `{u₁, u₂} ∉ E and |N₂(v)>1|`"

# ╔═╡ Cell order:
# ╠═c8c14900-b7f5-11eb-3dbe-4390dd0fcc0a
# ╠═3d08e981-b041-4cbf-a808-76532ccc9359
# ╠═83fb52ab-634e-48e5-8231-6196c240fb1e
# ╠═10737ca0-5419-432b-ba10-6e63d6ddddf0
# ╟─57033470-7163-4c10-bf2b-5126f059a1eb
# ╠═327414f7-7e8d-4d74-883c-25e35adf1a9a
# ╠═1891feb3-4e3a-4fdd-aca9-5b14fef13c30
# ╠═ecd0f2af-338c-432a-9ff8-5c7bdfb410d2
# ╠═545961ab-8145-470e-bc62-e7ea7ddacbea
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
