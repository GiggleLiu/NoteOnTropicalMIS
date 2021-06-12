### A Pluto.jl notebook ###
# v0.14.7

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

# ╔═╡ 4e4610ca-99e9-43e5-906b-c4f86e17ecd4
using OMEinsum, Viznet, TropicalNumbers, OMEinsumContractionOrders, Compose, PlutoUI

# ╔═╡ 18124cde-bdf0-11eb-382a-a727d1fbe17a
md"# Tensor network playground"

# ╔═╡ 7be402c5-bc70-49d9-aeb2-2cc4783660cd
md"## Generate a unit-disk graph"

# ╔═╡ d704acbf-12c0-4468-b122-3f25406e2e5a
n = 20;

# ╔═╡ 5527cdcb-dc41-4237-aeb0-f38b1e218729
unit_distance = 0.4;

# ╔═╡ 883443b5-8e7f-49e9-8234-6854f0609c8b
locations = [(rand(), rand()) for _ = 1:n];

# ╔═╡ 6d22f9be-a17d-47f5-8dc8-9a3b75368ed9
let
	xlocs = getindex.(locations, 1)
	ylocs = getindex.(locations, 2)
	xmin, ymin = minimum(xlocs), minimum(ylocs)
	Δx = maximum(xlocs) - xmin
	Δy = maximum(ylocs) - ymin
	Δ = max(Δx, Δy) / 0.8
	xlocs .-= xmin
	ylocs .-= ymin
	xlocs .= xlocs ./ Δ .+ 0.1
	ylocs .= ylocs ./ Δ .+ 0.1
	u = 0.2/Δ
	lattice = UnitDisk([zip(xlocs, ylocs)...], unit_distance/Δ)
	Compose.set_default_graphic_size(8cm, 8cm)
	viz(lattice; line_style=bondstyle(:default, stroke("black"), linewidth(2*u*mm)),
        node_style=nodestyle(:default, r=u*0.15, stroke("black"), fill("white"), linewidth(2*u*mm)),
        text_style=textstyle(:default, fontsize(u*50pt)), labels=1:n)
end

# ╔═╡ 2c23b19d-ddab-4332-b250-9e1356ef1717
md"## Generate the contraction code"

# ╔═╡ f8333fda-c181-47d2-9c32-e5f991c0c489
openlegs = [2,3,1,4]

# ╔═╡ 609b848e-ba20-4dac-85e5-c06bead5fc9c
code = let
	# tensor labels
	ixs = [(i,j) for i=1:n, j=1:n if j>i && sum(abs2, locations[i] .- locations[j]) < unit_distance^2]
	code = EinCode((ixs..., [(i,) for i in 1:n]...), (openlegs...,))
end

# ╔═╡ c7e48886-3cc4-44cf-a9dd-0fa54f05c860
md"""
contraction order optimizer $(@bind order_optimizer Select(["Greedy", "KaHyPar"]))
"""

# ╔═╡ 4480b2b7-4676-4980-a0a9-13831e223e27
let
	if order_optimizer == "Greedy"
		@doc optimize_greedy
	else
		md"""
$(@doc optimize_kahypar)
target space complexity = $(@bind sc_target NumberField(2:32; default=20))
"""
	end
end

# ╔═╡ db4fdd62-7fc9-4427-b072-2b992f59e3fd
optimized_code = if order_optimizer == "Greedy"
	optimize_greedy(code, Dict([i=>2 for i=1:n]))
else
	optimize_kahypar(code, Dict([i=>2 for i=1:n]), sc_target=sc_target, max_group_size=50)
end

# ╔═╡ f2b12677-2613-4b56-96ff-450ea388a929
md"## Generate tensors and do the contraction"

# ╔═╡ 092cb056-2c0e-40af-ba1b-8879460c6dd1
function mis_contract(f, ::Type{T}, code::OMEinsum.NestedEinsum) where T
   	tensors = map(OMEinsum.getixs(OMEinsum.flatten(code))) do ix
	   	@assert length(ix) == 1 || length(ix) == 2
		length(ix) == 1 ? [one(T), convert(T, f(ix[1]))] : [one(T) one(T); one(T) zero(T)]
   	end
   	code(tensors...)
end

# ╔═╡ 5d68b4ad-575f-46fe-a248-0710aa64fd98
element_type = TropicalF64   # the element type for computing

# ╔═╡ 1aba7973-9a14-4715-b801-8438779e2810
xfunc = ix->TropicalF64(1.0) # the function generates `x`, take vertex index as input

# ╔═╡ 42628cd3-5a34-46c6-a1a2-325463fafe45
mis_size = mis_contract(xfunc, element_type, optimized_code)

# ╔═╡ 0d4bb3b5-4eeb-41a9-951b-c05f21c6df55
md"## Notes

1. The `KaHyPar` based contraction order optimizer can generate much better contraction order.
2. Package `TropicalGEMM` provides BLAS level speed for tropical tensor contraction. To use this feature, just install and import the package.
"

# ╔═╡ Cell order:
# ╟─18124cde-bdf0-11eb-382a-a727d1fbe17a
# ╠═4e4610ca-99e9-43e5-906b-c4f86e17ecd4
# ╟─7be402c5-bc70-49d9-aeb2-2cc4783660cd
# ╠═d704acbf-12c0-4468-b122-3f25406e2e5a
# ╠═5527cdcb-dc41-4237-aeb0-f38b1e218729
# ╠═883443b5-8e7f-49e9-8234-6854f0609c8b
# ╟─6d22f9be-a17d-47f5-8dc8-9a3b75368ed9
# ╟─2c23b19d-ddab-4332-b250-9e1356ef1717
# ╠═f8333fda-c181-47d2-9c32-e5f991c0c489
# ╠═609b848e-ba20-4dac-85e5-c06bead5fc9c
# ╟─c7e48886-3cc4-44cf-a9dd-0fa54f05c860
# ╟─4480b2b7-4676-4980-a0a9-13831e223e27
# ╠═db4fdd62-7fc9-4427-b072-2b992f59e3fd
# ╟─f2b12677-2613-4b56-96ff-450ea388a929
# ╠═092cb056-2c0e-40af-ba1b-8879460c6dd1
# ╠═5d68b4ad-575f-46fe-a248-0710aa64fd98
# ╠═1aba7973-9a14-4715-b801-8438779e2810
# ╠═42628cd3-5a34-46c6-a1a2-325463fafe45
# ╟─0d4bb3b5-4eeb-41a9-951b-c05f21c6df55
