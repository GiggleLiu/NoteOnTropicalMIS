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
begin
	using Revise
	using NoteOnTropicalMIS
	using NoteOnTropicalMIS.OMEinsum, NoteOnTropicalMIS.TropicalNumbers, NoteOnTropicalMIS.Polynomials
	using OMEinsumContractionOrders, PlutoUI
end

# ╔═╡ 18124cde-bdf0-11eb-382a-a727d1fbe17a
md"# Tensor network playground
-- Random Regular graphs"

# ╔═╡ 7be402c5-bc70-49d9-aeb2-2cc4783660cd
md"## Generate a random-regular graph"

# ╔═╡ d704acbf-12c0-4468-b122-3f25406e2e5a
n = 80;

# ╔═╡ 5527cdcb-dc41-4237-aeb0-f38b1e218729
k = 3;

# ╔═╡ aa4fec11-2e61-4dd1-9484-3aca470f07bf
code = random_regular_eincode(n, k);

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
end;

# ╔═╡ 4df90029-1950-4701-a5ea-3b432992de61
let
	tc, sc = OMEinsum.timespace_complexity(optimized_code, uniformsize(code, 2))
	Text("time complexity = 2^$(round(tc; sigdigits=3)), space complexity = 2^$(sc)")
end

# ╔═╡ f2b12677-2613-4b56-96ff-450ea388a929
md"## Generate tensors and do the contraction"

# ╔═╡ 092cb056-2c0e-40af-ba1b-8879460c6dd1
function mis_contract(f, ::Type{T}, code::OMEinsum.NestedEinsum) where T
   	tensors = map(OMEinsum.getixs(Iterators.flatten(code))) do ix
	   	@assert length(ix) == 1 || length(ix) == 2
		length(ix) == 1 ? [one(T), convert(T, f(ix[1]))] : [one(T) one(T); one(T) zero(T)]
   	end
   	code(tensors...)
end

# ╔═╡ 7461f34d-39d0-4d99-a717-9f76815f1273
md"""
$(@bind element_type Select(["Select Task"=>"Select Task", "Regular"=>"number of independent sets", "Tropical"=>"maximum independent set size", "Polynomial"=>"Independence polynomial", "Tropical-Counting"=>"maximum independent set degeneracy", "Tropical-Config"=>"single optimal configuration", "Tropical-Config-All"=>"all optimal configurations"]))
"""

# ╔═╡ 42628cd3-5a34-46c6-a1a2-325463fafe45
if element_type !== "Select Task" 
	let
		T = Dict(
			"Regular"=>Float64,
			"Tropical"=>TropicalF64,
			"Polynomial"=>Polynomial,
			"Tropical-Counting"=>CountingTropicalF64,
			"Tropical-Config"=>ConfigTropical{Float64,n,TropicalNumbers._nints(n)},
			"Tropical-Config-All"=>CountingTropical{Float64,ConfigEnumerator{n,TropicalNumbers._nints(n)}},
			)[element_type]

		xfunc = if element_type == "Tropical-Config"
			c = TropicalNumbers._nints(n)
			ix->T(1, TropicalNumbers.onehot(StaticBitVector{n,c}, ix))
		elseif element_type == "Tropical-Config-All"
			c = TropicalNumbers._nints(n)
			ix->T(1, ConfigEnumerator([TropicalNumbers.onehot(StaticBitVector{n,c}, ix)]))
		elseif element_type == "Polynomial"
			ix->Polynomial([0.0, 1.0])
		else
			ix->T(1)
		end

		mis_contract(xfunc, T, optimized_code)[]
	end
end

# ╔═╡ 0d4bb3b5-4eeb-41a9-951b-c05f21c6df55
md"## Notes
Package `TropicalGEMM` provides BLAS level speed for tropical tensor contraction. To use this feature, just install and import the package.
"

# ╔═╡ Cell order:
# ╟─18124cde-bdf0-11eb-382a-a727d1fbe17a
# ╠═4e4610ca-99e9-43e5-906b-c4f86e17ecd4
# ╟─7be402c5-bc70-49d9-aeb2-2cc4783660cd
# ╠═d704acbf-12c0-4468-b122-3f25406e2e5a
# ╠═5527cdcb-dc41-4237-aeb0-f38b1e218729
# ╠═aa4fec11-2e61-4dd1-9484-3aca470f07bf
# ╟─c7e48886-3cc4-44cf-a9dd-0fa54f05c860
# ╟─4480b2b7-4676-4980-a0a9-13831e223e27
# ╠═db4fdd62-7fc9-4427-b072-2b992f59e3fd
# ╟─4df90029-1950-4701-a5ea-3b432992de61
# ╟─f2b12677-2613-4b56-96ff-450ea388a929
# ╠═092cb056-2c0e-40af-ba1b-8879460c6dd1
# ╟─7461f34d-39d0-4d99-a717-9f76815f1273
# ╠═42628cd3-5a34-46c6-a1a2-325463fafe45
# ╟─0d4bb3b5-4eeb-41a9-951b-c05f21c6df55
