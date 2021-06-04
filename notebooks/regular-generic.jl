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
	using NoteOnTropicalMIS.OMEinsum, NoteOnTropicalMIS.TropicalNumbers, NoteOnTropicalMIS.Polynomials, NoteOnTropicalMIS.TropicalGEMM
	using OMEinsumContractionOrders, PlutoUI
end

# ╔═╡ c5db88b3-f294-458f-aaf4-f4c50abc8e4e
using ProfileSVG

# ╔═╡ 525ecde1-d55c-4bde-833e-e5a2bd591165
using TensorOperations, LinearAlgebra

# ╔═╡ 18124cde-bdf0-11eb-382a-a727d1fbe17a
md"# Tensor network playground
-- Random Regular graphs"

# ╔═╡ 7be402c5-bc70-49d9-aeb2-2cc4783660cd
md"## Generate a random-regular graph"

# ╔═╡ d704acbf-12c0-4468-b122-3f25406e2e5a
n = 200;

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
	optimize_kahypar(code, Dict([i=>2 for i=1:n]), sc_target=sc_target, max_group_size=40)
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
$(@bind element_type Select(["Select Task"=>"Select Task", "Regular"=>"number of independent sets", "Tropical"=>"maximum independent set size", "Polynomial"=>"Independence polynomial", "Polynomial-finitefield"=>"Independence polynomial (finite field)", "Polynomial-fft"=>"Independence polynomial (fft)", "Tropical-Counting"=>"maximum independent set degeneracy", "Tropical-Config"=>"single optimal configuration", "Tropical-Config-All"=>"all optimal configurations"]))
"""

# ╔═╡ b80080d6-a178-4b4a-9897-4a17ddee5b50
let
	_r = @bind r NumberField(0.01:0.01:5.0; default=3.0)
	if element_type == "Polynomial-fft"
		md"r = $_r"
	end
end

# ╔═╡ 42628cd3-5a34-46c6-a1a2-325463fafe45
let
	element_dict = Dict(
		"Regular"=>Float64,
		"Tropical"=>TropicalF64,
		"Polynomial"=>Polynomial,
		"Tropical-Counting"=>CountingTropicalF64,
		"Tropical-Config"=>ConfigTropical{Float64,n,TropicalNumbers._nints(n)},
		"Tropical-Config-All"=>CountingTropical{Float64,ConfigEnumerator{n,TropicalNumbers._nints(n)}},
		)
	if haskey(element_dict, element_type)
		T = element_dict[element_type]

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
	elseif element_type == "Polynomial-finitefield"
		independence_polynomial(Val(:finitefield), optimized_code)
	elseif element_type == "Polynomial-fft"
		independence_polynomial(Val(:fft), optimized_code; r=r)
	end
end

# ╔═╡ d792602d-3f7d-4812-bbb2-6ef16560563a
function LinearAlgebra.permutedims!(C::Array{T,N}, A::StridedArray{T,N}, perm) where {T,N}
    TensorOperations.tensorcopy!(A, ntuple(identity,N), C, perm)
end

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
# ╟─b80080d6-a178-4b4a-9897-4a17ddee5b50
# ╠═42628cd3-5a34-46c6-a1a2-325463fafe45
# ╠═c5db88b3-f294-458f-aaf4-f4c50abc8e4e
# ╠═525ecde1-d55c-4bde-833e-e5a2bd591165
# ╠═d792602d-3f7d-4812-bbb2-6ef16560563a
