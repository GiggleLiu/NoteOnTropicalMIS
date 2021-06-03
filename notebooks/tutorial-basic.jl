### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 4e4610ca-99e9-43e5-906b-c4f86e17ecd4
using OMEinsum, LightGraphs, TropicalNumbers

# ╔═╡ f67a7abb-d96c-4167-962f-b9afa23300cd
using Polynomials

# ╔═╡ 0fd7b695-fd29-4570-b33f-a9261565b06b
using ForwardDiff

# ╔═╡ 496e51c1-6189-4cd2-8f82-8f0f6767f3b6
using Zygote

# ╔═╡ 18124cde-bdf0-11eb-382a-a727d1fbe17a
md"# Solving MIS with generic einsum

*A simple introduction that even a gorilla can understand*

Jinguo Liu"

# ╔═╡ 7be402c5-bc70-49d9-aeb2-2cc4783660cd
md"## Generate a graph"

# ╔═╡ 21165778-a9c7-452d-8330-451550c4176f
n, k = 16, 3

# ╔═╡ 6329fc8b-2956-4b90-9323-b21704e3c08f
g = LightGraphs.random_regular_graph(n, k)

# ╔═╡ b326c7f2-2caf-4cbd-9749-76307a984ace
md"## Generate the einsum code for the graph"

# ╔═╡ 47c852b5-9cc4-4691-8264-4b314359efc6
md"`einsum` is the so called Einstein's notation, it loops over indices, multiplies elements, and accumulates the results to the output tensor. e.g. the following statement does the matrix multiplication"

# ╔═╡ aa743568-a195-4f0f-aff7-8a8f77c125a6
ein"ij,jk->ik"(randn(2,2), randn(2,2))

# ╔═╡ 8f1198cd-b342-4bad-8158-576e2f01ca21
md"The following statement does the trace"

# ╔═╡ f043399f-6ef8-41bf-bdb5-8899109faad1
ein"ii->"(randn(2,2))

# ╔═╡ 9e3adcec-c012-43c5-8527-04e4e0674462
md"The following statement does the MPS contraction"

# ╔═╡ 7e8076e2-44a6-4e38-9ffa-d7684fd6a5cd
ein"(((αj,jβk),kγl),lδm),mξ->αβγδξ"(randn(2,2), randn(2,2,2), randn(2,2,2), randn(2,2,2), randn(2,2))

# ╔═╡ 12ab7ae5-c703-49d2-bb89-3bc4a0dced59
md"Here, we use `()` to specify the contraction orders, otherwise, it will use the naive looping (can be very slow). In the following, we contruct the contraction for the above graph."

# ╔═╡ 609b848e-ba20-4dac-85e5-c06bead5fc9c
ixs = [minmax(e.src,e.dst) for e in LightGraphs.edges(g)]

# ╔═╡ 099bc38d-6811-4fc6-83ec-96f3a4f54257
md"`ixs` is a vector of labels for bond tensors, the number of bond tensor is equivalent to the number of edges: $(ne(g))"

# ╔═╡ 8d1417eb-c062-44bf-8868-a9b5572d3035
code = EinCode((ixs..., [(i,) for i in LightGraphs.vertices(g)]...), ())

# ╔═╡ 9a36bfd4-5f57-4438-ae37-25bfa0b8640b
md"The `EinCode` contructor takes two input arguments. The first argument is a tuple of tuple, each tuple contains labels for one of the input tensors, while the second argument is a tuple of output tensor labels. One can not use this code directly for computing, because it is just naive looping and can be very slow. To optimize the contraction order, one can call the greedy optimizer in OMEinsum."

# ╔═╡ db4fdd62-7fc9-4427-b072-2b992f59e3fd
optimized_code = optimize_greedy(code, Dict([i=>2 for i=1:nv(g)]))

# ╔═╡ 4480b2b7-4676-4980-a0a9-13831e223e27
@doc optimize_greedy

# ╔═╡ e9e9f83b-5c1f-4cbf-883d-6e1f345e6355
md"Here, the size of each bond is 2."

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

# ╔═╡ 42628cd3-5a34-46c6-a1a2-325463fafe45
mis_size = mis_contract(ix->TropicalF64(1.0), TropicalF64, optimized_code)

# ╔═╡ 4fe151d5-cc5a-4447-9ae6-bb6a1fbbddf8
mis_size_with_degeneracy = mis_contract(ix->CountingTropicalF64(1.0), CountingTropicalF64, optimized_code)

# ╔═╡ 112f8104-88b2-4968-8c7a-2d96ef6a52b8
md"For independence polynomial"

# ╔═╡ 3ac40524-1f0f-4105-83b9-561a0810d475
OMEinsum.asarray(x, ::AbstractArray) = fill(x)

# ╔═╡ 87bcf018-3c56-427c-a8c3-24e600062f06
independence_polynomial = mis_contract(ix->Polynomial([0,1.0]), Polynomial, optimized_code)

# ╔═╡ ddacba43-20e2-4575-9ba1-35e26371f723
gadget = mis_contract(ix->TropicalF64(1), TropicalF64, optimize_greedy(ein"a,b,c,d,e,f,g,h,i,j,k,ae,ag,be,bf,ci,ck,dj,dk,ef,eg,fg,fh,gh,hi,hj,ij,ik,jk->abcd", Dict([(x=>2) for x in "abcdefghijk"])))

# ╔═╡ 7fb0099b-511b-496b-a411-839f72b62556
content.(gadget) .- 3

# ╔═╡ 03849452-3b16-43d9-b42b-80136b024227
md"### Optimal configuration"

# ╔═╡ 0b4509ba-a05a-43ab-900e-7586a0429223
md"One can use the `ConfigTropical` to compute optimal configurations. It has two fields, a tropical number and a static bit string vector. The bit string vector can be constructed with the following functions."

# ╔═╡ eb83a776-9629-4ae9-98ad-943426713a70
TropicalNumbers.staticfalses(StaticBitVector{11,1})  # 0s

# ╔═╡ f8a22082-35d7-4045-a210-b0ee6e02ab82
TropicalNumbers.statictrues(StaticBitVector{11,1})   # 1s

# ╔═╡ 0b99d453-220e-4b35-8b20-86dbdf687e07
TropicalNumbers.onehot(StaticBitVector{11,1}, 3)     # one hot vector

# ╔═╡ 327bc673-3490-4af0-b236-1af0519244ad
TropicalNumbers.onehot(StaticBitVector{11,1}, 3) | TropicalNumbers.onehot(StaticBitVector{11,1}, 5)     # bit-wise or

# ╔═╡ 74a47e93-00d9-4926-b5f1-2cf8a47e9503
mis_contract(ix->ConfigTropical(1.0, TropicalNumbers.onehot(StaticBitVector{11,1}, ix-'a'+1)), ConfigTropical{Float64,11,1}, optimize_greedy(ein"a,b,c,d,e,f,g,h,i,j,k,ae,ag,be,bf,ci,ck,dj,dk,ef,eg,fg,fh,gh,hi,hj,ij,ik,jk->abcd", Dict([(x=>2) for x in "abcdefghijk"])))

# ╔═╡ 7c0aa88b-b451-4210-aa30-d24c1e6f78ff
md"## Autodiff"

# ╔═╡ 6bceb180-095b-41e0-b6ca-e9c02e4117d6
md"#### Forward mode autodiff"

# ╔═╡ 4f2e4d03-2f57-4aa8-8875-0aa4f8001f14
tropical_loss(x) = mis_contract(_->Tropical(x), Tropical{eltype(x)}, optimized_code)[].n

# ╔═╡ 775f36e6-3663-4b68-9dc0-c3bb1b63f184
ForwardDiff.gradient(x->tropical_loss(x[]), [1.0])

# ╔═╡ 90c8611b-6c44-4398-955a-9e7af6d096fa
md"#### Reverse mode autodiff"

# ╔═╡ 4b1aae0b-9dfe-4978-9de8-3a3a232b2413
regular_loss(xs...) = ein"((ij,jk),klb),l->"(xs...)[]

# ╔═╡ e0536815-0e5e-472d-a9a4-e21f04fa09c8
Zygote.gradient(regular_loss, (randn(2,2), randn(2,2), randn(2,2,2), randn(2))...)

# ╔═╡ 0d4bb3b5-4eeb-41a9-951b-c05f21c6df55
md"## Notes

1. The `KaHyPar` based contraction order optimizer can generate much better contraction order.
2. Package `TropicalGEMM` provides BLAS level speed for tropical tensor contraction. To use this feature, just install and import the package.
"

# ╔═╡ Cell order:
# ╟─18124cde-bdf0-11eb-382a-a727d1fbe17a
# ╠═4e4610ca-99e9-43e5-906b-c4f86e17ecd4
# ╟─7be402c5-bc70-49d9-aeb2-2cc4783660cd
# ╠═21165778-a9c7-452d-8330-451550c4176f
# ╠═6329fc8b-2956-4b90-9323-b21704e3c08f
# ╟─b326c7f2-2caf-4cbd-9749-76307a984ace
# ╟─47c852b5-9cc4-4691-8264-4b314359efc6
# ╠═aa743568-a195-4f0f-aff7-8a8f77c125a6
# ╟─8f1198cd-b342-4bad-8158-576e2f01ca21
# ╠═f043399f-6ef8-41bf-bdb5-8899109faad1
# ╟─9e3adcec-c012-43c5-8527-04e4e0674462
# ╠═7e8076e2-44a6-4e38-9ffa-d7684fd6a5cd
# ╟─12ab7ae5-c703-49d2-bb89-3bc4a0dced59
# ╠═609b848e-ba20-4dac-85e5-c06bead5fc9c
# ╟─099bc38d-6811-4fc6-83ec-96f3a4f54257
# ╠═8d1417eb-c062-44bf-8868-a9b5572d3035
# ╟─9a36bfd4-5f57-4438-ae37-25bfa0b8640b
# ╠═db4fdd62-7fc9-4427-b072-2b992f59e3fd
# ╟─4480b2b7-4676-4980-a0a9-13831e223e27
# ╟─e9e9f83b-5c1f-4cbf-883d-6e1f345e6355
# ╟─f2b12677-2613-4b56-96ff-450ea388a929
# ╠═092cb056-2c0e-40af-ba1b-8879460c6dd1
# ╠═42628cd3-5a34-46c6-a1a2-325463fafe45
# ╠═4fe151d5-cc5a-4447-9ae6-bb6a1fbbddf8
# ╟─112f8104-88b2-4968-8c7a-2d96ef6a52b8
# ╠═f67a7abb-d96c-4167-962f-b9afa23300cd
# ╠═3ac40524-1f0f-4105-83b9-561a0810d475
# ╠═87bcf018-3c56-427c-a8c3-24e600062f06
# ╠═ddacba43-20e2-4575-9ba1-35e26371f723
# ╠═7fb0099b-511b-496b-a411-839f72b62556
# ╟─03849452-3b16-43d9-b42b-80136b024227
# ╟─0b4509ba-a05a-43ab-900e-7586a0429223
# ╠═eb83a776-9629-4ae9-98ad-943426713a70
# ╠═f8a22082-35d7-4045-a210-b0ee6e02ab82
# ╠═0b99d453-220e-4b35-8b20-86dbdf687e07
# ╠═327bc673-3490-4af0-b236-1af0519244ad
# ╠═74a47e93-00d9-4926-b5f1-2cf8a47e9503
# ╟─7c0aa88b-b451-4210-aa30-d24c1e6f78ff
# ╟─6bceb180-095b-41e0-b6ca-e9c02e4117d6
# ╠═4f2e4d03-2f57-4aa8-8875-0aa4f8001f14
# ╠═0fd7b695-fd29-4570-b33f-a9261565b06b
# ╠═775f36e6-3663-4b68-9dc0-c3bb1b63f184
# ╟─90c8611b-6c44-4398-955a-9e7af6d096fa
# ╠═4b1aae0b-9dfe-4978-9de8-3a3a232b2413
# ╠═496e51c1-6189-4cd2-8f82-8f0f6767f3b6
# ╠═e0536815-0e5e-472d-a9a4-e21f04fa09c8
# ╟─0d4bb3b5-4eeb-41a9-951b-c05f21c6df55
