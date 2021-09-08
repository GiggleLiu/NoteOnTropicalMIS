### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ ad82840f-f908-49f8-9bc0-091b0ceaafad
using Pkg; Pkg.activate()

# ╔═╡ 50862912-0ed6-11ec-3db0-01dce2a8a110
using OMEinsum, SymEngine

# ╔═╡ 67f68723-cb66-44c1-8c04-83eacf6eaaf0
using Test

# ╔═╡ 809e5737-3dd3-4d1f-817a-76b9d61dbd25
using Roots

# ╔═╡ 0077b8a4-207a-4dff-b4f1-4ece3055e9f5
using Zygote

# ╔═╡ 611e081b-f99c-4304-8e15-276523345ad3
using Optim

# ╔═╡ 5a96edc0-b8a1-4da8-885a-5e1d3e7a74fe
using StochasticOptimizers

# ╔═╡ ffb21993-9d10-4a66-b0a7-453670c80503
b, g, ξ, η, z = Basic.((:b, :g, :ξ, :η, :z))

# ╔═╡ b0c4c062-8a59-45a6-aecb-f122a4c5e878
A = [1, b]

# ╔═╡ a1b3cde7-61a4-4ea1-95ff-d24299afd6ab
F = [1 g; g 0]

# ╔═╡ f7c88204-8a60-4dca-8914-9ce998fe4236
w = let
	x = zeros(Basic, 2,2,2)
	x[1,1,1] = one(Basic)
	x[1,1,2] = x[2,1,1] = x[1,2,1] = Basic(z)^(1//6)
	x
end

# ╔═╡ 62e6be7b-4705-4255-a3ef-5833f268f24c
md"4a"

# ╔═╡ 9a277a55-feb4-4ec9-b8b0-06ee7cdfee4f
ein"(ab,(b,b)),ba->a"(F,A,A,F)

# ╔═╡ 8fca3104-230b-4ad2-bcfa-eb652b7a1e6b
ein"(a,a),(a,a)->a"(A,A,A,A)*ξ

# ╔═╡ 85af48ad-a493-4e10-b489-42ce6470612a
md"4b"

# ╔═╡ 8a5f7449-83e8-48ad-96bc-08b9b97df72a
ein"(abc,ac),(c,cb)->ab"(w,F,A,F)

# ╔═╡ 20e37afb-bc0c-49fd-adc8-7cba7455d06c
ein"(a,ab),b->ab"(A,F,A) * sqrt(η)

# ╔═╡ c54e3d4b-30e5-4aae-a4d2-8e2a6fff7350
rand_hermitian(::Type{T}, D::Int) where T = (x=randn(T, D,D); x+x')

# ╔═╡ 0e940c95-7837-4e1a-9a55-5919b8203829
D = 2

# ╔═╡ afd96431-154f-41c4-9172-78f1b87320d5
function rand_Atensor(::Type{T}, D::Int) where T
	A = zeros(T, D, D, 2)
	A[:,:,1] = rand_hermitian(T, D)
	A[:,:,2] = rand_hermitian(T, D)
	A
end;

# ╔═╡ 5e6dbfe5-a220-4a30-9993-d58bf50fb366
function rand_Ftensor(::Type{T}, D::Int) where T
	F = zeros(T, D, D, 2, 2)
	F[:,:,1,1] = rand_hermitian(T, D)
	F[:,:,2,2] = zeros(T, D, D)
	F[:,:,1,2] = randn(Float64, D, D)
	F[:,:,2,1] = F[:,:,1,2]
	F
end;

# ╔═╡ 2cb035c0-46bb-4a6a-aff8-a3a4e9ce527e
function wtensor(z::T) where T
	x = zeros(T,2,2,2)
	x[1,1,1] = one(T)
	x[1,1,2] = x[2,1,1] = x[1,2,1] = z^(1//6)
	x
end

# ╔═╡ dd13a7ba-aa18-49ff-ba3b-ca8977606885
function getρ(A::AbstractArray)
	A6 = ein"((((ija,jka),kla),lma),mna),noa->ioa"(A,A,A,A,A,A)
	Ap = ein"ii->"(A6[:,:,1])
	An = ein"ii->"(A6[:,:,2])
	return An[]/(An+Ap)[]
end

# ╔═╡ 8025bb4c-b1c5-432f-8ed4-859bd00b6218
let
	x = 0.23960969574206287
	x/(1+2x)
end

# ╔═╡ 7cae9004-9aeb-491b-8b46-e49ae6b40c29
function loss(A, F, w, sqrtη, ξ)
	x1 = sum(abs.(ein"(ijab,(jkb,klb)),lmba->ima"(F,A,A,F) .- ein"(ija,jka),(kla,lma)->ima"(A,A,A,A) .* ξ))
	x2 = sum(abs.(ein"(abc,ijac),(jkc,klcb)->ilab"(w,F,A,F) .- ein"(ija,jkab),klb->ilab"(A,F,A) .* sqrtη))
	x1 + x2
end

# ╔═╡ d051399d-3b99-494c-8c14-5de9e68b5a5c
let
	@testset "rank 1 loss" begin
		x = Roots.find_zero(x->(1-x)^5*(1-x^2)-x, 0.5)
		b = (x/(1+x))^(1//6)
		z = x/(1-x)^5/(1-x^2)
		g = sqrt(b^4/(1-b^6))
		sqrtη = (1+z^(1//6)*g*b*g)
		ξ = 1+g*b^2*g
		A2 = rand_Atensor(Float64,1)
		F2 = rand_Ftensor(Float64,1)
		w2 = wtensor(z)
		A2[:,:,1] .= 1
		A2[:,:,2] .= b
		F2[:,:,2,1] .= g
		F2[:,:,1,2] .= g
		F2[:,:,1,1] .= 1
		@test isapprox(loss(A2, F2, wtensor(z), sqrtη, ξ), 0; atol=1e-8)
		ρ = getρ(A2)
		@test isapprox(ρ, 0.161984; atol=1e-5)
		@test isapprox(log(sqrtη^2/ξ), 0.333050; atol=1e-5)
	end
end

# ╔═╡ 3dff825f-0722-4e1f-808c-d0a8eb7a1532
gs = let
	z = 1.0
	D = 2
	A2 = rand_Atensor(Float64,D)
	F2 = rand_Ftensor(Float64,D)
	w2 = wtensor(z)
	sqrtη = rand()
	ξ = rand()
	loss(A2, F2, wtensor(z), sqrtη, ξ)
	Zygote.gradient(loss, A2, F2, wtensor(z), sqrtη, ξ)
end

# ╔═╡ d407426b-b0ef-4d21-8a5f-d2976be7f1d6
function symmetrize!(A::AbstractMatrix)
	A .= (A .+ A') ./ 2
	return A
end

# ╔═╡ 0db38e20-23e6-498f-b139-e89949763f00
function pack!(y, gA2, gF2, gsqrtη, gξ)
	n = size(gA2, 1)
	y[1:n^2*2] .= vec(gA2)
	y[n^2*2+1:n^2*6] .= vec(gF2)
	y[n^2*6+1] = gsqrtη
	y[n^2*6+2] = gξ
	return y
end

# ╔═╡ 701934b7-22e9-4c2b-9d39-b3e0a15bb263
function unpack(x)
	n = Int(sqrt((length(x) - 2)÷6))
	reshape(x[1:n^2*2],n,n,2), reshape(x[n^2*2+1:n^2*6],n,n,2,2), x[n^2*6+1], x[n^2*6+2]
end

# ╔═╡ e777422f-339a-4217-a216-96ce6f16f920
function g!(z)
	function (y, x)
		A2, F2, sqrtη, ξ = unpack(x)
		gA2, gF2, gw, gsqrtη, gξ = Zygote.gradient(loss, A2, F2, wtensor(z), sqrtη, ξ)
		gA2[1] = 0
		gF2[1] = 0
		symmetrize!(view(gA2,:,:,1))
		symmetrize!(view(gA2,:,:,2))
		symmetrize!(view(gF2,:,:,1,1))
		gf = view(gF2,:,:,2,1) .+ view(gF2,:,:,1,2)'
		view(gF2,:,:,1,2) .= gf'
		view(gF2,:,:,2,1) .= gf
		pack!(y, gA2, gF2, gsqrtη, gξ)
	end
end

# ╔═╡ 288bb3a3-3049-4034-9260-7e57c041ad69
function lossf(z)
	function (x)
		A2, F2, η, ξ = unpack(x)
		#A2[1] = 1
		#F2[1] = 1
		#symmetrize!(view(A2,:,:,1))
		#symmetrize!(view(A2,:,:,2))
		#symmetrize!(view(F2,:,:,1,1))
		#gf = view(F2,:,:,2,1) .+ view(F2,:,:,1,2)'
		#view(F2,:,:,1,2) .= gf'
		#view(F2,:,:,2,1) .= gf
		@show loss(A2, F2, wtensor(z), η, ξ)
	end
end

# ╔═╡ ff8fd5dc-e176-42c4-9ac8-978625c2f556
res = let
	z=1.0
	n=2
	N = n^2*6+2
	x0 = rand(N)
	x0[1] = 1.0
	x0[2n^2+1] = 1.0
	#Optim.optimize(lossf(z), g!(z), x0, LBFGS(), Optim.Options(iterations=200))
end

# ╔═╡ f3fbd242-2b72-4c77-93b7-52e3b504ae75
getρ(unpack(res.minimizer)[1])

# ╔═╡ b7d70407-8b1a-489e-82d0-62293456b3fc
unpack(res.minimizer)[1]

# ╔═╡ cbcb46b5-bcde-48eb-985f-72115029b89d
res_adam = let
	z=1.0
	n=2
	N = n^2*6+2
	A2 = rand_Atensor(Float64,D)
	F2 = rand_Ftensor(Float64,D)
	A2[1] = 1
	F2[1] = 1
	symmetrize!(view(A2,:,:,1))
	symmetrize!(view(A2,:,:,2))
	symmetrize!(view(F2,:,:,1,1))
	gf = view(F2,:,:,2,1) .+ view(F2,:,:,1,2)'
	view(F2,:,:,1,2) .= gf'
	view(F2,:,:,2,1) .= gf
	x0 = pack!(zeros(N), A2,F2,rand(),rand())
	local res
	for (k,it) in enumerate(adam(lossf(z), x->g!(z)(zero(x), x), x0))
		@show lossf(z)(StochasticOptimizers.minimizer(it))
		if k>20000
			res = StochasticOptimizers.minimizer(it)
			break
		else
			unpack(StochasticOptimizers.minimizer(it))[1]
		end
	end
	res
end

# ╔═╡ b5685983-a93c-4f42-9ce5-d71f28ede9c2
unpack(res_adam)[1]

# ╔═╡ a4a7ac30-8a85-4387-bfc1-606b7446e73f
getρ(unpack(res_adam)[1])

# ╔═╡ 191b9fd1-e099-4ff9-bf5c-53ab03c12a2b
let
	sqrtη, ξ = unpack(res_adam)[3:4]
	sqrtη^2/ξ
end

# ╔═╡ Cell order:
# ╠═ad82840f-f908-49f8-9bc0-091b0ceaafad
# ╠═50862912-0ed6-11ec-3db0-01dce2a8a110
# ╠═ffb21993-9d10-4a66-b0a7-453670c80503
# ╠═b0c4c062-8a59-45a6-aecb-f122a4c5e878
# ╠═a1b3cde7-61a4-4ea1-95ff-d24299afd6ab
# ╠═f7c88204-8a60-4dca-8914-9ce998fe4236
# ╟─62e6be7b-4705-4255-a3ef-5833f268f24c
# ╠═9a277a55-feb4-4ec9-b8b0-06ee7cdfee4f
# ╠═8fca3104-230b-4ad2-bcfa-eb652b7a1e6b
# ╟─85af48ad-a493-4e10-b489-42ce6470612a
# ╠═8a5f7449-83e8-48ad-96bc-08b9b97df72a
# ╠═20e37afb-bc0c-49fd-adc8-7cba7455d06c
# ╠═c54e3d4b-30e5-4aae-a4d2-8e2a6fff7350
# ╠═0e940c95-7837-4e1a-9a55-5919b8203829
# ╠═afd96431-154f-41c4-9172-78f1b87320d5
# ╠═5e6dbfe5-a220-4a30-9993-d58bf50fb366
# ╠═2cb035c0-46bb-4a6a-aff8-a3a4e9ce527e
# ╠═67f68723-cb66-44c1-8c04-83eacf6eaaf0
# ╠═809e5737-3dd3-4d1f-817a-76b9d61dbd25
# ╠═dd13a7ba-aa18-49ff-ba3b-ca8977606885
# ╠═d051399d-3b99-494c-8c14-5de9e68b5a5c
# ╠═8025bb4c-b1c5-432f-8ed4-859bd00b6218
# ╠═0077b8a4-207a-4dff-b4f1-4ece3055e9f5
# ╠═7cae9004-9aeb-491b-8b46-e49ae6b40c29
# ╠═3dff825f-0722-4e1f-808c-d0a8eb7a1532
# ╠═d407426b-b0ef-4d21-8a5f-d2976be7f1d6
# ╠═e777422f-339a-4217-a216-96ce6f16f920
# ╠═288bb3a3-3049-4034-9260-7e57c041ad69
# ╠═611e081b-f99c-4304-8e15-276523345ad3
# ╠═0db38e20-23e6-498f-b139-e89949763f00
# ╠═701934b7-22e9-4c2b-9d39-b3e0a15bb263
# ╠═ff8fd5dc-e176-42c4-9ac8-978625c2f556
# ╠═f3fbd242-2b72-4c77-93b7-52e3b504ae75
# ╠═b7d70407-8b1a-489e-82d0-62293456b3fc
# ╠═5a96edc0-b8a1-4da8-885a-5e1d3e7a74fe
# ╠═cbcb46b5-bcde-48eb-985f-72115029b89d
# ╠═b5685983-a93c-4f42-9ce5-d71f28ede9c2
# ╠═a4a7ac30-8a85-4387-bfc1-606b7446e73f
# ╠═191b9fd1-e099-4ff9-bf5c-53ab03c12a2b
