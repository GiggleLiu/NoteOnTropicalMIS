using OMEinsum: NestedEinsum
export misb, misv, compress!

function misb(::Type{T}) where T
	res = ones(T, 2, 2)
	res[2, 2] = zero(T)
	return res
end

function misv(::Type{T}, n::Int, val) where T
	res = zeros(T, fill(2, n)...)
	res[1] = one(T)
	res[end] = T(val)
	return res
end

⪰(a, b) = b ⪯ a
⪯(a, b) = (a & b) == a

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

function cross(::Type{T}) where T
	ein"((((a,c),b),d),bd,ac)->abcd"([misv(T, 1, 1.0) for i=1:4]..., misb(T), misb(T))
end

_auto_mistensor(::Type{T}, ix::NTuple{2}) where T = misb(T)
_auto_mistensor(::Type{T}, ix::NTuple{1}) where T = misv(T, 1, 1.0)
function generate_xs!(f, ::Type{T}, code::NestedEinsum, xs) where {T}
    for (ix, arg) in zip(OMEinsum.getixs(code.eins), code.args)
		if arg isa Integer
			xs[arg] = f(T, ix)
		else
        	generate_xs!(f, T, arg, xs)
		end
    end
	return xs
end

function generate_xs!(f, ::Type{T}, code::EinCode, xs) where {T}
    for (i,ix) in enumerate(OMEinsum.getixs(code))
		xs[i] = f(T, ix)
    end
	return xs
end

ninput(::EinCode{ixs}) where ixs = length(ixs)
ninput(ne::Int) = 1
function ninput(ne::NestedEinsum)
    mapreduce(ninput, +, ne.args, init=0)
end

export mis_solve, mis_count
function mis_contract(::Type{T}, code) where {T}
	xs = generate_xs!(_auto_mistensor, T, code, Vector{Any}(undef, ninput(code)))
	code(xs...)
end
mis_solve(code) = mis_contract(TropicalF64, code)
mis_count(code) = mis_contract(CountingTropical{Float64,Float64}, code)

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

