using OMEinsum: NestedEinsum
export misb, misv

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

function generate_xs!(::Type{T}, codeInt, out, xs) where {T, ixs, iy}
    if length(ix) == 1
        push!(xs, misv(T, 1, 1.0))
    elseif length(ix)==2
        push!(xs, misb(T))
    else
        error("")
    end
end

function generate_xs!(::Type{T}, code::NestedEinsum, out, xs) where {T}
    for (ix, arg) in zip(OMEinsum.getixs(code.eins), code.args)
        generate_xs!(T, arg, ix, xs)
    end
end

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

