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

function widget1(::Type{T}, x) where T
	code = ein"a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,ae,bf,cg,dh,ei,ej,ek,el,fi,fm,gm,gn,go,gp,hl,hp,ij,im,in,jk,jm,jn,jo,kl,kn,ko,kp,lo,lp,mn,no,op->abcd"
    code = greedy_order(code)
    code([v(T, 1, xi) for xi in x]..., [b(T) for i=1:32]...)
end

