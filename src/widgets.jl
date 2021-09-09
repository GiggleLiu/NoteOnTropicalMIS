using GraphTensorNetworks: misv, misb
export widget16, widget12, unitdisk_edges, uniformsize, compress!

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

function cross(x::T) where T
	ein"((((a,c),b),d),bd,ac)->abcd"([misv(x) for i=1:4]..., misb(T), misb(T))
end

struct UnitdiskWidget{LT}
    nodes::Vector{Pair{LT,Tuple{Float64,Float64}}}
    radius::Float64
end

function get_edges(ud::UnitdiskWidget)
	edges = Tuple{Int,Int}[]
	for (i, p) in enumerate(ud.nodes)
		for (j,p2) in enumerate(ud.nodes)
			if i<j && sqrt(sum(abs2, p2.second .- p.second)) < ud.radius
				push!(edges, (i,j))
			end
		end
	end
	edges
end

function widget16()
	a = 0.12
	ymid = xmid = 0.5
	X = 0.45
	Y = 0.33
	D = 0.15
	y = [ymid-Y, ymid-Y+D, ymid-a/2, ymid+a/2, ymid+Y-D, ymid+Y]
	x = [xmid-X, xmid-X+D, xmid-1.5a, xmid-a/2, xmid+a/2, xmid+1.5a, xmid+X-D, xmid+X]
	xmin, xmax, ymin, ymax = x[1], x[end], y[1], y[end]
	locs = ["a"=>(xmid, y[1]), "b"=>(xmin, ymid), "c"=>(xmid, ymax), "d"=>(xmax, ymid),
		"e"=>(xmid, y[2]), "f"=>(x[2], ymid), "g"=>(xmid, y[end-1]),
		"h"=>(x[end-1], ymid), "i"=>(x[3], y[3]), "j"=>(x[4], y[3]),
		"k"=>(x[5], y[3]), "l"=>(x[6], y[3]), "m"=>(x[3], y[4]),
		"n"=>(x[4], y[4]), "o"=>(x[5], y[4]), "p"=>(x[6], y[4])]
    UnitdiskWidget(locs, 0.23)
end

function widget12()
	a = 0.12
	ymid = xmid = 0.5
	X = 0.33
	Y = 0.17
	D = 0.15
	y = [ymid-Y, ymid-Y+D, ymid-a/2, ymid+a/2, ymid+Y-D, ymid+Y]
	x = [xmid-X, xmid-X+D, xmid-1.5a, xmid-a/2, xmid+a/2, xmid+1.5a, xmid+X-D, xmid+X]
	xmin, xmax, ymin, ymax = x[1], x[end], y[1], y[end]
	locs = ["a"=>(xmid, y[1]), "b"=>(xmin, ymid), "c"=>(xmid, ymax), "d"=>(xmax, ymid),
		"i"=>(x[3], y[3]), "j"=>(x[4], y[3]),
		"k"=>(x[5], y[3]), "l"=>(x[6], y[3]), "m"=>(x[3], y[4]),
		"n"=>(x[4], y[4]), "o"=>(x[5], y[4]), "p"=>(x[6], y[4])]
    UnitdiskWidget(locs, 0.23)
end

function contract16(x::AbstractVector{T}) where T
	code = ein"a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,ae,bf,cg,dh,ei,ej,ek,el,fi,fm,gm,gn,go,gp,hl,hp,ij,im,in,jk,jm,jn,jo,kl,kn,ko,kp,lo,lp,mn,no,op->abcd"
    code = optimize_greedy(code, uniformsize(code, 2))
    code([misv(xi) for xi in x]..., [misb(T) for i=1:32]...)
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
