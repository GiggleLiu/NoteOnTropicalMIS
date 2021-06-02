using Polynomials
export independence_polynomial

using FFTW
function independence_polynomial(::Val{:fft}, code; mis_size=Int(mis_solve(code)[].n), r=1.0)
	ω = exp(-2im*π/(mis_size+1))
	xs = r .* collect(ω .^ (0:mis_size))
	ys = [mis_contract(x, code)[] for x in xs]
	Polynomial(real.(ifft(ys) ./ (r .^ (0:mis_size))))
end

function independence_polynomial(::Val{:fitting}, code; mis_size=Int(mis_solve(code)[].n))
	xs = (0:mis_size)
	ys = [mis_contract(x, code)[] for x in xs]
	fit(xs, ys, mis_size)
end

function independence_polynomial(::Val{:polynomial}, code)
    mis_contract(Polynomial([0, 1.0]), code)[]
end

using Mods, Primes

# pirate
Base.abs(x::Mod) = x
Base.isless(x::Mod{N}, y::Mod{N}) where N = mod(x.val, N) < mod(y.val, N)

function _independance_polynomial(::Type{T}, code, mis_size::Int) where T
	xs = 0:mis_size
	ys = [mis_contract(T(x), code)[] for x in xs]
	A = zeros(T, mis_size+1, mis_size+1)
	for j=1:mis_size+1, i=1:mis_size+1
		A[j,i] = T(xs[j])^(i-1)
	end
	A \ T.(ys)
end

function independence_polynomial(::Val{:finitefield}, code; mis_size=Int(mis_solve(code)[].n), max_order=100)
    N = typemax(Int)
    YS = []
    local res
    for k = 1:max_order
	    N = prevprime(N-1)
        T = Mods.Mod{N,Int}
        rk = _independance_polynomial(T, code, mis_size)
        push!(YS, rk)
        if max_order==1
            return Polynomial(Mods.value.(YS[1]))
        elseif k != 1
            ra = improved_counting(YS[1:end-1])
            res = improved_counting(YS)
            ra == res && return Polynomial(res)
        end
    end
    @warn "result is potentially inconsistent."
    return Polynomial(res)
end

function improved_counting(sequences)
    map(yi->Mods.CRT(yi...), zip(sequences...))
end