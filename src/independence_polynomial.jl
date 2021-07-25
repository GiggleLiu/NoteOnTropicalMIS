using Polynomials
using OMEinsum: NestedEinsum, getixs, getiy

export independence_polynomial
export misb, misv, mis_size, mis_count, mis_contract

# MIS bond tensor
misb(::Type{T}) where T = [one(T) one(T); one(T) zero(T)]
# MIS vertex tensor
misv(::Type{T}, val) where T = [one(T), convert(T, val)]

function neighbortensor(x::T, d::Int) where T
    t = zeros(T, fill(2, d)...)
    for i = 2:1<<(d-1)
        t[i] = one(T)
    end
    t[1<<(d-1)+1] = x
    return t
end

function mis_contract(x::T, code; usecuda=false) where {T}
	tensors = map(getixs(flatten(code))) do ix
        # if the tensor rank is 1, create a vertex tensor.
        # otherwise the tensor rank must be 2, create a bond tensor.
        t = length(ix)==1 ? misv(T, x) : misb(T)
        usecuda ? CuArray(t) : t
    end
	code(tensors...)
end

mis_size(code; usecuda=false) = Int(asscalar(mis_contract(TropicalF64(1.0), code; usecuda=usecuda)).n)
mis_count(code; usecuda=false) = asscalar(mis_contract(CountingTropical{Float64,Float64}(1.0, 1.0), code; usecuda=usecuda)).c

using FFTW
function independence_polynomial(::Val{:fft}, code; mis_size=mis_size(code), r=1.0, usecuda=false)
	ω = exp(-2im*π/(mis_size+1))
	xs = r .* collect(ω .^ (0:mis_size))
	ys = [asscalar(mis_contract(x, code; usecuda=usecuda)) for x in xs]
	Polynomial(ifft(ys) ./ (r .^ (0:mis_size)))
end

function independence_polynomial(::Val{:fitting}, code; mis_size=mis_size(code), usecuda=false)
	xs = (0:mis_size)
	ys = [asscalar(mis_contract(x, code; usecuda=usecuda)) for x in xs]
	fit(xs, ys, mis_size)
end

function independence_polynomial(::Val{:polynomial}, code; usecuda=false)
    @assert !usecuda "Polynomial type can not be computed on GPU!"
    mis_contract(Polynomial([0, 1.0]), code)
end

using Mods, Primes

# pirate
Base.abs(x::Mod) = x
Base.isless(x::Mod{N}, y::Mod{N}) where N = mod(x.val, N) < mod(y.val, N)

function _independance_polynomial(::Type{T}, code, mis_size::Int; usecuda) where T
	xs = 0:mis_size
	ys = [asscalar(mis_contract(T(x), code; usecuda=usecuda)) for x in xs]
	A = zeros(T, mis_size+1, mis_size+1)
	for j=1:mis_size+1, i=1:mis_size+1
		A[j,i] = T(xs[j])^(i-1)
	end
	A \ T.(ys)
end

function independence_polynomial(::Val{:finitefield}, code; mis_size=mis_size(code), max_order=100, usecuda=false)
    TI = Int32  # Int 32 is faster
    N = typemax(TI)
    YS = []
    local res
    for k = 1:max_order
	    N = prevprime(N-TI(1))
        T = Mods.Mod{N,TI}
        rk = _independance_polynomial(T, code, mis_size; usecuda=usecuda)
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

export maximal_contract, maximal_code

function idp_code(g::SimpleGraph; method=:kahypar, sc_target=17, max_group_size=40, nrepeat=10, imbalances=0.0:0.01:0.2)
    code = EinCode(([minmax(e.src,e.dst) for e in LightGraphs.edges(g)]..., # labels for edge tensors
                    [(i,) for i in LightGraphs.vertices(g)]...), ())        # labels for vertex tensors
    size_dict = Dict([s=>2 for s in symbols(code)])
    optcode = if method == :kahypar
        optimize_kahypar(code, size_dict; sc_target=sc_target, max_group_size=max_group_size, imbalances=imbalances)
    elseif method == :greedy
        optimize_greedy(code, size_dict; nrepeat=nrepeat)
    else
        ArgumentError("optimizer `$method` not defined.")
    end
    println("time/space complexity is $(OMEinsum.timespace_complexity(optcode, size_dict))")
    return optcode
end

function maximal_code(g::SimpleGraph; method=:kahypar, sc_target=17, max_group_size=40, nrepeat=10, imbalances=0.0:0.01:0.5)
    code = EinCode(([(LightGraphs.neighbors(g, v)..., v) for v in LightGraphs.vertices(g)]...,), ())
    size_dict = Dict([s=>2 for s in symbols(code)])
    optcode = if method == :kahypar
        optimize_kahypar(code, size_dict; sc_target=sc_target, max_group_size=max_group_size, imbalances=imbalances)
    elseif method == :greedy
        optimize_greedy(code, size_dict; nrepeat=nrepeat)
    else
        ArgumentError("optimizer `$method` not defined.")
    end
    println("time/space complexity is $(OMEinsum.timespace_complexity(optcode, size_dict))")
    return optcode
end

function maximal_contract(optcode::NestedEinsum, x::T; usecuda=false) where T
    ixs = OMEinsum.getixs(OMEinsum.flatten(optcode))
	tensors = map(ixs) do ix
        t = neighbortensor(x, length(ix))
        usecuda ? CuArray(t) : t
    end
	optcode(tensors...)
end

function maximal_polynomial(::Val{:fft}, g; usecuda=false, mis_size=run_task(g, :maxsize; usecuda=usecuda)[].n, r=1.0, kwargs...)
	ω = exp(-2im*π/(mis_size+1))
	xs = r .* collect(ω .^ (0:mis_size))
    optcode = maximal_code(g; kwargs...)
	ys = [OMEinsum.asscalar(maximal_contract(optcode, x; usecuda=usecuda)) for x in xs]
	Polynomial(ifft(ys) ./ (r .^ (0:mis_size)))
end

function maximal_polynomial(::Val{:polynomial}, g; usecuda=false, kwargs...)
    @assert !usecuda "Polynomial type can not be computed on GPU!"
    optcode = maximal_code(g; kwargs...)
    maximal_contract(optcode, Polynomial([0, 1.0]))
end