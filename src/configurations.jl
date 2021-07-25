export mis_config, ConfigEnumerator, ConfigTropical

struct ConfigEnumerator{N,C}
    data::Vector{StaticBitVector{N,C}}
end

Base.length(x::ConfigEnumerator{N}) where N = length(x.data)
Base.:(==)(x::ConfigEnumerator{N,C}, y::ConfigEnumerator{N,C}) where {N,C} = x.data == y.data

function Base.:+(x::ConfigEnumerator{N,C}, y::ConfigEnumerator{N,C}) where {N,C}
    length(x) == 0 && return y
    length(y) == 0 && return x
    return ConfigEnumerator{N,C}(vcat(x.data, y.data))
end

function Base.:*(x::ConfigEnumerator{L,C}, y::ConfigEnumerator{L,C}) where {L,C}
    M, N = length(x), length(y)
    M == 0 && return x
    N == 0 && return y
    z = Vector{StaticBitVector{L,C}}(undef, M*N)
    @inbounds for j=1:N, i=1:M
        z[(j-1)*M+i] = x.data[i] | y.data[j]
    end
    return ConfigEnumerator{L,C}(z)
end

Base.zero(::Type{ConfigEnumerator{N,C}}) where {N,C} = ConfigEnumerator{N,C}(StaticBitVector{N,C}[])
Base.one(::Type{ConfigEnumerator{N,C}}) where {N,C} = ConfigEnumerator{N,C}([TropicalNumbers.staticfalses(StaticBitVector{N,C})])

function symbols(::EinCode{ixs}) where ixs
    res = []
    for ix in ixs
        for l in ix
            if l ∉ res
                push!(res, l)
            end
        end
    end
    return res
end

function mis_config(code; all=false, bounding=true, usecuda=false)
    if all && usecuda
        throw(ArgumentError("ConfigEnumerator can not be computed on GPU!"))
    end
    flatten_code = flatten(code)
    vertex_index = Dict([s=>i for (i, s) in enumerate(symbols(flatten_code))])
    N = length(vertex_index)
    C = TropicalNumbers._nints(N)
    xs = map(getixs(flatten_code)) do ix
        T = all ? CountingTropical{Float64, ConfigEnumerator{N,C}} : ConfigTropical{Float64, N, C}
        if length(ix) == 2
            return misb(T)
        else
            s = TropicalNumbers.onehot(StaticBitVector{N,C}, vertex_index[ix[1]])
            if all
                misv(T, T(1.0, ConfigEnumerator([s])))
            else
                misv(T, T(1.0, s))
            end
        end
    end
    if bounding
        ymask = trues(fill(2, length(getiy(flatten_code)))...)
        xst = map(getixs(flatten_code)) do ix
            length(ix) == 1 ? misv(TropicalF64,Tropical(1.0)) : misb(TropicalF64)
        end
        if usecuda
            ymask = CuArray(ymask)
            xst = CuArray.(xst)
        end
        if all
            return bounding_contract(code, xst, ymask, xs)
        else
            @assert ndims(ymask) == 0
            t, res = mis_config_ad(code, xst, ymask)
            return fill(ConfigTropical(asscalar(t).n, StaticBitVector(map(l->res[l], 1:N))))
        end
    else
        if usecuda
            xs = CuArray.(xs)
        end
	    return code(xs...)
    end
end