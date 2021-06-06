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

function enumerator_t(::Type{T}, ix::NTuple{2}, vertex_index) where {T1, N, C, T<:CountingTropical{T1,ConfigEnumerator{N,C}}}
    [one(T) one(T); one(T) zero(T)]
end
function enumerator_t(::Type{T}, ix::NTuple{1}, vertex_index) where {T1, N, C, T<:CountingTropical{T1,ConfigEnumerator{N,C}}}
    [one(T), CountingTropical(one(T1), ConfigEnumerator([TropicalNumbers.onehot(StaticBitVector{N,C}, vertex_index[ix[1]])]))]
end
function enumerator_t(::Type{T}, ix::NTuple{2}, vertex_index) where {T1,N,C, T<:ConfigTropical{T1,N,C}}
    [one(T) one(T); one(T) zero(T)]
end
function enumerator_t(::Type{T}, ix::NTuple{1}, vertex_index) where {T1,N,C, T<:ConfigTropical{T1,N,C}}
    [one(T), T(one(T1), TropicalNumbers.onehot(StaticBitVector{N,C}, vertex_index[ix[1]]))]
end

symbols(::EinCode{ixs}) where ixs = unique(Iterators.flatten(ixs))
symbols(ne::OMEinsum.NestedEinsum) = symbols(Iterators.flatten(ne))
single_symbols(::EinCode{ixs}) where ixs = unique(Iterators.flatten(filter(x->length(x)==1,ixs)))
single_symbols(ne::OMEinsum.NestedEinsum) = single_symbols(Iterators.flatten(ne))
ninput(::EinCode{ixs}) where ixs = length(ixs)
ninput(ne::OMEinsum.NestedEinsum) = ninput(Iterators.flatten(ne))
_getiy(code::EinCode) = OMEinsum.getiy(code)
_getiy(code::NestedEinsum) = OMEinsum.getiy(code.eins)
function mis_config(code; all=false, usemask=true)
    nvertex = ninput(code)
    vertex_index = Dict([s=>i for (i, s) in enumerate(single_symbols(code))])
    N = length(vertex_index)
    C = TropicalNumbers._nints(N)
    T = all ? CountingTropical{Float64, ConfigEnumerator{N,C}} : ConfigTropical{Float64, N, C}
    xs = generate_xs!((T, ix)->enumerator_t(T, ix, vertex_index), T, code, Vector{Any}(undef, nvertex))
    if usemask
        ymask = trues(fill(2, length(_getiy(code)))...)
        xst = generate_xs!(_auto_mistensor, Tropical(1.0), code, Vector{Any}(undef, nvertex))
        if all
            return bounding_contract(code, xst, ymask, xs)
        else
            @assert ndims(ymask) == 0
            t, res = mis_config_ad(code, xst, ymask)
            return fill(ConfigTropical(t[].n, StaticBitVector(map(l->res[l], 1:N))))
        end
    else
	    return code(xs...)
    end
end
function Base.Matrix(x::ConfigEnumerator{N,C}) where {N,C}
    m = zeros(UInt64, C, length(x))
    for i=1:length(x)
        m[:,i] .= x.data[i].data
    end
    return m
end