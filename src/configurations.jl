export mis_config, ConfigEnumerator, ConfigTropical

struct ConfigEnumerator{N,C}
    data::Vector{StaticBitVector{N,C}}
end

Base.length(x::ConfigEnumerator{N}) where N = length(x.data)
Base.:(==)(x::ConfigEnumerator{N,C}, y::ConfigEnumerator{N,C}) where {N,C} = x.data == y.data

function Base.:+(x::ConfigEnumerator{N,C}, y::ConfigEnumerator{N,C}) where {N,C}
    res = ConfigEnumerator{N,C}(vcat(x.data, y.data))
    return res
end

function Base.:*(x::ConfigEnumerator{L,C}, y::ConfigEnumerator{L,C}) where {L,C}
    M, N = length(x), length(y)
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
ninput(::EinCode{ixs}) where ixs = length(ixs)
ninput(ne::OMEinsum.NestedEinsum) = ninput(Iterators.flatten(ne))
function mis_config(code; all=false)
    vertex_index = Dict([s=>i for (i, s) in enumerate(symbols(code))])
    N = length(vertex_index)
    C = TropicalNumbers._nints(N)
    T = all ? CountingTropical{Float64, ConfigEnumerator{N,C}} : ConfigTropical{Float64, N, C}
	xs = generate_xs!((T, ix)->enumerator_t(T, ix, vertex_index), T, code, Vector{Any}(undef, ninput(code)))
	code(xs...)
end