export mis_config, ConfigEnumerator, ConfigTropical

struct ConfigEnumerator{N}
    data::Vector{BitVector}
end

Base.length(x::ConfigEnumerator{N}) where N = length(x.data)
Base.:(==)(x::ConfigEnumerator, y::ConfigEnumerator) = x.data == y.data

function Base.:+(x::ConfigEnumerator{N}, y::ConfigEnumerator{N}) where N
    res = ConfigEnumerator{N}(vcat(x.data, y.data))
    return res
end

function Base.:*(x::ConfigEnumerator{L}, y::ConfigEnumerator{L}) where L
    M, N = length(x), length(y)
    z = Vector{BitVector}(undef, M*N)
    for j=1:N, i=1:M
        z[(j-1)*M+i] = x.data[i] .| y.data[j]
    end
    return ConfigEnumerator{L}(z)
end

Base.zero(::Type{ConfigEnumerator{N}}) where N = ConfigEnumerator{N}(BitVector[])
Base.one(::Type{ConfigEnumerator{N}}) where N = ConfigEnumerator{N}([falses(N)])

function onehot(::Type{ConfigEnumerator{N}}, i::Int) where N
    res = one(ConfigEnumerator{N})
    res.data[1][i] |= true
    return res
end

include("configtropical.jl")

function enumerator_t(::Type{T}, ix::NTuple{2}, vertex_index) where {T1, T2, T<:CountingTropical{T1, T2}}
    [one(T) one(T); one(T) zero(T)]
end
function enumerator_t(::Type{T}, ix::NTuple{1}, vertex_index) where {T1, T2, T<:CountingTropical{T1, T2}}
    [one(T), CountingTropical(one(T1), onehot(T2, vertex_index[ix[1]]))]
end
function enumerator_t(::Type{T}, ix::NTuple{2}, vertex_index) where {T1, T2, T<:ConfigTropical{T1, T2}}
    [one(T) one(T); one(T) zero(T)]
end
function enumerator_t(::Type{T}, ix::NTuple{1}, vertex_index) where {T1, T2, T<:ConfigTropical{T1, T2}}
    [one(T), onehot(T, vertex_index[ix[1]])]
end

symbols(::EinCode{ixs}) where ixs = unique(Iterators.flatten(ixs))
symbols(ne::NestedEinsum) = symbols(Iterators.flatten(ne))
function mis_config(code; all=false)
    vertex_index = Dict([s=>i for (i, s) in enumerate(symbols(code))])
    T = all ? CountingTropical{Float64, ConfigEnumerator{length(vertex_index)}} : ConfigTropical{Float64, length(vertex_index)}
	xs = generate_xs!((T, ix)->enumerator_t(T, ix, vertex_index), T, code, Vector{Any}(undef, ninput(code)))
	code(xs...)
end