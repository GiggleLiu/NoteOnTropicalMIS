struct ConfigTropical{T,N}
    n::T
    config::BitVector
end

Base.:(==)(x::ConfigTropical{T,N}, y::ConfigTropical{T,N}) where {T,N} = x.n == y.n && x.config == y.config

function Base.:+(x::ConfigTropical{T,N}, y::ConfigTropical{T,N}) where {T,N}
    return x.n > y.n ? x : y
end

function Base.:*(x::ConfigTropical{T,L}, y::ConfigTropical{T,L}) where {T,L}
    return ConfigTropical{T,L}(x.n+y.n, x.config .| y.config)
end

Base.zero(::Type{ConfigTropical{T, N}}) where {T,N} = ConfigTropical{T,N}(zero(Tropical{T}).n, trues(N))
Base.one(::Type{ConfigTropical{T, N}}) where {T,N} = ConfigTropical{T,N}(one(Tropical{T}).n, falses(N))
function onehot(::Type{ConfigTropical{T,N}}, i::Int) where {T,N}
    c = falses(N)
    c[i] |= true
    return ConfigTropical{T,N}(one(T), c)
end