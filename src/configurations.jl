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

symbols(::EinCode{ixs}) where ixs = unique(Iterators.flatten(ixs))
symbols(ne::OMEinsum.NestedEinsum) = symbols(Iterators.flatten(ne))
single_symbols(::EinCode{ixs}) where ixs = unique(Iterators.flatten(filter(x->length(x)==1,ixs)))
single_symbols(ne::OMEinsum.NestedEinsum) = single_symbols(Iterators.flatten(ne))
ninput(::EinCode{ixs}) where ixs = length(ixs)
ninput(ne::OMEinsum.NestedEinsum) = ninput(Iterators.flatten(ne))
_getiy(code::EinCode) = OMEinsum.getiy(code)
_getiy(code::NestedEinsum) = OMEinsum.getiy(code.eins)
function mis_config(code; all=false, usemask=true)
    vertex_index = Dict([s=>i for (i, s) in enumerate(single_symbols(code))])
    N = length(vertex_index)
    C = TropicalNumbers._nints(N)
    xs = generate_vertextensors(code) do ix
        T = all ? CountingTropical{Float64, ConfigEnumerator{N,C}} : ConfigTropical{Float64, N, C}
        if length(ix) == 2
            return [one(T) one(T); one(T) zero(T)]
        else
            s = TropicalNumbers.onehot(StaticBitVector{N,C}, vertex_index[ix[1]])
            if all
                [one(T), T(1.0, ConfigEnumerator([s]))]
            else
                [one(T), T(1.0, s)]
            end
        end
    end
    if usemask
        ymask = trues(fill(2, length(_getiy(code)))...)
        xst = generate_vertextensors(ix->length(ix)==1 ? misv(TropicalF64,1,Tropical(1.0)) : misb(TropicalF64), code)
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