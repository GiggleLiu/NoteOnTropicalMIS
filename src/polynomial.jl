export mis_polysolve

_auto_mispolytensor(x::T, ix::NTuple{2}) where T = T[1 1; 1 0]
_auto_mispolytensor(x::T, ix::NTuple{1}) where T = T[1, x]
function generate_polyxs!(x::T, code::OMEinsum.NestedEinsum, xs=Vector{Any}(undef, ninput(code))) where {T}
    for (ix, arg) in zip(OMEinsum.getixs(code.eins), code.args)
		if arg isa Integer
			xs[arg] = _auto_mispolytensor(x, ix)
		else
        	generate_polyxs!(x, arg, xs)
		end
    end
	return xs
end

function generate_polyxs!(x::T, code::EinCode, xs=Vector{Any}(undef, ninput(code))) where {T}
    for (i,ix) in enumerate(OMEinsum.getixs(code))
		xs[i] = _auto_mispolytensor(x, ix)
    end
	return xs
end

function mis_polysolve(code, x::T) where {T}
	xs = generate_polyxs!(x, code, Vector{Any}(undef, NoteOnTropicalMIS.ninput(code)))
	code(xs...)
end