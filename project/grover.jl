using Yao, LightGraphs
using Yao.ConstGate: P1

# U*nᵢnⱼ + Ω*X + Δ*Z
function rydberg_hamiltonian(g::SimpleGraph, U, Ω, Δ)
    n = nv(g)
    res = AbstractBlock{n}[]
    U != 0 && push!(res, U * +([put(n, (e.src, e.dst)=>kron(P1, P1)) for e in edges(g)]...))
    Ω != 0 && push!(res, Ω * +([put(n, v=>X) for v in vertices(g)]...))
    Δ != 0 && push!(res, Δ * +([put(n, v=>(-P1)) for v in vertices(g)]...))
    return Add(res)
end

rydberg_hamiltonian(random_regular_graph(10, 3), 100.0, 0.0, 1.0)