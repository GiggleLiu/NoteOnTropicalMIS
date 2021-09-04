using NoteOnTropicalMIS, LightGraphs
using StatsBase, TropicalNumbers

function random_graph(n::Int, m::Int)
    g = SimpleGraph(n)
    for i=1:m
        (a,b) = sample(1:n, 2; replace=true)
        add_edge!(g, a, b)
    end
    return g
end

function Base.count_ones(x::StaticBitVector)
    sum(v->count_ones(v),x.data)
end

function phase_space_reachable(g)
    code = NoteOnTropicalMIS.idp_code(g; method=:greedy)
    s, c = (res = mis_contract(CountingTropical{Int,Int}(1,1), code)[]; (res.n, res.c))
    println("Graph size = $(nv(g)), MIS size = $s, degeneracy = $c")
    s2, config1, config0 = (res = mis_max2_config(code)[]; (res.maxorder, res.a, res.b))
    @assert s == s2
    @assert length(config0) == c
    configs = config0.data ∪ config1.data
    @show config0, config1
    G = SimpleGraph(length(configs))
    for (ic,c) in enumerate(configs)
        for (id,d) in enumerate(configs)
            diffpos = c ⊻ d
            ndiff = count_ones(diffpos)
            if ndiff == 1
                add_edge!(G, ic, id)
            elseif ndiff == 2
                a = 0
                b = 0
                for (ip,p) in enumerate(diffpos)
                    if p == 1
                        if a == 0
                            a = ip
                        else
                            b = ip
                        end
                    end
                end
                if has_edge(g, a, b)
                    add_edge!(G, ic, id)
                end
            end
        end
    end
    @show connected_components(G)
    all(sg->any(v->v<=length(config0.data), sg), connected_components(G))
end

#=
for s in [
        (:diamond        , 4, 5, false),
        (:bull           , 5, 5, false),
        (:chvatal        , 12, 24, false),
        (:cubical        , 8, 12, false),
        (:desargues      , 20, 30, false),
        (:dodecahedral   , 20, 30, false),
        (:frucht         , 12, 18, false),
        (:heawood        , 14, 21, false),
        (:house          , 5,6, false),
        (:housex         , 5, 8, false),
        (:icosahedral    , 12, 30, false),
        (:krackhardtkite , 10, 18, false),
        (:moebiuskantor  , 16, 24, false),
        (:octahedral     , 6, 12, false),
        (:pappus         , 18, 27, false),
        (:petersen       , 10, 15, false),
        (:sedgewickmaze  , 8, 10, false),
        (:tetrahedral    , 4, 6, false),
        (:truncatedcube  , 24, 36, false),
        (:truncatedtetrahedron , 12, 18, false),
        (:tutte          , 46, 69, false)]
    @show s
    g = LightGraphs.smallgraph(s[1])
    @assert phase_space_reachable(g)
end
=#

using Random
for i=1:100
    Random.seed!(i)
    println(i)
    g = random_graph(5, 8)
    if !phase_space_reachable(g)
        break
    end
end