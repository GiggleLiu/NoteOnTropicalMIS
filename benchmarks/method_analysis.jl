using MethodAnalysis
using OMEinsum, GenericTensorNetworks, LightGraphs

for i=1:10
    print(i)
    g = LightGraphs.random_regular_graph(100, 3)
    id = Independence(g; optmethod=:tree, niters=10, ntrials=1)
    solve(id, "counting all (finitefield)")
end

meths = []
visit(OMEinsum) do item
    isa(item, Method) && push!(meths, item)
    true   # walk through everything
end
argmax(x->length(methodinstances(x)), meths)
sorted_meths = sort(meths, by=x->length(methodinstances(x)), rev=true)
println(length.(methodinstances.(sorted_meths[1:10])))