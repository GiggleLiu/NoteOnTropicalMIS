julia> using GraphTensorNetworks, Random, Graphs

julia> graph = (Random.seed!(2); Graphs.smallgraph(:petersen))
{10, 15} undirected simple Int64 graph

julia> problem = Independence(graph; optimizer=TreeSA(sc_target=0, sc_weight=1.0, ntrials=10, βs=0.01:0.1:15.0, niters=20, rw_weight=0.2));
┌ Warning: target space complexity not found, got: 4.0, with time complexity 7.965784284662087, read-right complexity 8.661778097771988.
└ @ OMEinsumContractionOrders ~/.julia/dev/OMEinsumContractionOrders/src/treesa.jl:71
time/space complexity is (7.965784284662086, 4.0)

# maximum independent set size
julia> solve(problem, "size max")
0-dimensional Array{TropicalNumbers.TropicalF64, 0}:
4.0ₜ

# all independent sets
julia> solve(problem, "counting sum")
0-dimensional Array{Float64, 0}:
76.0

# counting maximum independent sets
julia> solve(problem, "counting max")
0-dimensional Array{TropicalNumbers.CountingTropicalF64, 0}:
(4.0, 5.0)ₜ

# counting independent sets of max two sizes
julia> solve(problem, "counting max2")
0-dimensional Array{Max2Poly{Float64, Float64}, 0}:
30.0*x^3 + 5.0*x^4

# using `Polynomial` type
julia> solve(problem, "counting all")
0-dimensional Array{Polynomial{Float64, :x}, 0}:
Polynomial(1.0 + 10.0*x + 30.0*x^2 + 30.0*x^3 + 5.0*x^4)

# using the finitefield approach
julia> solve(problem, "counting all (finitefield)")
0-dimensional Array{Polynomial{BigInt, :x}, 0}:
Polynomial(1 + 10*x + 30*x^2 + 30*x^3 + 5*x^4)

# using the fourier approach
julia> solve(problem, "counting all (fft)", r=1.0)
0-dimensional Array{Polynomial{ComplexF64, :x}, 0}:
Polynomial(1.0000000000000029 + 2.664535259100376e-16im + (10.000000000000004 - 1.9512435398857492e-16im)x + (30.0 - 1.9622216671393801e-16im)x^2 + (30.0 + 1.1553104311877194e-15im)x^3 + (5.0 - 1.030417436395244e-15im)x^4)

# one of MISs
julia> solve(problem, "config max")
0-dimensional Array{CountingTropical{Float64, ConfigSampler{10, 1, 1}}, 0}:
(4.0, ConfigSampler{10, 1, 1}(1010000011))ₜ

julia> solve(problem, "config max (bounded)")
0-dimensional Array{CountingTropical{Float64, ConfigSampler{10, 1, 1}}, 0}:
(4.0, ConfigSampler{10, 1, 1}(1010000011))ₜ

# enumerate all MISs
julia> solve(problem, "configs max")  # not recommended
0-dimensional Array{CountingTropical{Float64, ConfigEnumerator{10, 1, 1}}, 0}:
(4.0, {1010000011, 0100100110, 1001001100, 0010111000, 0101010001})ₜ

julia> solve(problem, "configs max (bounded)")
0-dimensional Array{CountingTropical{Int64, ConfigEnumerator{10, 1, 1}}, 0}:
(4, {1010000011, 0100100110, 1001001100, 0010111000, 0101010001})ₜ

# enumerate all configurations of independent sets of size |MIS| and |MIS|-1
julia> solve(problem, "configs max2")
0-dimensional Array{Max2Poly{ConfigEnumerator{10, 1, 1}, Float64}, 0}:
{0010101000, 0101000001, 0100100010, 0010100010, 0100000011, 0010000011, 1001001000, 1010001000, 1001000001, 1010000001, 1010000010, 1000000011, 0100100100, 0000101100, 0101000100, 0001001100, 0000100110, 0100000110, 1001000100, 1000001100, 1000000110, 0100110000, 0000111000, 0101010000, 0001011000, 0010110000, 0010011000, 0001010001, 0100010001, 0010010001}*x^3 + {1010000011, 0100100110, 1001001100, 0010111000, 0101010001}*x^4

# enumerate all IS configurations
julia> solve(problem, "configs all")
0-dimensional Array{Polynomial{ConfigEnumerator{10, 1, 1}, :x}, 0}:
Polynomial({0000000000} + {0010000000, 0000100000, 0001000000, 0100000000, 0000001000, 0000000001, 0000000010, 1000000000, 0000000100, 0000010000}*x + {1000000010, 0010100000, 0010001000, 0100100000, 0000101000, 0101000000, 0001001000, 0001000001, 0100000001, 0010000001, 0000100010, 0100000010, 0010000010, 0000000011, 1001000000, 1000001000, 1010000000, 1000000001, 0000000110, 0000100100, 0001000100, 0100000100, 0000001100, 1000000100, 0010010000, 0000110000, 0001010000, 0100010000, 0000011000, 0000010001}*x^2 + {1010000010, 1000000011, 0010101000, 0101000001, 0100100010, 0010100010, 0100000011, 0010000011, 1001001000, 1010001000, 1001000001, 1010000001, 0000100110, 0100000110, 0100100100, 0000101100, 0101000100, 0001001100, 1001000100, 1000001100, 1000000110, 0010110000, 0010011000, 0100110000, 0000111000, 0101010000, 0001011000, 0001010001, 0100010001, 0010010001}*x^3 + {1010000011, 0100100110, 1001001100, 0010111000, 0101010001}*x^4)
