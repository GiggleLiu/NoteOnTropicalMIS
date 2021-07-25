using NoteOnTropicalMIS, Polynomials, Random, OMEinsum, FFTW

Random.seed!(2)
g = diagonal_coupled_graph(rand(13,13) .<= 0.8)
@time maximal_polynomial(Val(:fft), g; sc_target=24, imbalances=0.0:0.00151:1.0, max_group_size=50)
@time maximal_polynomial(Val(:polynomial), g; sc_target=19, imbalances=0.0:0.00151:1.0, max_group_size=50)