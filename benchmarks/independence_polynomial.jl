using NoteOnTropicalMIS, DelimitedFiles
using LightGraphs
using BenchmarkTools
using CUDA

CASE_SET = :r3
cases = if CASE_SET == :r3
    [case_r3(n, 3; seed=2, sc_target=s) for (n, s) in [
        (10, 3), (20, 4), (30, 6), (40, 7), (50, 8), (60, 10), (70, 11), (80, 11), (90, 15), (100, 15),
        (110, 15), (120, 18), (130, 17), (140, 16), (150, 21), (160, 20), (170, 22), (180, 24), (190, 26), (200, 25),
    ]]
else
    [case_dc(L, 0.8; seed=2, sc_target=s) for (L, s) in [
        (4, 5), (6, 7), (8, 9), (10, 8), (12, 12), (14, 13), (16, 17), (18, 18), (20, 18), (22, 23), (24, 23),
    ]]
end

function run_benchmarks(cases; output_file)
    suite = BenchmarkGroup()
    for (case, f) in cases
        suite[case] = @benchmarkable $f()
    end

    tune!(suite)
    res = run(suite)

    times = zeros(length(cases))
    for (k, (case, f)) in enumerate(cases)
        times[k] = minimum(res[case].times)
    end

    println("Writing benchmark results to file: $output_file.")
    mkpath(dirname(output_file))
    writedlm(output_file, times)
end

run_benchmarks([("n$(10*i)", ()->run_task(case, :maxsize)) for (i, case) in enumerate(cases[1:3])],
            output_file=NoteOnTropicalMIS.project_relative_path("benchmarks", "data", "maxsize-r3.dat"))