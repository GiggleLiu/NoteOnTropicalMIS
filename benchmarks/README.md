# Run Benchmarks

Open a bash terminal and go to the folder this RAEDME is in. Then type

```bash
$ julia --project -e "using Pkg; Pkg.instantiate()"

$ julia --project benchmarks.jl cpu 1
```

Here, the 3rd argument is the benchmark group name, please check the code for details. For GPU benchmarks, just type

```bash
$ julia --project benchmarks.jl gpu 1
```

For maximal independent sets, change filename to `benchmarks_maximal.jl`.
