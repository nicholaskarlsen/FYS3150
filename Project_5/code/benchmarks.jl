using PyCall
using BenchmarkTools
include("SIRS_Basic.jl")
include("SIRS_MC.jl")

println("SIRS_basic_pylike:")
@btime SIRS_basic_pylike(S0=300, I0=100, R0=0, a=4, b=1, c=0.5, stop_time=10)
println("SIRS_basic:")
@btime SIRS_basic(S0=300, I0=100, R0=0, a=4, b=1, c=0.5, stop_time=10)

#= Output in terminal (LAPTOP):

~/D/u/F/P/code   master ±  julia benchmarks.jl 
SIRS_basic_pylike:
  90.920 μs (11 allocations: 94.41 KiB)
SIRS_basic:
  519.826 μs (16061 allocations: 947.94 KiB)
=#