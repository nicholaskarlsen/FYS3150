using PyCall
using BenchmarkTools
include("SIRS_MC_Basic.jl")
include("SIRS_MC_vax.jl")

println("SIRS_basic:")
@btime SIRS_basic(S0=300, I0=100, R0=0, a=4, b=1, c=0.5, stop_time=10, trials=100)
println("SIRS_vax:")
@btime SIRS_vax(S0=300, I0=100, R0=0, a=4, b=1, c=0.5, f=200, stop_time=10, trials=100)


