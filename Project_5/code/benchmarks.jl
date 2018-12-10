using PyCall
using BenchmarkTools
include("SIRS_Basic.jl")
include("SIRS_MCMC.jl")


@btime SIRS_basic_pylike(S0=300, I0=100, R0=0, a=4, b=1, c=0.5, stop_time=10)
@btime SIRS_basic(S0=300, I0=100, R0=0, a=4, b=1, c=0.5, stop_time=10)
