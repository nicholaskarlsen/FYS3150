# Check how much PyJulia affects execution time

import numpy as np
import time
from julia import Main as jcall

t0 = time.time()
jcall.include("SIRS_MC_Basic.jl")
print "including both files took:",  time.time() - t0
timings = []
for i in range(100):
    t0 = time.time()
    jcall.eval("SIRS_basic(S0=300, I0=100, R0=0, a=4, b=1, c=0.5, stop_time=10, trials=100)")
    timings.append(time.time() - t0)

print "mean time = %.3e" % np.mean(timings)