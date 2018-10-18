#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Written by candidate 15205
"""
Solves exercise 2
"""

import numpy as np
import matplotlib.pyplot as plt
import solver as fys

# Initialize solver with init conditions
x0 = 1
v0 = 0
t0 = 0
m = 500e-3
tn = 20
dt = 1e-2
k = 1
b = 0.1

sys1 = fys.diffeqsolver(x_init=x0,
                        v_init=v0,
                        t_init=t0)

# Create function for the Diff eq


def diffEq1(x, v, t):
    "a(t) = (-kx(t) - bv(t))/m"
  #  return(- k * x - b * v) / m
    return(- k * x) / m


# Solve Diffeq
pos, vel, t = sys1.solver(tn, dt, sys1.rk4, diffEq1)

# Plot phase space
plt.figure(figsize=(8, 4), dpi=100)
plt.subplot(121)
plt.plot(pos, vel)
plt.xlabel("$x(t)\enspace[m]$")
plt.axhline(0, color="black", linewidth=0.75, linestyle="-")
plt.axvline(0, color="black", linewidth=0.75, linestyle="-")
plt.ylabel("$\dot x(t)\enspace[ms^{-1}]$")
plt.axis("equal")
plt.title("Phase diagram of  $m\ddot x(t) + b\dot x(t) + kx(t)=0$")
plt.subplot(122)
plt.plot(t, pos)
plt.xlabel("time [s]")
plt.ylabel("x(t) [m]")
plt.title("Time domain of $m\ddot x(t) + b\dot x(t) + kx(t)=0$")
plt.tight_layout()
plt.show()
