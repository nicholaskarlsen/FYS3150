# This script solves the 2-body problem using a simple approach in two dimensions
# Meant mostly to familiarize myself with the system, hence the spaghetti code.
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from main import figsetup

_G = 4 * np.pi**2    # Gravitational Constant [AU^3 yr^-2 M_sun^-1]


def gravity(pos, mass):
    return - pos * _G * mass / np.linalg.norm(pos) ** 3


def velocityverlet(diffeq):
    for i in xrange(N - 1):
        a1 = diffeq(pos[i], mass=1)   # mass of sun = 1
        pos[i + 1] = pos[i] + vel[i] * h + 0.5 * a1 * h**2
        a2 = diffeq(pos[i + 1], mass=1)
        vel[i + 1] = vel[i] + 0.5 * h * (a1 + a2)


def eulercromer(diffeq):
    for i in xrange(N - 1):
        acc = diffeq(pos[i], mass=1)   # mass of sun = 1
        vel[i + 1] = vel[i] + h * acc
        pos[i + 1] = pos[i] + h * vel[i + 1]


if __name__ == '__main__':

    tn = 10          # Timespan (1 year)
    N = int(1e6)    # Number of integration points
    h = tn / N      # Step size
    print "Solving 2 body system using Euler cromer & velocityverlet, h=", h

    mass = 3.00348959632E-6  # Mass of earth [Solar mass]

    time = np.linspace(0, tn, N)
    pos = np.zeros([N, 2])    # (t)(x, y)
    vel = np.zeros([N, 2])    # (t)(x, y)
    pos[0] = [1, 0]           # [AU]
    vel[0] = [0, 2 * np.pi]   # [AU/yr]

    velocityverlet(gravity)
    E = np.zeros(N)
    U = np.zeros(N)
    for i in range(N):
        E[i] = 0.5 * mass * np.linalg.norm(vel[i])**2
        U[i] = - mass * _G / np.linalg.norm(pos[i])

    plt.figure(figsize=(5, 3))
    plt.plot(time, U + E)
    figsetup(" ", "Time [Years]", "Energy", "exb_energy_verlet", show=False)

    plt.figure(figsize=(5, 5))
    plt.plot(pos[:, 0], pos[:, 1])
    figsetup(" ", "x [AU]", "y [AU]", "exb_orbit_verlet", show=False)

    # re setting data
    time = np.linspace(0, tn, N)
    pos = np.zeros([N, 2])    # (t)(x, y)
    vel = np.zeros([N, 2])    # (t)(x, y)
    pos[0] = [1, 0]           # [AU]
    vel[0] = [0, 2 * np.pi]   # [AU/yr]

    eulercromer(gravity)
    E = np.zeros(N)
    U = np.zeros(N)
    for i in range(N):
        E[i] = 0.5 * mass * np.linalg.norm(vel[i])**2
        U[i] = - mass * _G / np.linalg.norm(pos[i])

    plt.figure(figsize=(5, 3))
    plt.plot(time, U + E)
    figsetup(" ", "Time [Years]", "Energy", "exb_energy_euler", show=False)

    plt.figure(figsize=(5, 5))
    plt.plot(pos[:, 0], pos[:, 1])
    figsetup(" ", "x [AU]", "y [AU]", "exb_orbit_euler", show=False)
