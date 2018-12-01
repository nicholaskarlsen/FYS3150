from __future__ import division


def solve(S0, I0, R0, a, b, c, N, tN):

    t = np.linspace(0, tN, N)

    S = np.zeros(N)  # Suceptible
    I = np.zeros(N)  # Infected
    R = np.zeros(N)  # Recovered
    population = np.zeros(N)
    S[0] = S0
    I[0] = I0
    R[0] = R0
    population[0] = sum([S[0], I[0], R[0]])

    
