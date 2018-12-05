function main(S0, I0, R0, N, tN; a=0, b=0, c=0, d=0, d_I=0, e=0, f=0)
    # note: only args after ; can, and must be called as kwargs, because julia.

    h = tN / N

    S = zeros(N)        # Susceptible
    I = zeros(N)        # Infected
    R = zeros(N)        # Recovered
    Popul = zeros(N)    # Population
    t = zeros(N)

    S[1] = S0
    I[1] = I0
    R[1] = R0
    Popul[1] = S0 + I0 + R0

    dSdt(t, S, I, R) = c * R - a * S * I / (S + I + R)
    dIdt(t, S, I, R) = a * S * I / (S + I + R) - b * I
    dRdt(t, S, I, R) = b * I - c * R

    w = [0, .5*h, .5*h, h]  # Weights of k

    for i in 1:(N-1);
        k_S = zeros(4)
        k_I = zeros(4)
        k_R = zeros(4)

        k_S = dSdt(t[i], S[i], I[i], R[i])
        k_I = dIdt(t[i], S[i], I[i], R[i])
        k_R = dRdt(t[i], S[i], I[i], R[i])
        for j in 2:4
            k_S[j] = dSdt(
                t[j] * w[j], 
                S[j] + w[j] * k_S[j - 1],
                I[j] + w[j] * k_S[j - 1],
                R[j] + w[j] * k_S[j - 1])

            k_S[j] = dSdt(
                t[j] * w[j], 
                S[j] + w[j] * k_S[j - 1],
                I[j] + w[j] * k_S[j - 1],
                R[j] + w[j] * k_S[j - 1])

            k_S[j] = dSdt(
                t[j] * w[j], 
                S[j] + w[j] * k_S[j - 1],
                I[j] + w[j] * k_S[j - 1],
                R[j] + w[j] * k_S[j - 1])
        end

        S[i + 1] = S[i] + (h / 6.0) * (k_S[1] + 2k_S[2] + 2k_S[3] + k_S[4])
        I[i + 1] = I[i] + (h / 6.0) * (k_I[1] + 2k_I[2] + 2k_I[3] + k_I[4])
        R[i + 1] = R[i] + (h / 6.0) * (k_R[1] + 2k_R[2] + 2k_R[3] + k_R[4])
    end
end



main(300, 100, 0, 100, 10, a=4, b=1, c=0.5)