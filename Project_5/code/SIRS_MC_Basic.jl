# tested to work on Julia v1.0.2

using PyPlot

function main()
    # check functionality, not report material
    println("Running main()")
    ntr = 100
    t, S, I, R = SIRS_basic(S0=300, I0=100, R0=0, a=4, b=1, c=0.5, stop_time=10, trials=ntr)
    for i in 1:ntr
        plot(t, S[i, :], color="blue", alpha=.1)
        plot(t, I[i, :], color="red", alpha=.1)
        plot(t, R[i, :], color="green", alpha=.1)
    end
    savefig("../figs/julia_sirs_basic_pylike.png")
end


function SIRS_basic(;S0::Int64, I0::Int64, R0::Int64, a, b, c, stop_time, trials::Int64)
    #= Simulates the simplest case of the SIRS model with a Monte Carlo method

    Note: Written to mirror the program in SIRS_MCMC.py

    Parameters
    ----------

        S0        : Initial no. susceptibles
        I0        : Initial no. infected
        R0        : Initial no. recovered
        a         : rate of transmission (S -> I)
        b         : rate of recovery (I -> R)
        c         : rate of immunity loss (R -> S)
        stop_time : Time at which the loop is terminated
        trials    : No. MC trials
    Returns
    -------
        t, S, I, R : Arrays, containing the resulting data
    =#

    N = S0 + I0 + R0 # Population number, for computing step size & probabilities

    Δt = minimum([4.0 / (a*N), 1.0 / (b * N), 1.0 / (c * N)]) 
    t = range(0, stop_time, step=Δt)
    len = length(t)
    # Initialize Arrays for storing population number
    S = Array{Int}(undef, trials, len)
    I = Array{Int}(undef, trials, len)
    R = Array{Int}(undef, trials, len)

    for i in 1:trials
        S[i, 1] = S0
        I[i, 1] = I0
        R[i, 1] = R0
    end


    # Begin Monte Carlo Cycle
    for j in 1:trials
        for i in 1:(len-1)    
            # Initialize [i + 1]'th elements
            S[j, i + 1] = S[j, i]
            I[j, i + 1] = I[j, i]
            R[j, i + 1] = R[j, i]
            # Compute probabilities
            S_I = (a * S[j, i] * I[j, i] * Δt) / N   # P(S->I)
            I_R = b * I[j, i] * Δt                # P(I->R)
            R_S = c * R[j, i] * Δt                # P(R->S)
            # Evaluate probabilities against random number, [0, 1)
            if rand(Float64) < S_I
                # Accept & update arrays accordingly
                if S[j, i + 1] > 0 
                    S[j, i + 1] -= 1
                    I[j, i + 1] += 1
                end
            end

            if rand(Float64) < I_R
                if I[j, i + 1] > 0
                    I[j, i + 1] -= 1
                    R[j, i + 1] += 1
                end
            end

            if rand(Float64) < R_S
                if R[j, i + 1] > 0
                    R[j, i + 1] -= 1
                    S[j, i + 1] += 1
                end
            end
        end
    end

    return t, S, I, R
end
