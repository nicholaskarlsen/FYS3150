#= Contains the an alternate implementation of the SIRS basic MC algoritim implemented 
in a similar fashion to the python code in SIRS_MCMC.py 

Written for, and tested to work on Julia v1.0.2
=#

using PyPlot

function main()
    for i in 1:100
        t, S, I, R = SIRS_basic_pylike(S0=300, I0=100, R0=0, a=4, b=1, c=0.5, stop_time=10)
        plot(t, S, color="blue", alpha=.1)
        plot(t, I, color="red", alpha=.1)
        plot(t, R, color="green", alpha=.1)
    end
    savefig("../figs/julia_sirs_basic_pylike.png")
end


function SIRS_basic_pylike(;S0::Int64, I0::Int64, R0::Int64, a, b, c, stop_time)
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
    Returns
    -------
        t, S, I, R : Arrays, containing the resulting data
    =#

    N = S0 + I0 + R0 # Population number, for computing step size & probabilities

    Δt = minimum([4.0 / (a*N), 1.0 / (b * N), 1.0 / (c * N)]) 
    t = range(0, stop_time, step=Δt)
    len = length(t)
    # Initialize Arrays for storing population number
    S = Array{Int}(undef, len)
    I = Array{Int}(undef, len)
    R = Array{Int}(undef, len)

    S[1] = S0
    I[1] = I0
    R[1] = R0

    # Begin Monte Carlo Cycle
    for i in 1:(len-1)    
        # Initialize [i+1]'th elements
        S[i + 1] = S[i]
        I[i + 1] = I[i]
        R[i + 1] = R[i]
        # Compute probabilities
        S_I = (a * S[i] * I[i] * Δt) / N   # P(S->I)
        I_R = b * I[i] * Δt                                   # P(I->R)
        R_S = c * R[i] * Δt                                   # P(R->S)
        # Evaluate probabilities against random number, [0, 1)
        if rand(Float64) < S_I
            # Accept & update arrays accordingly
            S[i+1] -= 1
            I[i+1] += 1
        end

        if rand(Float64) < I_R
            I[i+1] -= 1
            R[i+1] += 1
        end

        if rand(Float64) < R_S
            R[i+1] -= 1
            S[i+1] += 1
        end
    end

    return t, S, I, R
end


if PROGRAM_FILE == @__FILE__
    main()
end