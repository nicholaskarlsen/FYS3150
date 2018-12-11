# tested to work on Julia v1.0.2

using PyPlot

function main()
    # Test functionality, not report material
    # Note : import guard @ bottom of file
    for i in 1:100
        t, S, I, R = SIRS_vitdyn(S0=300, I0=0, R0=0, a=4, b=1, c=0.5, d=4, d_I=0, e=0.1, stop_time=1,
                                limit=100000)
        plot(t, S, color="blue", alpha=.1)
        plot(t, I, color="red", alpha=.1)
        plot(t, R, color="green", alpha=.1)
        plot(t, S + I + R, color="black", alpha=.1)
    end
    savefig("../figs/julia_sirs_vitdyn3.png")
end

function SIRS_vitdyn(;S0::Int64, I0::Int64, R0::Int64, a, b, c, d, d_I, e, stop_time, 
                     limit=100000)
    #= Simulates ONE cycle of the SIRS mode with vital dynamics.
    Elaborated in report, but tl;dr differing number of steps per cycle
    => need to interpolate to produce averages. this is done in main.py

    Parameters
    ----------
        S0        : Initial no. susceptible
        I0        : Initial no. infected
        R0        : Initial no. recovered
        a         : rate of transmission (S -> I)
        b         : rate of recovery (I -> R)
        c         : rate of immunity loss (R -> S)
        d         : death rate (Any Group -> Dead)
        d_I       : death rate of the infected, additive to d. (Any group -> )
        e         : Birth Rate (-> S)
        stop_time : Time at which the loop is terminated
        limit     : limit number of loops incase N grows too large (causing Δt to become smaller)
                    Note: when producing averages, RAISE limit or lower stop_time untill 
                    there are no early breaks!!
    Returns
    -------
        t, S, I, R : Arrays, containing the resulting data
    =#

    # Initialize Arrays for storing population number
    S = Int64[]
    I = Int64[]
    R = Int64[]
    N = Int64[]
    t = Float64[]

    append!(S, S0)
    append!(I, I0)
    append!(R, R0)
    append!(N, S0 + I0 + R0)
    append!(t, 0)

    # Initialize time array
    # Begin Monte Carlo Cycle
    i=1 # Initialize itteration variable

    while t[i] <= stop_time
        # birthrate may cause N -> large, => Δt -> small. Terminate at limit.
        if i > limit
            println("Broke SIRS_MC_vitdyn early: i > $limit")
            break
        end


        Δt = minimum([4.0 / (a*N[i]), 1.0 / (b * N[i]), 1.0 / (c * N[i]), 
                    1.0 / (e * N[i]), (1.0 / d * N[i]), 1.0 / ((d + d_I) * N[i])]) 
        # Initialize [i + 1]'th elements
        append!(S, 0)
        append!(I, 0)
        append!(R, 0)
        append!(t, 0)
        append!(N, 0)
        # Appending directly caused issues with scope inside of IF statements???
        S[i + 1] = S[i]
        I[i + 1] = I[i]
        R[i + 1] = R[i]
        t[i + 1] = t[i] + Δt

        # Compute transition probabilities
        S_I = (a * S[i] * I[i] * Δt) / (N[i])   # P(S->I)
        I_R = b * I[i] * Δt                     # P(I->R)
        R_S = c * R[i] * Δt                     # P(R->S)
        # Probability of dying
        S_D = d * S[i] * Δt                     # P(S->DEAD) 
        I_D = (d + d_I) * I[i] * Δt             # P(I->DEAD)
        R_D = d * R[i] * Δt                     # P(R->DEAD)
        # Probability of birth
        B_S = e * N[i] * Δt                     # P(BIRTH->S)

        # Evaluate probabilities against random number, (0, 1)
        if rand(Float64) < S_I  # Transition S -> I  
            if S[i + 1] > 0
                # Accept & update arrays accordingly
                S[i + 1] -= 1
                I[i + 1] += 1
            end
        end

        if rand(Float64) < I_R  # Transition I -> R
            if I[i + 1] > 0
                I[i + 1] -= 1
                R[i + 1] += 1
            end
        end

        if rand(Float64) < R_S  # Transition R -> S
            if R[i + 1] > 0
                R[i + 1] -= 1
                S[i + 1] += 1
            end
        end

        if rand(Float64) < S_D  # Death in S group
            if S[i + 1] > 0
                S[i + 1] -= 1
            end
        end

        if rand(Float64) < I_D; # Death in I group
            if I[i + 1] > 0;
                I[i + 1] -= 1
            end
        end

        if rand(Float64) < R_D;  # Death in R group
            if R[i + 1] > 0
                R[i + 1] -= 1
            end
        end

        if rand(Float64) < B_S;  # Birth to S group
            S[i + 1] += 1
        end

        N[i+1] = S[i+1] + I[i+1] + R[i+1]

        i+=1
    end

    return t, S, I, R

end

# Import guard
if PROGRAM_FILE == @__FILE__
    main()
end