# tested to work on Julia v1.0.2

function SIRS_vitdyn(;S0::Int64, I0::Int64, R0::Int64, a, b, c, d, d_I, e, stop_time)
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

        # Ensure we dont get negative population groups
        for list in [S, I, R]
            if list[i] < 0
                list[i] = 0
            end
        end

        # Update current population after population check due to deaths etc.
        N[i] = S[i] + I[i] + R[i]

        Δt = minimum([4.0 / (a*N[i]), 1.0 / (b * N[i]), 1.0 / (c * N[i])]) 
        # Initialize [i + 1]'th elements
        append!(S, S[i])
        append!(I, I[i])
        append!(R, R[i])
        append!(N, N[i])
        append!(t, t[i] + Δt)
        # Compute transition probabilities
        S_I = (a * S[i] * I[i] * Δt) / (N[i])   # P(S->I)
        I_R = b * I[i] * Δt                     # P(I->R)
        R_S = c * R[i] * Δt                     # P(R->S)
        # Probability of dying
        S_D = d * S[i] * Δt                     # P(S->DEAD) 
        I_D = (d + d_I) * I[i] * Δt             # P(I->DEAD)
        R_D = d * R[i] * Δt                     # P(R->DEAD)
        # Probability of birth
        B_S = e * N[i]                          # P(B->S)

        # Evaluate probabilities against random number, [0, 1)
        # NOTE: Mersenne Twister
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

        if rand(Float64) < I_D  # Death in I group
            if I[i + 1] > 0
                I[i + 1] -= 1
            end
        end

        if rand(Float64) < R_D  # Death in R group
            if R[i + 1] > 0
                R[i + 1] -= 1
            end
        end

        if rand(Float64) < B_S  # Birth to S group
            S[i + 1] += 1
        end

        i+=1
    end

    return t, S, I, R

end

