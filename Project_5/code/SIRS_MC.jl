#=
Contains all different versions of the SIRS MC implementations.

Written for, and tested to work on Julia v1.0.2.
=#


function SIRS_basic(;S0::Int64, I0::Int64, R0::Int64, a, b, c, stop_time)
    #= Simulates the simplest case of the SIRS model with a Monte Carlo method

    Note: This method is inefficient by design. It was used to test how i might implement
          the algorithm with dynamic time-steps as required in some of the later systems.
          a more efficient Julia implementation is found in SIRS_Basic.jl

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

    # Begin Monte Carlo Cycle
    i=1 # Initialize itteration variable
    while t[i] <= stop_time
        Δt = minimum([4.0 / (a*N[i]), 1.0 / (b * N[i]), 1.0 / (c * N[i])]) 
        # Initialize [i+1]'th elements
        append!(S, S[i])
        append!(I, I[i])
        append!(R, R[i])
        append!(N, N[i])
        append!(t, t[i] + Δt)
        # Compute probabilities
        S_I = (a * S[i] * I[i] * Δt) / (N[i])   # P(S->I)
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

        i+=1
    end

    return t, S, I, R
end



function SIRS_vitdyn(;S0::Int64, I0::Int64, R0::Int64, a, b, c, d, d_I, e, stop_time)
    #= Simulates the SIRS model with added vital dynamics, birth and death, with
    different death-rate for the infected population using a Monte-Carlo method.

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
        # Initialize [i+1]'th elements
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
            # Accept & update arrays accordingly
            S[i+1] -= 1
            I[i+1] += 1
        end

        if rand(Float64) < I_R  # Transition I -> R
            I[i+1] -= 1
            R[i+1] += 1
        end

        if rand(Float64) < R_S  # Transition R -> S
            R[i+1] -= 1
            S[i+1] += 1
        end

        if rand(Float64) < S_D  # Death in S group
            S[i+1] -= 1
        end

        if rand(Float64) < I_D  # Death in I group
            I[i+1] -= 1
        end

        if rand(Float64) < R_D  # Death in R group
            R[i+1] -= 1
        end

        if rand(Float64) < B_S  # Birth to S group
            S[i+1] += 1
        end

        i+=1
    end

    return t, S, I, R

end


function SIRS_svar(;S0::Int64, I0::Int64, R0::Int64, a0, A, omega, b, c, stop_time)
    #= Simulates the SIRS model with seasonal variation

    Parameters
    ----------
        S0        : Initial no. susceptible
        I0        : Initial no. infected
        R0        : Initial no. recovered
        a0        : base rate of transmission (S -> I)
        A         : maximum deviation from base rate of transmission 
        ω         : frequency of deviation from base rate of transmission
        b         : rate of recovery (I -> R)
        c         : rate of immunity loss (R -> S)
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

    # Begin Monte Carlo Cycle
    i=1 # Initialize iteration variable
    while t[i] <= stop_time
        if i > 1000000
            println("Broke early")
            break
        end
        # Ensure that we don't get negative population
        for list in [S, I, R]
            if list[i] < 0
                list[i] = 0
            end
        end
        a = A * cos(omega * t[i]) + a0
        Δt = minimum([4.0 / (a*N[i]), 1.0 / (b * N[i]), 1.0 / (c * N[i])]) 
        # Initialize [i+1]'th elements
        append!(S, S[i])
        append!(I, I[i])
        append!(R, R[i])
        append!(N, N[i])
        append!(t, t[i] + Δt)
        # Compute transition probabilities
        S_I = (a * S[i] * I[i] * Δt) / (N[i])   # P(S->I)
        I_R = b * I[i] * Δt                     # P(I->R)
        R_S = c * R[i] * Δt                     # P(R->S)
        # Evaluate probabilities against random number, [0, 1)
        if rand(Float64) < S_I  # Transition S -> I
            # Accept & update arrays accordingly
            S[i+1] -= 1
            I[i+1] += 1
        end

        if rand(Float64) < I_R  # Transition I -> R
            I[i+1] -= 1
            R[i+1] += 1
        end

        if rand(Float64) < R_S  # Transition R -> S
            R[i+1] -= 1
            S[i+1] += 1
        end
        i+=1
    end

    return t, S, I, R

end


function SIRS_vax(;S0::Int64, I0::Int64, R0::Int64, a, b, c, f, stop_time)
    #= Simulates the SIRS model with the added property of vaccinations!

    Parameters
    ----------
        S0        : Initial no. susceptibles
        I0        : Initial no. infected
        R0        : Initial no. recovered
        a         : rate of transmission (S -> I)
        b         : rate of recovery (I -> R)
        c         : rate of immunity loss (R -> S)
        f         : rate of vaccination (S -> R)
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

    append!(S, S0)
    append!(I, I0)
    append!(R, R0)
    append!(N, S0 + I0 + R0)

    # Initialize time array
    t = Float64[]
    Δt = minimum([4.0 / (a*N[1]), 1.0 / (b * N[1]), 1.0 / (c * N[1])]) 
    append!(t, Δt)

    # Begin Monte Carlo Cycle
    i=1 # Initialize itteration variable
    while t[i] <= stop_time
        Δt = minimum([4.0 / (a*N[i]), 1.0 / (b * N[i]), 1.0 / (c * N[i]), f]) 
        # Initialize [i+1]'th elements
        append!(S, S[i])
        append!(I, I[i])
        append!(R, R[i])
        append!(N, N[i])
        append!(t, t[i] + Δt)
        # Compute probabilities
        S_I = (a * S[i] * I[i] * Δt) / (N[i])   # P(S->I)
        I_R = b * I[i] * Δt                     # P(I->R)
        R_S = c * R[i] * Δt                     # P(R->S)
        S_R = f * Δt                            # P(S->R)
        # Evaluate probabilities against random number, [0, 1]
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

        if rand(Float64) < S_R
            S[i+1] -= 1
            R[i+1] += 1
        end

        i+=1
    end

    return t, S, I, R
end


function main()
    #function body
end

if PROGRAM_FILE == @__FILE__
    main()
end
