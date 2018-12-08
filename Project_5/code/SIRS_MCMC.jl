using NPZ  # Used for saving arrays as numpy array.

function main(;trials::Int)
    S0 = 300
    I0 = 100
    R0 = 0
    a = 4
    c = 0.5
    stop_time = 10
    for i in 1:4
        for j in 1:trials
            t, S, I, R = SIRS_basic(S0=S0, I0=I0, R0=R0, a=a, b=i, c=c,
                stop_time=stop_time)
            #= In some of the MCMC simulations, where the population N varies over time,
            the step sizes and by extension, length of arrays vary. So can not store
            as a square array in a single data file. =#
            npzwrite("../data/prob_b/b_$i/S_$j.npy", S)
            npzwrite("../data/prob_b/b_$i/I_$j.npy", I)
            npzwrite("../data/prob_b/b_$i/R_$j.npy", R)
            npzwrite("../data/prob_b/b_$i/t_$j.npy", t)
        end
    end
end



function SIRS_basic(;S0::Int64, I0::Int64, R0::Int64, a, b, c, stop_time)
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



function SIRS_vitdyn(;S0::Int64, I0::Int64, R0::Int64, a, b, c, e, d, d_I, stop_time)
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

        if rand(Float64) < S_D
            S[i+1] -= 1
        end

        if rand(Float64) < I_D
            I[i+1] -= 1
        end

        if rand(Float64) < R_D
            R[i+1] -= 1
        end

        i+=1
    end

    return t, S, I, R

end


if PROGRAM_FILE == @__FILE__
    main(trials=1000)
end