# tested to work on Julia v1.0.2

function main()
    # for checking functionality & getting error messages.
    # Note : Import guard @ bottom of file
    SIRS_svar(S0=300, I0=100, R0=0, a0=4, A=1, omega=1, b=1, c=0.5, stop_time=10, trials=100)
end


function SIRS_svar(;S0::Int64, I0::Int64, R0::Int64, a0, A, omega, b, c, stop_time, trials)
    #= Simulates the SIRS model with seasonal variation

    NOTE : no unicode omega, python says no.

    Parameters
    ----------
        S0        : Initial no. susceptible
        I0        : Initial no. infected
        R0        : Initial no. recovered
        a0        : base rate of transmission (S -> I)
        A         : maximum deviation from base rate of transmission 
        omega     : frequency of deviation from base rate of transmission
        b         : rate of recovery (I -> R)
        c         : rate of immunity loss (R -> S)
        stop_time : Time at which the loop is terminated
        trials    : No. Monte-carlo cycles to run
    Returns
    -------
        t, S, I, R : Arrays, containing the resulting data
    =#

    N = S0 + I0 + R0  # Population must be convserved in this model

    # Pre-compute a, t, Δt. -> same for each solution & required to compute length of S, I, R
    a = Float64[]
    t = Float64[]
    Δt = Float64[]
    append!(t, 0)
    i = 1
    while t[i] <= stop_time;
        append!(a, A * cos(omega * t[i]) + a0)
        append!(Δt, minimum([4.0 / (a[i]*N), 1.0 / (b * N), 1.0 / (c * N)]))
        append!(t, t[i] + Δt[i])
        i += 1
    end

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

    # Begin Monte Carlo Cycles
    for j in 1:trials
        for i in 1:(len - 1)
            if i > 1000000 # TODO : remove this once bugs are fixed
                println("Broke early")
                break
            end
            # Initialize [i + 1]'th elements
            S[j, i + 1] = S[j, i]
            I[j, i + 1] = I[j, i]
            R[j, i + 1] = R[j, i]
            # Compute transition probabilities
            S_I = (a[i] * S[j,i] * I[j,i] * Δt[i]) / (N)   # P(S->I)
            I_R = b * I[j, i] * Δt[i]                     # P(I->R)
            R_S = c * R[j, i] * Δt[i]                     # P(R->S)
            # Evaluate probabilities against random number, [0, 1)

            if rand(Float64) < S_I  # Transition S -> I
                # Accept & update arrays accordingly
                if S[j, i + 1] > 0
                    S[j, i + 1] -= 1
                    I[j, i + 1] += 1
                end
            end

            if rand(Float64) < I_R  # Transition I -> R
                if I[j, i + 1] > 0
                    I[j, i + 1] -= 1
                    R[j, i + 1] += 1
                end
            end

            if rand(Float64) < R_S  # Transition R -> S
                if R[j, i + 1] > 0
                    R[j, i + 1] -= 1
                    S[j, i + 1] += 1
                end
            end
        end
    end

    return t, S, I, R

end


if PROGRAM_FILE == @__FILE__ 
    main()
end
