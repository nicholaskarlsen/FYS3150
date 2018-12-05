using NPZ   # Allows for reading and writing numpy arrays

function ODE_SIRS(;S0, I0, R0, a, b, c, N, tN)
    #=
    Note: ALL kwargs for clarity when calling the code (-> needs to be called 
    as such, because julia does not allow to pass an arg as a kwarg.)
    =#

    t = range(0, stop = tN, length= N)
    h = tN / N

    S = zeros(N)
    I = zeros(N)
    R = zeros(N)
    Population = zeros(N)

    S[1] = S0
    I[1] = I0
    R[1] = R0

    Population[1] = sum([S0, I0, R0]) 

    for i in 1:(N-1)
        dSdt = c * R[i] - a * S[i] * I[i] / Population[i]
        dIdt = a * S[i] * I[i] / Population[i] - b * I[i]
        dRdt = b * I[i] - c * R[i]

        S[i + 1] = S[i] + h * dSdt
        I[i + 1] = I[i] + h * dIdt
        R[i + 1] = R[i] + h * dRdt
    end
    return S, I, R, t
end


function main()
    S, I, R, t = ODE_SIRS(S0=300, I0=100, R0=0, a=4, b=1, c=0.5, N=1000, tN=10)
    npzwrite("../data/S.npy", S)
    npzwrite("../data/I.npy", I)
    npzwrite("../data/R.npy", R)
end


if PROGRAM_FILE == @__FILE__
    main()
end