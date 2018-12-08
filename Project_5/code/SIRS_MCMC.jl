
function SIRS_basic(;S0::Int64, I0::Int64, R0::Int64, a, b, c, trials::Int64, stop_time)

    S = Int64[]
    I = Int64[]
    R = Int64[]
    N = Int64[]

    for i in 1:trials
        append!(S, S0)
        append!(I, I0)
        append!(R, R0)
        append!(N, S0 + I0 + R0)
    end

    Δt = minimum([4.0 / (a*N[1]), 1.0 / (b * N[1]), 1.0 / (c * N[1])]) 
    t = range(0, stop=stop_time, step=Δt)

    for i in 1:length(t)-1
        tmp = zeros(trials)
        S = cat(dims=2, S, tmp)   
        I = cat(dims=2, I, tmp)   
        R = cat(dims=2, R, tmp)   
        for j in 1:trials
            println(S[i][j])
        end
    end

    return size(S)

end





println(SIRS_basic(S0=300, I0=100, R0=0, a=4, b=1, c=0.5, trials=10, stop_time=1))
