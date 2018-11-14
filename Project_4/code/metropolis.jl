
function periodic(i, lim, add)
    return (i + lim + add) % lim + 1
end

function montecarlo(N, T, cycles)
    println("starting montecarlo")

    spins = ones((N, N));  # [N, N] array of ones
    E = 0;


    for j in 1:N
        for i in 1:N
            E-= spins[i, j] * (spins[periodic(i, N, -1), j] + spins[i, periodic(j, N, 1)])
        end
    end
end



function main()
    montecarlo(2, 0, 1)
end

main()