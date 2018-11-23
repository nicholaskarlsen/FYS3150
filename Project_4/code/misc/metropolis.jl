
function periodic(i, lim, add)
    #=
    Accounts for periodic bound. conditions when choosing lattice indices.
    (+1 because arrays start at 1, in julia)
    =#
    return (i + lim + add) % lim + 1
end

function latticeEnergy(lattice)
    # Computes the energy of square spin lattice
    energy = 0                      # initialize energy variable
    N = length(lattice[1])  # length of first row
    for j in 1:N      # Compute sum of all couplings with periodic bounds.
        for i in 1:N
            energy -= lattice[i, j] * (lattice[periodic(i, N, -1), j] + lattice[i, periodic(j, N, 1)])
        end
    end
    return energy
end


function montecarlo(N, T, cycles)
    println("starting montecarlo")

    spins = ones((N, N)) # [N, N] array of ones
    E = lattinceEnergy(spins)
    print(E)
    
end


function testFunc()
    # Expect energy of   
    testLattice = ones((2,2))
    println(testLattice)
    testEnergy = latticeEnergy(testLattice)
    println(testEnergy)

end




function main()
    testFunc()
end

main()