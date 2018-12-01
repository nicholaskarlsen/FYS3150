# Julia Version 1.0.2

function ODE_SIRS(S0, I0, R0, a, b, c, N, tN)
    initPopulation = S0 + I0 + R0

    t = LinRange(0, tN, N)
end

function main()
   ODE_SIRS(300, 100, 0,4, 1, 0.5, 1e3, 10) 
end


if PROGRAM_FILE == @__FILE__
    main()
end
