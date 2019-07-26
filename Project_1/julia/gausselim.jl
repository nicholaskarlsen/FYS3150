#=
Written for Julia v1.1.0
@author Nicholas Karlsen
=#

function source(x)
    100 * exp(-100*x)
end

function general(n, a, b, c)
    #=
    Solves for lenght n vector u in system Au = f(x) where A, [nxn] sparse matrix
    and f(x) lenght n vector, function of x.

    Parameters:
        n - Number of steps
        a, b, c - diagonal entries of a sparse matrix A. Below diagonal, diagonal and above diagonal
            respectively. 
    Returns:
        x - vector
        u - eigenvalue
    =#

    x = linspace(0, 1.0, n)
    h = 1.0 / (n - 1)

    f = source(x)
    u = zeros(n)
    
    # TODO: UPDATE INDICES FROM 0 Based to 1 Based 
    # Forward Substitution 
    for i in 2:n-1
        b[i] = b[i] - ((a[i] * c[i - 1]) / b[i - 1])
        f[i] = f[i] - ((a[i] * f[i - 1]) / b[i - 1])
    end
    u[-2] = f[-2] / b[-2]

    # Backward substitution
    for i in n-2:0
        u[i] = (f[i] - c[i] * u[i + 1]) / b[i]
    end

    return x, u
end
