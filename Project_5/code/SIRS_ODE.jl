function RK4(;t, S, I, R, h, diffEq::Function)
        k1 = diffEq(
            t,
            S,
            I,
            R
        )
        k2 = diffEq(
            t + h / 2.0,
            S + h / 2.0 * k1,
            I + h / 2.0 * k1,
            R + h / 2.0 * k1
        )
        k3 = diffEq(
            t + h / 2.0,
            S + h / 2.0 * k2,
            I + h / 2.0 * k2,
            R + h / 2.0 * k2
        )
        k4 = diffEq(
            t + h,
            S + h * k3,
            I + h * k3,
            R + h * k3
        )

        current = [S, I, R]
        next = zeros(3)

        for i in 1:3
            next[i] = current[i] + (h / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
        end

        return next
end





function SIRS_Basic(;S0, I0, R0, a, b, c,
        N, stop)
     
    t = range(0, stop=1, length=N)
    h = stop / N

    S = zeros(N)
    I = zeros(N)
    R = zeros(N)

    function diffeqs(t, S_, I_, R_)
        num_people = S_ + I_ + R_
        dSdt = c * R_ - a * S_ * I_ / num_people   
        dIdt = a * S_ * I_ / num_people - b * I_
        dRdt = b * I_ - c * R_
        return [dSdt, dIdt, dRdt]
    end

    for i in 1:N-1
        RK4(t=0., S=0., I=0., R=0., h=0., diffEq=diffeqs)
    end




    return S
end


println(SIRS_Basic(S0=300., I0=100., R0=0., a=4., b=1., c=.5, N=10, stop=15.))