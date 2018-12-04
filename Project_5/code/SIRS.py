from __future__ import division  # No on expects the integer division
import numpy as np
import matplotlib.pyplot as plt


class SIRS:
    def __init__(self, S0, I0, R0, N, tN, a=0, b=0, c=0, d=0, d_I=0, e=0, f=0):
        """
        Solves the SIRS model numerically as a differential equation using Runge-Kutta 4
        or as a Monte-Carlo simulation using a Markov Chain Monte-Carlo method.

        (Originally planned to contian both ODE & MCMC Solvers, but turned out
        to be impractical, hence 1 Class-based solution and one Function-based)

        Parameters
        ----------
        S0: Initial number of Suceptible people
        I0:               ... Infected
        R0:               ... Recovered people
        a:  Rate of Transmission                    [1/Time]
        b:      ... Recovery                        [1/Time]
        c:      ... Immunity loss                   [1/Time]
        N: Number of integration points
        tN: Simulation time                         [Time]
        """

        self.N = N
        self.tN = tN
        self.h = self.tN / self.N

        # Initialize arrays for storing data
        self.S = np.zeros(N)
        self.I = np.zeros(N)
        self.R = np.zeros(N)
        self.t = np.linspace(0, tN, N)
        # Rate of; Transmission, Recovery and Immunity
        self.a, self.b, self.c = a, b, c

        # Initial Conditions
        self.S[0] = S0
        self.I[0] = I0
        self.R[0] = R0

        return

    def diffeqs(self, t, S, I, R):
        """
        Compute differentials for the SIRS system at a time t. Kept as a nested function
        to access the local namespace for a, b, c constants.
        ("Wasted" FLOPS by multiplying by variables which =0, but requires less methods.
        If runtime becomes issue, add new method with only the speciffic trains)

        Parameters
        ----------
        t: Time
        S: Suceptibles
        I: Infected
        R: Recovered

        Returns
        -------
        Numpy array
            (S', I', R')
        """
        dSdt = self.c * R - self.a * S * I / (S + I + R)
        dIdt = self.a * S * I / (S + I + R) - self.b * I
        dRdt = self.b * I - self.c * R

        return np.array([dSdt, dIdt, dRdt])

    def euler_fw(self, i, diffEq):
        """
        Computes time-step using the Euler-Forward method
        (used to test implementation of RK4. if both yield same result -> things 
        are working correclty, probably.)

        Parameters
        ----------
        i : Current index
        diffEq : Function, must take arguements as following func(t, S, I, R)

        Returns
        -------
        Numpy Array
            (Next S, Next I, Next R)
        """

        current_val = np.array([self.S[i], self.I[i], self.R[i]])
        next_val = current_val + self.h * diffEq(self.t[i], self.S[i], self.I[i], self.R[i])

        return next_val

    def rk4(self, i, diffEq):
        """
        Computes time-step using the 4th order Runge Kutta method

        Parameters
        ----------
        i : Current index
        diffEq : Function, must take arguements as following func(t, S, I, R)

        Returns
        -------
        Numpy Array
            (Next S, Next I, Next R)
        """
        current_val = np.array([self.S[i], self.I[i], self.R[i]])

        k2 = np.zeros(3)
        k3 = np.zeros(3)
        k4 = np.zeros(3)

        # returns size=(3) array
        k1 = diffEq(
            self.t[i],
            self.S[i],
            self.I[i],
            self.R[i]
        )
        k2 = diffEq(
            self.t[i] + self.h / 2.0,
            self.S[i] + self.h / 2.0 * k1[0],
            self.I[i] + self.h / 2.0 * k1[1],
            self.R[i] + self.h / 2.0 * k1[2]
        )
        k3 = diffEq(
            self.t[i] + self.h / 2.0,
            self.S[i] + self.h / 2.0 * k2[0],
            self.I[i] + self.h / 2.0 * k2[1],
            self.R[i] + self.h / 2.0 * k2[2]
        )
        k4 = diffEq(
            self.t[i] + self.h,
            self.S[i] + self.h * k3[0],
            self.I[i] + self.h * k3[1],
            self.R[i] + self.h * k3[2]
        )

        next_val = current_val + (self.h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

        return next_val

    def ODESolve(self, solverMethod, diffEq):
        """
        Parameters
        ----------
        solverMethod : Method which computes the next time step with input (i, diffEq)
                    and which outputs a len=3 array in the form (S, I, R)
        diffEq : Method which outputs derivative of system with input (t, S, I, R) and
                and outputs satisfying requirements for solverMethod.

        Returns
        -------
        N/A

        """
        for i in range(self.N - 1):
            self.S[i + 1], self.I[i + 1], self.R[i + 1] = solverMethod(i, diffEq)

        return

    def getData(self):
        return self.t, self.S, self.I, self.R

    def plotSIR(self):  # Added SIR at end to avoid confusion with matplot
        plt.plot(self.t, self.S, label="Suceptible", color="Blue")
        plt.plot(self.t, self.I, label="Infected", color="Red")
        plt.plot(self.t, self.R, label="Recovered", color="Green")

        plt.xlabel("Time")
        plt.ylabel("No. People")
        plt.legend(loc="best")
        plt.text(0.5, 0.5, "Population = %i" % sum([self.S[0], self.I[0], self.R[0]]))
        plt.show()

        return

    def __repr__(self):
        """
        A more informative __repr__ method. Useful when working with multiple instances
        in an interactive enviorment.
        """

        return "<SIRS Object: a=%.2f, b=%.2f, c=%.2f>" % (self.a, self.b, self.c)


def main():
    sys1 = SIRS(S0=300, I0=100, R0=0, a=4, b=4, c=0.5, N=100, tN=10)
    sys1.ODESolve(sys1.rk4, sys1.diffeqs)
    sys1.plotSIR()
    return


if __name__ == '__main__':
    main()
