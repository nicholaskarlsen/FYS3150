import numpy as np


class SIRS:
    def __init__(self, S0, I0, R0, a, b, c, N, tN):
        """
        Solves the SIRS model numerically as a differential equation using Runge-Kutta 4
        or as a Monte-Carlo simulation using a Markov Chain Monte-Carlo method.

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
        # Rate of; Transmission, Recovery and Immunity
        self.a, self.b, self.c = a, b, c

    def diffEqs(self, t, S, I, R):
        """
        Compute differentials for the SIRS system at a time t. Kept as a nested function
        to access the local namespace for a, b, c constants & allow for the usage of
        numba.jit.

        Parameters
        ----------
        t: Time
        S: Suceptibles
        I: Infected
        R: Recovered

        Returns
        -------
        Numpy array
            Array filled with differentials governing the system at time t.
        """
        dSdt = self.c * R - self.a * S * I / (S + I + R)
        dIdt = self.a * S * I / (S + I + R) - self.b * I
        dRdt = self.b * I - self.c * R

        return np.array([dSdt, dIdt, dRdt])

    def rk4(self, i, diffEq):
        """
        Performs a single step in a differential equation, or set of coupled
        differential equations. If the input function returns the derivative of multiple
        paramaters, the function will compute and return the next step for all of the
        parameters, negating the need for multiple function calls (and enjoying the
        slight performance increase of working with numpy arrays over loops)

        Parameters
        ----------
        i : Current index
        diffEq : Function, must take arguements as following func(t, S, I, R)

        Returns
        -------
        Number OR Numpy Array
            Output will match the output of the diffEq function.
        """
        current_val = np.array([self.S[i], self.I[i], self.R[i]])
        k1 = diffEq(
            self.t[i],
            self.S[i],
            self.I[i],
            self.R[i]
        )
        k2 = diffEq(
            self.t[i] + self.h / 2.0,
            self.S[i] + self.h / 2.0 * k1,
            self.I[i] + self.h / 2.0 * k1,
            self.R[i] + self.h / 2.0 * k1
        )
        k3 = diffEq(
            self.t[i] + self.h / 2.0,
            self.S[i] + self.h / 2.0 * k2,
            self.I[i] + self.h / 2.0 * k2,
            self.R[i] + self.h / 2.0 * k2
        )
        k4 = diffEq(
            self.t[i] + self.h,
            self.S[i] + self.h * k3,
            self.I[i] + self.h * k3,
            self.R[i] + self.h * k3
        )
        next_val = current_val + (self.h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        return next_val

    def ODESolve(self, method, diffEq):
        for i in range(self.N):
            pass
        return

    def __repr__(self):
        """
        A more informative __repr__ method. Useful when working with multiple instances
        in an interactive enviorment.
        """
        return "<SIRS Object: a=%.2f, b=%.2f, c=%.2f>" % (self.a, self.b, self.c)


def main():
    sys1 = SIRS(S0=300, I0=100, R0=0, a=4, b=1, c=0.5, N=1000, tN=10)

    print sys1
    return


if __name__ == '__main__':
    main()
