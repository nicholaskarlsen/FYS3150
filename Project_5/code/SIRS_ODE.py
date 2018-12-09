from __future__ import division  # Nobody expects the integer division
import numpy as np
import matplotlib.pyplot as plt


class SIRS:
    def __init__(self, S0, I0, R0, N, tN, a=0, b=0, c=0, d=0, d_I=0, e=0, f=0,
                 Amplitude=0, omega=0):
        """
        Solves the SIRS model numerically as a differential equation using Runge-Kutta 4,
        or Euler forward (for checking correct implementation).
        (Originally planned to contian both ODE & MCMC Solvers, but turned out
        to be impractical, and slow, hence two wildy different aproaches)

        Always initializing the a, b, ... variables makes it so that the class can be
        called, and solved in a very similar fashion for all the different problems,
        but detracts from the expandability.

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
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.d_I = d_I
        self.e = e
        self.f = f

        self.S[0] = S0
        self.I[0] = I0
        self.R[0] = R0

        return

    def sirs_basic(self, t, S, I, R):
        """
        Compute differentials for the SIRS system at a time t. Kept as a nested function
        to access the local namespace for a, b, c constants.
        ("Wasted" FLOPS by multiplying by variables which =0, but requires less methods.
        If runtime becomes issue, add new method with only the speciffic traits)

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
        N = S + I + R
        dSdt = self.c * R - self.a * S * I / N - self.d * S + self.e * N - self.f
        dIdt = self.a * S * I / N - self.b * I - self.d * I - self.d_I * I
        dRdt = self.b * I - self.c * R - self.d * R + self.f

        return np.array([dSdt, dIdt, dRdt])

    def sirs_svar(self, t, S, I, R):
        N = S + I + R
        avar = self.Amplitude * np.cos(self.omega * t) + self.a

        dSdt = self.c * R - avar * S * I / N
        dIdt = self.a * S * I / N - self.b * I
        dRdt = self.b * I - self.c * R

        return np.array([dSdt, dIdt, dRdt])

    def euler_fw(self, i, diffEq):
        """
        Computes time-step using the Euler-Forward method
        (used to check implementation of RK4. if both yield same result -> things
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

    def solve(self, diffEq):
        """
        Parameters
        ----------
        diffEq : Method which outputs derivative of system with input (t, S, I, R) and
                and outputs satisfying requirements for rk4 method.

        Returns
        -------
        N/A

        """
        for i in range(self.N - 1):
            self.S[i + 1], self.I[i + 1], self.R[i + 1] = self.rk4(i, diffEq)

        return

    def get(self):
        """
        Method for extracting the arrays
        Parameters
        ----------
        N/a

        Returns
        -------
        list of numpy arrays
        """

        return self.t, self.S, self.I, self.R

    def plot(self):  # Added SIR at end to avoid confusion with matplot
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
    sys1 = SIRS(S0=300, I0=100, R0=0, a=4, b=1, c=0.5, e=1, d=1, d_I=0, N=100, tN=20)
    sys1.solve(sys1.sirs_basic)
    sys1.plot()
    return


if __name__ == '__main__':
    main()
