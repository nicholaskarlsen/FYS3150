
class SIRS(object):
    def __init__(self, S0, I0, R0, dSdt, dIdt, dRdt, a=0, b=0, c=0, d=0, d_I=0, f=0):
        self.S0 = S0
        self.I0 = I0
        self.R0 = R0

        if isinstance(a, (float, int)):
            self.a = lambda t: a
        elif callable(a):
            self.a = a

        self.b = b
        self.c = c
        self.d = d
        self.d_I = d_I

        if isinstance(f, (float, int)):
            self.f = lambda t: f
        elif callable(f):
            self.f = f

    def __call__(self, u, t):
        return self.dSdt, self.dIdt, self.dRdt
