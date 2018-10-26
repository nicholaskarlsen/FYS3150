import numpy as np
from numba import jitclass          # import the decorator
from numba import float32    # import the types

spec = [
    ('a', float32),
    ('b', float32),
]


@jitclass(spec)
class test(object):
    def __init__(self, a, b):
        self.a = a
        self.b = b

    def f(self, val):
        return self.a * val

    def g(self, func):
        val = 3
        return self.b + func(val)

    def wrap(self):
        return self.g(self.f)


obj = test(1, 2)

c = obj.wrap()
print(c)  # prints '5' when @jitclass(spec) is commented out
