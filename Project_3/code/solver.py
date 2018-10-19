#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Written by Nicholas Karlsen
# Python version 2.7.15

"""
An ODE solver which solves for single and N-body systems for systems which is governed
by which are functions of their position.
"""

from __future__ import division  # Nobody expects the integer division
import numpy as np
from scipy import constants as const


class diffeqsolver:
    def __init__(self, x_init, v_init):
        self.x_init = x_init    # Initial position
        self.v_init = v_init    # Initial velocity
        self.G = const.G

    def eulercromer(self, current_pos, current_vel, current_t, dt, diffEq):
        next_vel = current_vel + dt * diffEq(current_pos)
        next_pos = current_pos + dt * next_vel
        return next_pos, next_vel

    def velocityverlet(self, current_pos, current_vel, current_t, dt, diffEq):
        a1 = diffEq(current_pos, 0, 0)   # a_i
        next_pos = current_pos + current_vel * dt + 0.5 * a1 * dt**2
        a2 = diffEq(next_pos, 0, 0)
        next_vel = current_vel + 0.5 * dt * (a1 + a2)
        return next_pos, next_vel

    def solve(self, tn, dt, method, diffEq):
        time = np.linspace(0, tn, int(tn) / dt)
        position = np.zeros([len(time), 2])
        velocity = np.zeros([len(time), 2])
        position[0] = self.x_init
        velocity[0] = self.v_init
        for i in xrange(len(time) - 1):
            position[i + 1], velocity[i + 1] = method(position[i],
                                                      velocity[i],
                                                      time[i],
                                                      dt, diffEq)

        return position, velocity, time
