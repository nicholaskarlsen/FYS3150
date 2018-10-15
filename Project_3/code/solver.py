#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Written by Nicholas Karlsen
# Python version 2.7.15
"""
A general purpose ODE-solver that i have written previously, and extended to include
the velocity-verlet algorithm for this project.
"""
from __future__ import division  # Nobody expects the integer division
import numpy as np


class diffeqsolver:
    def __init__(self, x_init, v_init, t_init):
        self.x_init = x_init    # Initial position
        self.v_init = v_init    # Initial velocity
        self.t_init = t_init    # Initial time

    def eulercromer(self, current_pos, current_vel, current_t, dt, diffEq):
        """ Solves one timestep of a differential equation using the Euler-cromer
        method. Mostly used as a test function, since it is easier
        to implement correctly than RK4.
        """
        next_vel = current_vel + dt * self.diffEq(current_pos,
                                                  current_vel,
                                                  current_t)
        next_pos = current_pos + dt * next_vel
        return next_pos, next_vel

    def rk4(self, current_pos, current_vel, current_t, dt, diffEq):
        """ 
        Solves one timestep of a differential equation using the
        Runge-Kutta4 method. Slightly modified version of
        what is written in the textbook.
        """
        a1 = diffEq(current_pos, current_vel, current_t)
        v1 = current_vel
        pos_half1 = current_pos + v1 * dt / 2.0
        vel_half1 = current_vel + a1 * dt / 2.0
        a2 = diffEq(pos_half1, vel_half1, current_t + dt / 2.0)
        v2 = vel_half1
        pos_half2 = current_pos + v2 * dt / 2.0
        vel_half2 = current_vel + a2 * dt / 2.0
        a3 = diffEq(pos_half2, vel_half2, current_t + dt / 2.0)
        v3 = vel_half2
        next_pos = current_pos + v3 * dt
        next_vel = current_vel + a3 * dt
        a4 = diffEq(next_pos, next_vel, current_t + dt)
        v4 = next_vel
        acc_mid = 1.0 / 6.0 * (a1 + 2 * a2 + 2 * a3 + a4)
        vel_mid = 1.0 / 6.0 * (v1 + 2 * v2 + 2 * v3 + v4)
        next_pos = current_pos + vel_mid * dt
        next_vel = current_vel + acc_mid * dt

        return next_pos, next_vel

    def velocityverlet(self, current_pos, current_vel, current_t, dt, diffEq):
        return

    def solver(self, tn, dt, method, diffEq):
        """ Solves a differential equation using a specified method
        For diffEq, either one of the built in functions can be called,
        or a self made function not contained in this class as long
        as it is structured similarly; i.e diff(pos, vel, t).
        Same goes for method.
        """
        time = np.linspace(0, tn, int(tn) / dt)
        position = np.zeros(len(time))
        velocity = np.zeros(len(time))
        position[0] = self.x_init
        velocity[0] = self.v_init
        for i in xrange(len(time) - 1):
            position[i + 1], velocity[i + 1] = method(position[i], velocity[i],
                                                      time[i], dt, diffEq)

        return position, velocity, time
