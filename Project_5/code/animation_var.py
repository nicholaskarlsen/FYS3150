# Python v2.7
# By Nicholas Karlsen

# Creates an animation, showing behaviour of function when varyiong one of the parameters

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from SIRS_ODE import SIRS

omegas_A = np.linspace(4, 1, 50)
omegas_B = np.linspace(1, 1e-5, 100)
omegas = np.append(omegas_A, omegas_B)
NO_FRAMES = len(omegas)
inst = SIRS(S0=300, I0=100, R0=0, a=4, b=1, c=0.5, N=10000, tN=204, Amplitude=4, omega=omegas[0])
inst.solve(inst.sirs_svar)
t, S, I, R = inst.get()

fig, ax = plt.subplots()

ax.set_ylim(0, 400)
ax.set_xlim(t[0], t[-1])

lineS, = ax.plot(t, S, label="S")
lineI, = ax.plot(t, I, label="I")
lineR, = ax.plot(t, R, label="R")
text = ax.text(0.25, (S[0] + I[0] + R[0]) * .95, "")

ax.legend(loc="upper right")


def init():  # only required for blitting to give a clean slate.
    lineS.set_ydata([np.nan] * len(t))
    lineI.set_ydata([np.nan] * len(t))
    lineR.set_ydata([np.nan] * len(t))

    return (lineS, lineI, lineR) + (text,)


def animate(i):
    inst = SIRS(S0=300, I0=100, R0=0, a=4, b=1, c=0.5, N=10000, tN=204, Amplitude=2, omega=omegas[i])
    inst.solve(inst.sirs_svar)
    t, S, I, R = inst.get()
    lineS.set_ydata(S)  # update the data.
    lineI.set_ydata(I)  # update the data.
    lineR.set_ydata(R)  # update the data.
    text.set_text("$\\omega$ = %.2f" % omegas[i])
    return (lineS, lineI, lineR) + (text,)


ani = animation.FuncAnimation(
    fig, animate, init_func=init, interval=100, blit=False, frames=np.arange(0, NO_FRAMES))

# To save the animation, use e.g.
#
# ani.save("movie.mp4")
#
# or
#
# from matplotlib.animation import FFMpegWriter
# writer = FFMpegWriter(fps=15, metadata=dict(artist='Me'), bitrate=1800)
# ani.save("movie.mp4", writer=writer)
ani.save("../figs/SIRS_svar.gif", dpi=80, writer='imagemagick')
plt.close()
