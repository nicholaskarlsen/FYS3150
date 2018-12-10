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
inst = SIRS(S0=300, I0=100, R0=0, a=4, b=1, c=0.5, N=10000, tN=204, Amplitude=1, omega=omegas[0])
inst.solve(inst.sirs_svar)
t, S, I, R = inst.get()

fig, (ax1, ax2) = plt.subplots(1,2)

ax1.set_ylim(0, 400)
ax1.set_xlim(t[0], t[-1])

lineS1, = ax1.plot(t, S, label="S")
lineI1, = ax1.plot(t, I, label="I")
lineR1, = ax1.plot(t, R, label="R")
text1 = ax1.text(0.25, (S[0] + I[0] + R[0]) * .95, "")

ax2.set_ylim(0, 400)
ax2.set_xlim(t[0], t[-1])

lineS2, = ax2.plot(t, S, label="S")
lineI2, = ax2.plot(t, I, label="I")
lineR2, = ax2.plot(t, R, label="R")
text2 = ax2.text(0.25, (S[0] + I[0] + R[0]) * .95, "")

ax1.legend(loc="upper right")
ax2.legend(loc="upper right")


def init():  # only required for blitting to give a clean slate.
    lineS1.set_ydata([np.nan] * len(t))
    lineI1.set_ydata([np.nan] * len(t))
    lineR1.set_ydata([np.nan] * len(t))

    lineS2.set_ydata([np.nan] * len(t))
    lineI2.set_ydata([np.nan] * len(t))
    lineR2.set_ydata([np.nan] * len(t))


    return (lineS1, lineI1, lineR1, lineS2, lineI2, lineR2) + (text1, text2)


def animate(i):
    inst = SIRS(S0=300, I0=100, R0=0, a=4, b=1, c=0.5, N=10000, tN=204, Amplitude=1, omega=omegas[i])
    inst.solve(inst.sirs_svar)
    t, S, I, R = inst.get()
    lineS1.set_ydata(S)  # update the data.
    lineI1.set_ydata(I)  # update the data.
    lineR1.set_ydata(R)  # update the data.
    text1.set_text("A=1, $\\omega$ = %.2f" % omegas[i])
    inst = SIRS(S0=300, I0=100, R0=0, a=4, b=1, c=0.5, N=10000, tN=204, Amplitude=4, omega=omegas[i])
    inst.solve(inst.sirs_svar)
    t, S, I, R = inst.get()
    lineS2.set_ydata(S)  # update the data.
    lineI2.set_ydata(I)  # update the data.
    lineR2.set_ydata(R)  # update the data.
    text2.set_text("A=4, $\\omega$ = %.2f" % omegas[i])
    return (lineS1, lineI1, lineR1) + (text1,)


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
ani.save("../figs/SIRS_seasonal_var.gif", dpi=80, writer='imagemagick')
plt.close()
