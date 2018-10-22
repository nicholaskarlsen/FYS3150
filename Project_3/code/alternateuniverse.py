# Simulating an alternate universe where the inverse square law creeps up to an inverse cube

from solarsystem import *


class alternateuniverse(solarsystem):
    def __init__(self, initPos, initVel, mass, N, tn, names, r_exp):
        solarsystem.__init__(self, initPos, initVel, mass, N, tn, names)
        self.r_exp = r_exp

        return

    def vargravity(self, planetIndex, timeIndex):
        accel = np.zeros(self.dim)
        for j in xrange(self.numPlanets):
            if j != planetIndex:  # Not interested in self-attraction
                relPos = self.pos[j, timeIndex] - self.pos[planetIndex, timeIndex]
                accel += relPos * self.G * self.mass[j] / (np.linalg.norm(relPos) ** (1 + self.r_exp))

        return accel


if __name__ == '__main__':
    m = np.array([3.00348959632E-6, 1])
    x0 = np.array([[1, 0], [0, 0]])
    v0 = np.array([[0, 2 * np.pi], [0, 0]])
    #planets = {"Earth": 399, "Sun": 10}
    #x01, v01, m1 = hori.fetch_data(jpl_id=planets)
    N = 1e5
    tn = 1
    names = ["Earth", "Sun"]

    for i in [2, 2.25, 2.5, 2.75, 3]:
        esys = alternateuniverse(r_exp=i, initPos=x0, initVel=v0, mass=m, N=N, tn=tn, names=names)
        esys.velocityverlet(diffeq=esys.vargravity)
        pos, vel = esys.get()
        epos = pos[0] - pos[1]   # Position of earth relative to sun, correcting for any movement non-fixed sun may have.
        plt.plot(epos[:, 0], epos[:, 1], label="r=%.2f"%i)
    plt.legend()
    plt.show()
