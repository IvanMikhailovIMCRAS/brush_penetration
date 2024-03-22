import numpy as np
from sfbox_api import Composition, Frame, Lat, Mol, Mon, Sys

from brush_penetration import brushes

data = []

if __name__ == "__main__":
    for d in range(40, 400):
        frame = brushes(Nb=1000, n=0, m=10, n_layers=d, sigma=0.01)
        # print(len(frame.profile['pol1']))
        phi1 = frame.profile["pol1"][1:-1]
        phi2 = frame.profile["pol2"][1:-1]
        z = frame.profile["layer"][1:-1]
        z_mean = 0.5 * (z[0] + z[-1])
        if len(z) % 2 == 0:
            z0 = len(z) // 2 + 1
        else:
            z0 = len(z) // 2
        L1 = 2.0 * np.sum(phi1[z0:] * (z[z0:] + 0.5 - z_mean)) / np.sum(phi1[z0:])

        L2 = np.sum(phi1[:] * phi2[:] * np.square(z[:] + 0.5 - z_mean)) / np.sum(
            phi1[:] * phi2[:]
        )
        L2 = 2.0 * L2**0.5

        phi_mean = phi1[z0] + phi2[z0]

        F_int = np.sum((1 - phi1[:] - phi2[:]) * np.log(1 - phi1[:] - phi2[:]))

        L_Omega = np.sum(phi1[:] * phi2[:])

        F = frame.stats["sys : name : free energy"]

        print(d, L1, L2, phi_mean, F_int, F)

        data.append([d, L1, L2, phi_mean, F_int, F])
    np.savetxt("N1000_sigma0_01_linear.txt", np.array(data).T)
