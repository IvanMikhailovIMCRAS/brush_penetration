import os

import matplotlib.pyplot as plt
import numpy as np
from sfbox_api import Composition, Frame, Lat, Mol, Mon, Sys


def brushes(Nb: int, n: int, m: int, n_layers, sigma: float) -> Frame:
    if n <= 0:
        comp1 = f"(X1)1((A1){m}){Nb//m-1}(A1){m}"
        comp2 = f"(X2)1((A2){m}){Nb//m-1}(A2){m}"
    else:
        comp1 = f"(X1)1((A1){m}[(A1){n}]){Nb//m-1}(A1){m}"
        comp2 = f"(X2)1((A2){m}[(A2){n}]){Nb//m-1}(A2){m}"
    theta = (Nb + 1 + (Nb // m - 1) * n) * sigma
    lat = Lat(
        **{
            "n_layers": n_layers,
            "geometry": "flat",
            "lowerbound": "surface",
            "upperbound": "surface",
        }
    )
    mons = [
        Mon(**{"name": "X1", "freedom": "pinned", "pinned_range": "1"}),
        Mon(**{"name": "A1", "freedom": "free"}),
        Mon(**{"name": "X2", "freedom": "pinned", "pinned_range": str(n_layers)}),
        Mon(**{"name": "A2", "freedom": "free"}),
        Mon(**{"name": "W", "freedom": "free"}),
    ]
    mols = [
        Mol(**{"name": "Water", "composition": "(W)1", "freedom": "solvent"}),
        Mol(
            **{
                "name": "pol1",
                "composition": comp1,
                "freedom": "restricted",
                "theta": theta,
            }
        ),
        Mol(
            **{
                "name": "pol2",
                "composition": comp2,
                "freedom": "restricted",
                "theta": theta,
            }
        ),
    ]
    sys = Sys()
    chi = 0.0
    chi_list = {"X1 W": chi, "A1 W": chi, "X2 W": chi, "A2 W": chi}

    frame = Frame(lat, sys, mols, mons, chi_list=chi_list)
    frame.run()

    return frame


if __name__ == "__main__":
    params = [[160, 56, 10], [136, 54, 8], [112, 8, 1]]

    sigma = 0.01
    for p in params:
        Nb, n, m = p
        data = []
        theta = (Nb + 1 + (Nb // m - 1) * n) * sigma

        for n_layers in range(42, 116, 2):
            frame = brushes(Nb, n, m, n_layers, sigma)
            omega_L = np.sum(frame.profile["pol1"] * frame.profile["pol2"])
            if len(frame.profile["pol1"]) % 2 == 0:
                first_layer = len(frame.profile["pol1"]) // 2 + 1
                H = 2 * (
                    np.sum(
                        frame.profile["pol1"][first_layer:]
                        * (frame.profile["layer"][first_layer:] + 0.5 - first_layer)
                    )
                    / np.sum(frame.profile["pol1"][first_layer:])
                )
                data.append(
                    [n_layers, H, omega_L, frame.stats["sys : name : free energy"]]
                )
            print(n_layers)

        data = np.array(data).T
        path = os.path.abspath(__file__)
        np.savetxt(f"{path[:path.rfind('/',2)]}/Nb{Nb}n{n}m{m}.txt", data)
