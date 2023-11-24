from sfbox_api import Frame, Lat, Mol, Mon, Sys
import matplotlib.pyplot as plt
import numpy as np

def comb_brush(Nb: int, n: int, m: int, n_layers, sigma: float) -> Frame:
    if n <= 0:
        if m > 1:
            comp = f"(X)1((A){m}){Nb//m-1}(A){m-1}(G)1"
        else:
            comp = f"(X)1((A){m}){Nb//m-1}(G)1"
    else:
        if m > 1:
            comp = f"(X)1((A){m}[(A){n}]){Nb//m-1}(A){m-1}(G)1"
        else:
            comp = f"(X)1((A){m}[(A){n}]){Nb//m-1}(G)1"
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
        Mon(**{"name": "X", "freedom": "pinned", "pinned_range": "1"}),
        Mon(**{"name": "A", "freedom": "free"}),
        Mon(**{"name": "G", "freedom": "free"}),
        Mon(**{"name": "W", "freedom": "free"}),
    ]
    mols = [
        Mol(**{"name": "Water", "composition": "(W)1", "freedom": "solvent"}),
        Mol(
            **{
                "name": "pol",
                "composition": comp,
                "freedom": "restricted",
                "theta": theta,
            }
        ),
    ]
    sys = Sys()
    chi = 0.0
    chi_list = {"X W": chi, "A W": chi}

    frame = Frame(lat, sys, mols, mons, chi_list=chi_list)
    frame.run()

    return frame

def eta(n: int, m: int) -> float:
    return (1.0 + n / m)**0.5

def H(n, m, sigma, N):
    return (4.0 * sigma / eta(n,m)**2 / np.pi**2)**(1/3) * N

def U0_predict(n, m, sigma, Nb):
    N = 1 + Nb + (Nb // m - 1) * n
    return 3/8 * np.pi**2 / N**2 * eta(n,m)**2 * H(n, m, sigma, N)**2
    

if __name__ == "__main__":
    params = [
    [1000, 0, 100],
    [256, 24, 8],
    [252, 44, 14],
    [240, 40, 12],
    [230, 35, 10],
    [190, 45, 10],
    [160, 56, 10],
    [136, 54, 8],
    [112, 8, 1],
            ]
    sigma = 0.01
    n_layers = 200
    plt.figure(figsize=(7.5, 5))
    for p in params:
        Nb, n, m = p
        frame = comb_brush(Nb, n, m, n_layers, sigma)
        z = frame.profile['layer'][1:-2]
        u = - np.log(frame.profile['Water'][1:-2])
        g = frame.profile['G'][1:-2]
        g = g / np.sum(g)
        u0 = U0_predict(n, m, sigma, Nb)
        N = 1 + Nb + (Nb // m - 1) * n
        h = H(n, m, sigma, N)
        # phi = frame.profile['pol'][1:-2]
        plt.plot(z/h, g*h, '-o', label=f'$\eta={round((1+n/m)**0.5, 2)}$')
    # plt.plot(z, phi, label='phi')
    # plt.plot(z/z[-1], 1 - z/z[-1], '--', color='black')
    plt.xlabel("$z/H$", fontsize=16)
    plt.ylabel("$g(z)H$", fontsize=16)
    plt.legend()
    plt.xlim(0, 1.1)
    plt.ylim(0, 1.8)
    plt.title(f"$\sigma = {sigma}$", fontsize=16)
    plt.savefig("g_z.jpg")