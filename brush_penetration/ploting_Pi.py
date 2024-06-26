import os

import matplotlib.pyplot as plt
import numpy as np

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
path = os.path.abspath(__file__)

plt.figure(figsize=(7.5, 5))
for p in params:
    Nb, n, m = p
    X = np.loadtxt(f"{path[:path.rfind('/',2)]}/Nb{Nb}n{n}m{m}.txt")
    plt.plot(
        X[0, 1:] - 0.5 * (X[0, 1] - X[0, 0]),
        -np.diff(X[3]) / np.diff(X[0]),
        "-o",
        label=f"eta = {round((1 + n/m)**0.5,2)}",
    )

# plt.plot(X[0], 50*X[0]**(-1/3), '--', color='black')
# plt.plot(X[0], 20*X[0]**(-1/3), '--', color='black')
plt.plot(X[0], 3000 * X[0] ** (-5 / 2), "--", color="black")

plt.semilogx()
plt.semilogy()
plt.legend()
plt.xlabel("$d$", fontsize=16)
plt.ylabel("$\Pi(d)$", fontsize=16)
plt.title("$\sigma = 0.01$, slope -5/2")
plt.savefig("Pi_d.jpg")
