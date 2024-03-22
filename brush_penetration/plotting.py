import matplotlib.pyplot as plt
import numpy as np

X0 = np.loadtxt("brush_penetration/Nb1000n0m100.txt")
X = np.loadtxt("brush_penetration/N1000_sigma0_01_linear.txt")
D0 = X0[0]
L = X0[1]

D = X[0]
L1 = X[1]
L2 = X[2]
phi_mean = X[3]
F_int = X[4]
F = X[5]

# plt.plot(D, L1, label='Leermakers')
# plt.plot(D, L2, label='Amoskov')
# plt.plot(D0, L, label='Our old version')
# plt.plot(D, 44*D**(-1/3), '--', color='black')
# plt.plot(D, 66*D**(-1/3), '--', color='black')
# plt.xlabel('$D$', fontsize=14)
# plt.ylabel('$L$', fontsize=14)
# plt.title("$N = 1000,~~\sigma = 0.01, ~~ \chi=0.0$")
# plt.legend()
# plt.loglog()
# plt.savefig('L_mesure.jpg')
# plt.close()

plt.plot(D, F + 20.0, label="total free energy")
plt.plot(D, F_int + 20.0, label="osmotic contribution")
plt.plot(D, 300 * D ** (-1), "--", color="black")
plt.title("$N = 1000,~~\sigma = 0.01, ~~ \chi=0.0$")
plt.legend()
plt.loglog()
# plt.savefig('L_mesure.jpg')
plt.show()
plt.close()
