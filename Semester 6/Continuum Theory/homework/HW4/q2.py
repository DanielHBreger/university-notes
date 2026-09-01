import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from pathlib import Path

eta_max = 6.0
s0 = 2 / np.sqrt(3)  # f''(0) from the first integral of g'' = g^2 - 1, g=f'

def rhs(eta, y):
    f, fp, fpp = y
    return [fp, fpp, fp**2 - 1]

eta_eval = np.linspace(0, eta_max, 600)
sol = solve_ivp(rhs, (0, eta_max), [0.0, 0.0, s0], t_eval=eta_eval, rtol=1e-10, atol=1e-12)

plt.figure(figsize=(6, 4))
plt.plot(sol.t, sol.y[1])
plt.xlabel(r"$\eta=y/\delta(x)$")
plt.ylabel(r"$u/U=f'(\eta)$")
# plt.title("Boundary-layer profile for sink flow")
plt.grid(True)
plt.tight_layout()
plt.show()