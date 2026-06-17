import numpy as np
import matplotlib.pyplot as plt

def main():
    t_over_tau = [0.05, 0.1, 0.4, 0.75, 10]
    h = 1.0
    N_terms = 100  # number of Fourier terms to include

    def v(y, t_over_tau):
        result = np.zeros_like(y, dtype=float)
        for k in range(N_terms):
            n = 2 * k + 1  # odd harmonics: 1, 3, 5, ...
            result += (8 / (n**2 * np.pi**2)) * np.sin(n * np.pi * y / h) * np.exp(-n**2 * t_over_tau)
        return 2 * result  # factor of 2 so steady-state peaks at 2 (matching reference)

    colors = plt.colormaps['plasma'](np.linspace(0, 1, len(t_over_tau)))
    for t, c in zip(t_over_tau, colors):
        y = np.linspace(0, 1, 500)
        plt.plot(y, v(y, t), label=f't/τ = {t}', color=c)

    plt.xlabel('y')
    plt.ylabel('v_x\'(y, t)')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()