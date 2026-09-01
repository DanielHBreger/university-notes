import numpy as np
from scipy.optimize import brentq

def main():
    times = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8])
    counts = np.array([997, 520, 265, 127, 70, 35, 16, 7, 3])
    count0 = counts[0]

    def equation(lam):
        return np.sum((counts - count0 * np.exp(-lam * times)) * times * np.exp(-lam * times))

    lam_solution: float = brentq(equation, 1e-6, 10)  # type: ignore[assignment]
    print(f"lambda = {lam_solution:.6f}")

    residuals = counts - count0 * np.exp(-lam_solution * times)
    chi2 = np.sum(residuals**2)
    dof = len(counts) - 1  # one free parameter (lambda), N0 fixed
    print(f"chi^2 = {chi2:.4f}, dof = {dof}, chi^2/dof = {chi2/dof:.4f}")

    halflife = np.log(2) / lam_solution
    print(f"half-life T = {halflife:.4f} hours")

    # Skip i=0 (t=0 gives division by zero)
    Ni = counts[1:]
    ti = times[1:]
    N0 = count0

    denom = np.log(N0) - np.log(Ni)  # ln(N0) - ln(Ni)

    dT_dNi = np.log(2) * ti / (denom**2 * Ni)
    dT_dN0 = -np.log(2) * ti / (denom**2 * N0)
    dT_dti = np.log(2) / denom

    delta_Ni = 1
    delta_N0 = 1
    delta_ti = 5/60

    delta_T = np.sqrt(
        (dT_dNi * delta_Ni)**2 +
        (dT_dN0 * delta_N0)**2 +
        (dT_dti * delta_ti)**2
    )

    T_i = np.log(2) / denom * ti
    print("\nPer-measurement half-life estimates and errors:")
    for i in range(len(ti)):
        print(f"  t={ti[i]}: T = {T_i[i]:.4f} ± {delta_T[i]:.4f} hours")


    weights = 1 / delta_T**2
    T_avg = np.sum(weights * T_i) / np.sum(weights)
    delta_T_avg = 1 / np.sqrt(np.sum(weights))
    print(f"\nWeighted average: T = {T_avg:.4f} ± {delta_T_avg:.4f} hours")

if __name__ == "__main__":
    main()