import numpy as np
import matplotlib.pyplot as plt

def main():
    ts = np.arange(4,15,1)
    x1 = [18.81,20.01,21.35,18.73,20.57,22.28,21.26,21.36,23.63,20.15,22.3]
    x2 = [17.66,19.05,20.25,21.15,21.84,22.15,22.4,22.27,21.9,21.35,20.43]
    x3 = [19.53,20.19,20.55,19.72,19.85,20.25,22.94,22.38,21.55,19.65,23.84]
    x4 = [19.9,20.09,21.21,20.31,20.43,20.56,20.7,20.8,20.95,25.31,21.19]
    sigma_x = 1
    v_estimator = lambda x,t: (np.average(t*x)-np.average(t)*np.average(x))/(np.average(t**2) - np.average(t)**2)
    v1 = v_estimator(x1, ts)
    v2 = v_estimator(x2, ts)
    v3 = v_estimator(x3, ts)
    v4 = v_estimator(x4, ts)
    print(f"v1 = {v1:.4f}, v2 = {v2:.4f}, v3 = {v3:.4f}, v4 = {v4:.4f}")

    v_estimator_error = lambda sigma, t, N: sigma**2/(N*(np.average(t**2) - np.average(t)**2))
    N = len(ts)
    error1 = np.sqrt(v_estimator_error(sigma_x, ts, N))
    error2 = np.sqrt(v_estimator_error(sigma_x, ts, N))
    error3 = np.sqrt(v_estimator_error(sigma_x, ts, N))
    error4 = np.sqrt(v_estimator_error(sigma_x, ts, N))
    print(f"error1 = {error1:.4f}, error2 = {error2:.4f}, error3 = {error3:.4f}, error4 = {error4:.4f}")
    chi_square = lambda sigma, x, t, v, x0: np.sum((x - (x0 + v*t))**2) / sigma**2
    x0_estimator = lambda x, t, v: np.average(x) - v*np.average(t)
    x0_1 = x0_estimator(x1, ts, v1)
    x0_2 = x0_estimator(x2, ts, v2)
    x0_3 = x0_estimator(x3, ts, v3)
    x0_4 = x0_estimator(x4, ts, v4)
    print(f"x0_1 = {x0_1:.4f}, x0_2 = {x0_2:.4f}, x0_3 = {x0_3:.4f}, x0_4 = {x0_4:.4f}")
    chi2_1 = chi_square(sigma_x, x1, ts, v1, x0_1)
    chi2_2 = chi_square(sigma_x, x2, ts, v2, x0_2)
    chi2_3 = chi_square(sigma_x, x3, ts, v3, x0_3)
    chi2_4 = chi_square(sigma_x, x4, ts, v4, x0_4)
    dof = N - 2
    print(f"chi2_1 = {chi2_1:.4f}, chi2_2 = {chi2_2:.4f}, chi2_3 = {chi2_3:.4f}, chi2_4 = {chi2_4:.4f}")
    print(f"dof = {dof}")

    plt.figure(figsize=(10, 6))
    plt.plot(ts, x1, 'o-', label='x1', color='orange')
    plt.plot(ts, x2, 's-', label='x2', color='green')
    plt.plot(ts, x3, '^-', label='x3', color='blue')
    plt.plot(ts, x4, 'd-', label='x4', color='red')
    # linear fits
    plt.plot(ts, x0_1 + v1*ts, 'r--', label=f'Fit x1: v={v1:.2f}', color='orange')
    plt.plot(ts, x0_2 + v2*ts, 'g--', label=f'Fit x2: v={v2:.2f}', color='green')
    plt.plot(ts, x0_3 + v3*ts, 'b--', label=f'Fit x3: v={v3:.2f}', color='blue')
    plt.plot(ts, x0_4 + v4*ts, 'm--', label=f'Fit x4: v={v4:.2f}', color='red')
    plt.xlabel('Time')
    plt.ylabel('Position')
    plt.title('Position vs Time')
    plt.legend()
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    main()