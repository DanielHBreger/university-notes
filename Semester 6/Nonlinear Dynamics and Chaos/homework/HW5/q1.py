import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def rossler_rhs(x, y, z, a, b, c):
    xdot = -y - z
    ydot = x + a * y
    zdot = b + z * (x - c)
    return xdot, ydot, zdot


def fixed_points(a, b, c):
    """Both fixed points of the Rossler system, found by solving F(x,y,z) = 0.

    xdot = 0  =>  y = -z
    ydot = 0  =>  x = -a*y = a*z
    Substituting into zdot = 0: b + z*(a*z - c) = 0, i.e. a*z^2 - c*z + b = 0.
    This quadratic's two roots give the inner (small z) and outer fixed point.
    """
    disc = c**2 - 4 * a * b
    z_inner = (c - np.sqrt(disc)) / (2 * a)
    z_outer = (c + np.sqrt(disc)) / (2 * a)

    def point(z):
        y = -z
        x = a * z
        assert np.allclose(rossler_rhs(x, y, z, a, b, c), 0)
        return np.array([x, y, z])

    return point(z_inner), point(z_outer)


def jacobian(a, x, y, z, c):
    return np.array([
        [0, -1, -1],
        [1, a, 0],
        [z, 0, x - c],
    ])

def part_a(a, b, c):
    inner, outer = fixed_points(a, b, c)
    print("Inner fixed point (near origin, small z):", inner)
    print("Outer fixed point:", outer)

    J = jacobian(a, inner[0], inner[1], inner[2], c)
    eigvals, eigvecs = np.linalg.eig(J)
    print("Jacobian at inner fixed point:\n", J)
    print("Eigenvalues:", eigvals)
    print("Eigenvectors (columns):\n", eigvecs)

def integrate_rossler(a, b, c, t_span=(0, 1000), n_points=100_000, ic=(0.1, 0.1, 0.1)):
    t_eval = np.linspace(t_span[0], t_span[1], n_points)
    sol = solve_ivp(
        lambda t, xyz: rossler_rhs(xyz[0], xyz[1], xyz[2], a, b, c),
        t_span, ic, method='RK45', rtol=1e-9, atol=1e-9, t_eval=t_eval,
    )
    return sol.t, sol.y[0], sol.y[1], sol.y[2]


def part_b(a, b):
    # (label, c) pairs targeting period-1, period-2, period-4, and chaos
    regimes = [
        ("period-1", 2.5),
        ("period-2", 3.5),
        ("period-4", 4.0),
        ("chaotic", 5.7),
    ]
    transient_frac = 0.5  # discard the first half of the integration

    n = len(regimes)
    fig = plt.figure(figsize=(10, 4 * n))

    for row, (label, c) in enumerate(regimes):
        t, x, y, z = integrate_rossler(a, b, c)
        keep = t > transient_frac * t[-1]
        x, y, z = x[keep], y[keep], z[keep]

        inner, outer = fixed_points(a, b, c)

        ax3d = fig.add_subplot(n, 2, 2 * row + 1, projection='3d')
        ax3d.plot(x, y, z, lw=0.5)
        ax3d.scatter(*inner, color='red', s=50, label='inner fixed point')
        ax3d.set_xlabel('x')
        ax3d.set_ylabel('y')
        ax3d.set_zlabel('z')
        ax3d.set_title(f"c = {c} ({label}) — 3D")
        ax3d.legend()

        ax2d = fig.add_subplot(n, 2, 2 * row + 2)
        ax2d.plot(x, y, lw=0.5)
        ax2d.scatter(inner[0], inner[1], color='red', s=50, label='inner fixed point')
        ax2d.set_xlabel('x')
        ax2d.set_ylabel('y')
        ax2d.set_title(f"c = {c} ({label}) — (x, y) projection")
        ax2d.legend()

    fig.suptitle("Rossler attractor across regimes")
    plt.tight_layout()
    plt.show()

def divergence(x, a, c):
    # div F = d(xdot)/dx + d(ydot)/dy + d(zdot)/dz = 0 + a + (x - c)
    return a + (x - c)


def part_c(a, b, c):
    t, x, y, z = integrate_rossler(a, b, c, t_span=(0, 2000), n_points=200_000)
    keep = t > 0.5 * t[-1]
    t, x, y, z = t[keep], x[keep], y[keep], z[keep]

    # (i) average contraction rate
    div = divergence(x, a, c)
    print(f"div F: mean = {div.mean():.4f}, min = {div.min():.4f}, max = {div.max():.4f}")

    # (ii) excursions from z = 0: linearizing zdot = b + z*(x-c) near z=0 gives
    # zdot ~ b + z*(x-c), which grows (pushes z away from 0) when x > c.
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    axes[0].plot(t, x, lw=0.6)
    axes[0].axhline(c, color='red', ls='--', label=f'x = c = {c}')
    axes[0].set_ylabel('x(t)')
    axes[0].legend()

    axes[1].plot(t, z, lw=0.6, color='green')
    axes[1].set_ylabel('z(t)')

    axes[2].plot(t, div, lw=0.6, color='purple')
    axes[2].axhline(0, color='black', ls='--', lw=0.8)
    axes[2].set_ylabel(r'$\nabla \cdot F$')
    axes[2].set_xlabel('t')

    fig.suptitle(f"Contraction and z-excursions at c = {c}")
    plt.tight_layout()
    plt.show()

    return div.mean()


def part_d(a, b, c, mean_div, t_span=(0, 200), n_points=200_000, ic=(0.1, 0.1, 0.1),
           perturbation=1e-8, fit_window=(9, 106)):
    # discard transient on the reference trajectory first, then perturb from
    # a point already on the attractor
    t0, x0, y0, z0 = integrate_rossler(a, b, c, t_span=(0, 200), n_points=200_000, ic=ic)
    on_attractor = np.array([x0[-1], y0[-1], z0[-1]])

    perturbed_ic = on_attractor + perturbation * np.array([1.0, 0.0, 0.0])

    t1, x1, y1, z1 = integrate_rossler(a, b, c, t_span=t_span, n_points=n_points, ic=on_attractor)
    t2, x2, y2, z2 = integrate_rossler(a, b, c, t_span=t_span, n_points=n_points, ic=perturbed_ic)

    delta = np.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
    log_delta = np.log(delta)

    # fit the region of clean exponential growth: skip the initial transient
    # (the perturbation direction must first align with the leading Lyapunov
    # direction) and stop before delta saturates at the attractor size
    fit_mask = (t1 >= fit_window[0]) & (t1 <= fit_window[1])
    slope, intercept = np.polyfit(t1[fit_mask], log_delta[fit_mask], 1)
    lambda1 = slope

    # |delta(t)| ~ delta_0 * exp(lambda1 * t), so reaching delta_max ~ 1 takes
    # T_pred ~ (1/lambda1) * ln(delta_max/delta_0)
    e_folding_time = 1.0 / lambda1
    prediction_horizon = np.log(1.0 / perturbation) / lambda1
    print(f"lambda1 ~ {lambda1:.4f}")
    print(f"e-folding time 1/lambda1 ~ {e_folding_time:.4f}")
    print(f"Prediction horizon ~ {prediction_horizon:.4f}")

    lambda2 = 0.0
    lambda3 = mean_div - lambda1 - lambda2
    print(f"lambda2 = {lambda2}, lambda3 = <div F> - lambda1 - lambda2 = {lambda3:.4f}")

    # Kaplan-Yorke dimension for lambda1 > lambda2 > 0 > lambda3:
    d_ky = 2 + (lambda1 + lambda2) / abs(lambda3)
    print(f"Kaplan-Yorke dimension D_KY ~ {d_ky:.4f}")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(t1, log_delta, lw=0.7, label=r'$\ln|\delta(t)|$')
    fit_t = t1[fit_mask]
    ax.plot(fit_t, slope * fit_t + intercept, 'r--',
            label=fr'fit: $\lambda_1$ = {lambda1:.3f}')
    ax.set_xlabel('t')
    ax.set_ylabel(r'$\ln|\delta(t)|$')
    ax.set_title(f'Divergence of nearby trajectories at c = {c}')
    ax.legend()
    plt.tight_layout()
    plt.show()

    return lambda1, lambda2, lambda3, d_ky


def find_x_maxima(t, x, z):
    """Successive local maxima of x(t): points where xdot = 0, xddot < 0.

    Found as interior samples that are larger than both neighbors, i.e. a
    sign change of the discrete derivative from + to -. Returns the maxima
    values x_n and the corresponding z at those times.
    """
    is_peak = (x[1:-1] > x[:-2]) & (x[1:-1] > x[2:])
    peak_idx = np.where(is_peak)[0] + 1
    return x[peak_idx], z[peak_idx]


def part_e(a, b, c):
    t, x, y, z = integrate_rossler(a, b, c, t_span=(0, 2000), n_points=400_000)
    keep = t > 0.5 * t[-1]
    t, x, y, z = t[keep], x[keep], y[keep], z[keep]

    x_max, z_at_max = find_x_maxima(t, x, z)
    print(f"Found {len(x_max)} maxima of x(t)")

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(x_max, z_at_max, s=8)
    ax.set_xlabel(r'$x_n$')
    ax.set_ylabel(r'$z_n$')
    ax.set_title(f'Poincare section: successive maxima of x(t), c = {c}')
    plt.tight_layout()
    plt.show()

    return x_max


def part_f(x_max, c):
    x_n = x_max[:-1]
    x_next = x_max[1:]

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(x_n, x_next, s=8)
    ax.set_xlabel(r'$x_n$')
    ax.set_ylabel(r'$x_{n+1}$')
    ax.set_title(f'Lorenz (next-maximum) map, c = {c}')
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.show()


def main():
    a, b, c = 0.2, 0.2, 5.7
    part_a(a, b, c)
    part_b(a, b)
    mean_div = part_c(a, b, c)
    part_d(a, b, c, mean_div)
    x_max = part_e(a, b, c)
    part_f(x_max, c)


if __name__ == "__main__":
    main()