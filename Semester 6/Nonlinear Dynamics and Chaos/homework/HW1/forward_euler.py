import sympy as sp
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import norm

ONE_ORDER_OF_MAGNITUDE = 10


def plot_trajectories(trajectories, title):
    for initial_condition, trajectory in trajectories.items():
        plt.plot(trajectory, label=f"x0={initial_condition}")
    plt.title(title)
    plt.xlabel('Time steps')
    plt.ylabel('x')
    plt.legend()
    plt.grid()
    max_length = min(len(trajectory) for trajectory in trajectories.values())
    plt.xlim(0, max_length)
    plt.show()


def relaxation_time(equation, fixed_point):
    x = sp.symbols('x')
    linearized = sp.diff(equation.rhs, x).subs(x, fixed_point)
    if linearized == 0:
        return float('inf')
    return sp.sympify(1) / sp.Abs(linearized)


def integration_time_prediction(tau):
    return tau * ONE_ORDER_OF_MAGNITUDE


def integrate_forward_euler(f, x0, dt, T):
    dt = float(dt)
    steps = int(T / dt)
    trajectory = [x0]
    for _ in range(steps):
        try:
            x_next = trajectory[-1] + dt * f(trajectory[-1])
        except OverflowError:
            break
        if not np.isfinite(x_next):
            break
        trajectory.append(x_next)
    return trajectory


def return_closest_dt(dts, x):
    closest_point = min(dts, key=lambda point: abs(point - x))
    return dts[closest_point]


def part_a(fixed_points, xdot, r):
    x = sp.symbols('x')
    deriv = sp.diff(xdot.rhs, x)
    dts = {}
    for point in fixed_points:
        linearized = deriv.subs(x, point)
        tau = relaxation_time(xdot, point)
        dts[point] = tau / ONE_ORDER_OF_MAGNITUDE
        print(f"Fixed point: {point}, Linearized: {linearized}, tau: {tau}, dt: {dts[point]}, T: {integration_time_prediction(tau)}")
    return dts


def part_b(dts, xdot, x, r):
    initial_conditions = [-2, -0.5, 0.01, 0.5, 2]
    rvals = [1, -1]
    x_sym = sp.symbols('x')
    for rval in rvals:
        dts_sub = {point.subs(r, rval): time.subs(r, rval) for point, time in dts.items()}
        xdot_sub = xdot.subs(r, rval)
        f = sp.lambdify(x_sym, xdot_sub.rhs, modules='numpy')
        solutions = {}
        for start in initial_conditions:
            T = 8 * integration_time_prediction(relaxation_time(xdot_sub, start))
            dt = return_closest_dt(dts_sub, start)
            trajectory = integrate_forward_euler(f, start, dt, T)
            solutions[start] = trajectory
        plot_trajectories(solutions, "Trajectories starting at fixed points with predicted dt")


def plot_part_c(trajectories):
    fig, axs = plt.subplots(3, 2, figsize=(10, 15))
    for i, (dt, trajectory) in enumerate(trajectories.items()):
        times = np.arange(len(trajectory)) * dt
        row, col = divmod(i, 2)
        axs[row, col].plot(times, trajectory)
        axs[row, col].set_title(f"dt: {dt}")
        axs[row, col].set_xlabel('Time')
        axs[row, col].set_ylabel('x')
        axs[row, col].grid()
    plt.tight_layout()
    plt.show()


def part_c(dts, xdot, x, r):
    rval = 1
    x_sym = sp.symbols('x')
    xdot_sub = xdot.subs(r, rval)
    dts_sub = {point.subs(r, rval): time.subs(r, rval) for point, time in dts.items()}
    x0 = 0.01
    dt = return_closest_dt(dts_sub, x0)
    T = integration_time_prediction(relaxation_time(xdot_sub, x0))
    f = sp.lambdify(x_sym, xdot_sub.rhs, modules='numpy')
    trajectories = {}
    for _ in range(6):
        trajectory = integrate_forward_euler(f, x0, dt, T)
        trajectories[dt] = trajectory
        dt *= 2
        T *= 1.5
    plot_part_c(trajectories)
    print("relaxation time (tau):", relaxation_time(xdot_sub, x0))


def plot_scatter(final_values, title='Final Values vs Parameter r'):
    plt.scatter(list(final_values.keys()), list(final_values.values()))
    plt.xlabel('r')
    plt.ylabel('x(T)')
    plt.title(title)
    plt.grid()
    plt.show()


def _integrate_over_rvals(rvals, r, x, xdot, x0):
    x_sym = sp.symbols('x')
    deriv_expr = sp.diff(xdot.rhs, x_sym)
    final_values = {}
    for rval in rvals:
        xdot_subs = xdot.subs(r, rval)
        tau = 1 / abs(deriv_expr.subs(x_sym, x0).subs(r, rval))
        T = float(integration_time_prediction(tau))
        fixed_points = sp.solve(xdot_subs, x)
        dts = part_a(fixed_points, xdot_subs, r)
        dt = float(return_closest_dt(dts, x0))
        f = sp.lambdify(x_sym, xdot_subs.rhs, modules='numpy')
        integrated = solve_ivp(lambda t, y: float(f(y[0])), [0, T], [x0], method='RK45', max_step=dt)
        final_values[rval] = integrated.y[0][-1]
    return final_values


def part_d(r, x, xdot, dts):
    rvals = np.concatenate([np.linspace(-1, -0.1, 25), np.linspace(0.1, 1, 25)])
    x0 = norm.rvs(loc=0, scale=0.01)
    plot_scatter(_integrate_over_rvals(rvals, r, x, xdot, x0))


def part_e(r, x, xdot, dts):
    rvals = np.linspace(-0.1, 0.1, 50)
    x0 = norm.rvs(loc=0, scale=0.01)
    plot_scatter(_integrate_over_rvals(rvals, r, x, xdot, x0))


def main():
    r, x = sp.symbols('r x')
    xdot = sp.Eq(x, r*x - x**3)
    fixed_points = sp.solve(xdot, x)
    dts = part_a(fixed_points, xdot, r)
    # part_b(dts, xdot, x, r)
    # part_c(dts, xdot, x, r)
    # part_d(r, x, xdot, dts)
    part_e(r, x, xdot, dts)


if __name__ == "__main__":
    main()
