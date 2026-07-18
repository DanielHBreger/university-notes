from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
import scienceplots as sp
from scipy.optimize import curve_fit
plt.style.use(['science', 'notebook', 'grid'])

def main():
    a=1
    b_c=2
    xdot = lambda a,b,x,y: a-(b+1)*x + x**2*y
    ydot = lambda b,x,y: b*x - x**2*y
    t_span = (0, 50)
    t_eval = np.linspace(0, 50, 1000)
    bs = [1.8, 2.3]
    initial_conditions = [[1,1],[1.6,1.8]]
    sols = []
    for b in bs:
        for ic in initial_conditions:
            sol = solve_ivp(lambda t, z: [xdot(a, b, z[0], z[1]), ydot(b, z[0], z[1])], t_span, ic, t_eval=t_eval)
            sols.append(sol)
    for i, sol in enumerate(sols):
        plt.plot(sol.y[0], sol.y[1], label=f'b={bs[i//len(initial_conditions)]}, IC={initial_conditions[i%len(initial_conditions)]}')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Phase Space')
    plt.legend()
    plt.show()

    mus = [0.05, 0.1, 0.2, 0.4]
    amplitudes = []
    t_long = np.linspace(0, 1000, 20000)
    for mu in mus:
        sol = solve_ivp(lambda t, z: [xdot(a, b_c+mu, z[0], z[1]), ydot(b_c+mu, z[0], z[1])], (0, 1000), [1, 1], t_eval=t_long, rtol=1e-8, atol=1e-10)
        x_tail = sol.y[0][int(0.8*len(sol.t)):]  # discard the transient before measuring
        amplitudes.append((x_tail.max()-x_tail.min())/2)

    func = lambda mu, c, d: c * (mu**d)
    # fit

    popt, _ = curve_fit(func, mus, amplitudes)
    xs = np.linspace(0.05, 0.4, 100)
    ys = func(xs, popt[0], popt[1])

    plt.plot(mus, amplitudes, marker='o', label='Data', color='blue')
    plt.plot(xs, ys, label=f'Fit: c={popt[0]:.3f} d={popt[1]:.3f}', color='red')
    plt.xlabel('mu')
    plt.ylabel('Amplitude of x')
    plt.title('Amplitude vs mu')
    plt.legend()
    plt.show()

    # measure oscillation period
    mu = mus[0]  # only consider the first value of mus
    sol = solve_ivp(lambda t, z: [xdot(a, b_c+mu, z[0], z[1]), ydot(b_c+mu, z[0], z[1])], (0, 1000), [1, 1], t_eval=t_long, rtol=1e-8, atol=1e-10)
    x_tail = sol.y[0][int(0.8*len(sol.t)):]  # discard the transient before measuring
    peaks = (x_tail[1:-1] > x_tail[:-2]) & (x_tail[1:-1] > x_tail[2:])  # find local maxima
    peak_times = sol.t[int(0.8*len(sol.t)) + np.where(peaks)[0] + 1]  # get the times of the peaks
    periods = np.diff(peak_times)  # calculate the periods
    average_period = np.mean(periods)
    print(f'Average period for mu={mu}: {average_period}')

if __name__ == "__main__":
    main()