import numpy as np
from scipy.stats import chi2 as chi2_dist, norm
import matplotlib.pyplot as plt
import pandas as pd

'''
Model Functions
'''
def y1(x,A,B,C,D):
    return A*np.sin(B*x)+C*np.sin(2*B*x)+D

def y2(x,A,B,C,D):
    return A*np.exp(-x/B)*np.sin(C*x)+D

def y3(x,A,B,C,D):
    return A*(x**3)+B*(x**2)+C*x+D


def load_excel(file_path, worksheet_name):
    data = pd.read_excel(file_path, sheet_name=worksheet_name)
    times = data.iloc[:, 0].values
    amplitudes = data.iloc[:, 1].values
    time_errors = data.iloc[:, 2].values
    amplitude_errors = data.iloc[:, 3].values
    data = {
        "times": times,
        "amplitudes": amplitudes,
        "time_errors": time_errors,
        "amplitude_errors": amplitude_errors
    }
    return data

def chi2(params, model, t, y, sigma):
    A, B, C, D = params
    model_values = model(t, A, B, C, D)
    chi_squared = np.sum(((y - model_values) / sigma) ** 2)
    return chi_squared

def model_slope(model, params, t):
    '''
    Calculate the time derivative of the model
    '''
    eps = 1e-6 * (t.max() - t.min()) # step size chosen arbitrarily
    return (model(t + eps, *params) - model(t - eps, *params)) / (2 * eps)

def check_fit_quality(chi_squared, n_points, n_params):
    dof = n_points - n_params
    reduced_chi_squared = chi_squared / dof
    p_value = chi2_dist.sf(chi_squared, dof)
    return reduced_chi_squared, p_value

def solve_linear(basis_columns, y, sigma):
    '''
    Weighted least squares
    '''
    # turn the list of basis columns into a matrix
    X = np.column_stack(basis_columns)
    # transpose and weight by the inverse variance
    XtW = X.T / sigma**2
    # solve the equation (X^T W X) p = X^T W y
    return np.linalg.solve(XtW @ X, XtW @ y)

def nonlinear_ranges(model, t):
    '''
    scan the nonlinear parameter(s) over a grid of reasonable values
    '''
    span = t.max() - t.min()
    freq = (2 * np.pi / (10 * span), 2 * np.pi / (span / 10))
    if model is y2:
        return [(span / 1000, 1000 * span), freq]  # (decay, frequency)
    return [freq]

def profile_params(model, nonlinear, t, y, sigma):
    '''
    Solve for the linear parameters given the nonlinear ones.
    '''
    ones = np.ones_like(t)
    if model is y1:
        (freq,) = nonlinear
        A, C, D = solve_linear([np.sin(freq * t), np.sin(2 * freq * t), ones], y, sigma)
        return np.array([A, freq, C, D])
    decay, freq = nonlinear  # y2: B is the decay time, C the frequency
    A, D = solve_linear([np.exp(-t / decay) * np.sin(freq * t), ones], y, sigma)
    return np.array([A, decay, freq, D])

def fit_model(model, t, y, sigma, num_points=100, steps=4):
    '''
    Fit the model to the data by scanning the nonlinear parameters, solving for the linear parameters, and returning the best-fit parameters and the corresponding chi-squared.
    '''
    if model is y3:
        # model is linear, just solve
        params = solve_linear([t**3, t**2, t, np.ones_like(t)], y, sigma)
        return params, chi2(params, model, t, y, sigma)
    axes = nonlinear_ranges(model, t)
    limits = list(axes)
    best_params, best_chi2 = None, np.inf
    for _ in range(steps):
        # create a grid of nonlinear parameters
        grids = np.meshgrid(*[np.linspace(lo, hi, num_points) for lo, hi in axes])
        best_nl = None
        # iterate over the grid of nonlinear parameters.
        # g.ravel() flattens the grid to a 1D array, and zip() combines the flattened arrays into tuples of nonlinear parameters.
        for nonlinear in zip(*(g.ravel() for g in grids)):
            params = profile_params(model, nonlinear, t, y, sigma)
            c = chi2(params, model, t, y, sigma)
            if c < best_chi2:
                best_params, best_chi2, best_nl = params, c, nonlinear
        if best_nl is None:
            break  # no improvement at this resolution, we're done
        # zoom in on the best nonlinear parameters
        step_sizes = [(hi - lo) / (num_points - 1) for lo, hi in axes]
        axes = [(max(v - 2 * s, lim[0]), min(v + 2 * s, lim[1])) for v, s, lim in zip(best_nl, step_sizes, limits)]
    if best_params is None:
        raise RuntimeError("Grid scan found no finite chi2")
    return best_params, best_chi2

def fit_with_t_errors(model, t, y, sigma_t, sigma_y, max_iter=10, rtol=1e-6):
    '''
    the time errors mean the simple chi2 doesn't apply, so we iterate: we start with the amplitude errors only, then calculate the effective sigma including the time errors, and repeat until convergence.
    '''
    params, chi_squared = fit_model(model, t, y, sigma_y)
    sigma_eff = sigma_y
    for _ in range(max_iter):
        slope = model_slope(model, params, t)
        sigma_eff = np.sqrt(sigma_y**2 + (slope * sigma_t)**2)
        new_params, new_chi2 = fit_model(model, t, y, sigma_eff)
        # if chi2 changes by a very small amount, we converged
        converged = abs(new_chi2 - chi_squared) < rtol * max(chi_squared, 1.0)
        params, chi_squared = new_params, new_chi2
        if converged:
            break
    return params, chi_squared, sigma_eff