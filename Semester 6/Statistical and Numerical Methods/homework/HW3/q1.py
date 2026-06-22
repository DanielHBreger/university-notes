import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.optimize as optimize
import matplotlib.pyplot as plt

def main():
    data = pd.read_excel("meas2.xlsx")
    xs = data.iloc[:, 2].values
    ys = data.iloc[:, 3].values

    # fit gaussian and lorentzian
    gaussian = lambda x, A, mu, sigma: A * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))
    lorentzian = lambda x, A, mu, gamma: A * gamma**2 / ((x - mu) ** 2 + gamma**2)
    popt_gaussian, _ = optimize.curve_fit(gaussian, xs, ys, p0=[max(ys), np.mean(xs), np.std(xs)])
    popt_lorentzian, _ = optimize.curve_fit(lorentzian, xs, ys, p0=[max(ys), np.mean(xs), np.std(xs)])


    # calcualte likelihood functions
    residuals_gaussian = ys - gaussian(xs, *popt_gaussian)
    sigma_est = np.std(residuals_gaussian)
    likelihood_gaussian = np.prod(stats.norm.pdf(residuals_gaussian, scale=sigma_est))

    residuals_lorentzian = ys - lorentzian(xs, *popt_lorentzian)
    gamma_est = np.median(np.abs(residuals_lorentzian - np.median(residuals_lorentzian)))
    likelihood_lorentzian = np.prod(stats.cauchy.pdf(residuals_lorentzian, scale=gamma_est))

    print(f"Gaussian Likelihood: {likelihood_gaussian}")
    print(f"Lorentzian Likelihood: {likelihood_lorentzian}")

    #ratio of likelihoods    likelihood_ratio = likelihood_gaussian / likelihood_lorentzian
    print(f"Likelihood Ratio: {likelihood_gaussian / likelihood_lorentzian}")


if __name__ == "__main__":
    main()