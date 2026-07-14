import numpy as np
import pandas as pd
import prettytable
from scipy.stats import chi2 as chi2_dist, norm
import matplotlib.pyplot as plt

def background_stats(baseline):
    return np.mean(baseline), np.std(baseline, ddof=1)

def check_background_gaussianity(baseline, plot=True):
    '''
    Check if the baseline noise is Gaussian using a chi-squared test
    '''
    mean, std_dev = background_stats(baseline)
    # histogram bin edges: 12 bins from mean-2.5*std to mean+2.5*std, plus two extra bins for the tails
    edges = np.concatenate(([-np.inf], np.linspace(mean - 2.5 * std_dev, mean + 2.5 * std_dev, 12), [np.inf]))
    observed, _ = np.histogram(baseline, bins=edges)
    expected = len(baseline) * np.diff(norm.cdf(edges, mean, std_dev))
    chi_squared = np.sum((observed - expected) ** 2 / expected)
    dof = len(observed) - 1 - 2  # bins - 1 - two estimated parameters (mean, std)
    reduced_chi_squared = chi_squared / dof
    p_value = chi2_dist.sf(chi_squared, dof)

    if plot:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(8, 5))
        plt.hist(baseline, bins=list(edges[1:-1]), density=True, alpha=0.5, label='Baseline Histogram')
        x = np.linspace(mean - 4 * std_dev, mean + 4 * std_dev, 100)
        plt.plot(x, norm.pdf(x, mean, std_dev), 'r-', label='Fitted Gaussian')
        plt.title('Baseline Histogram and Fitted Gaussian')
        plt.xlabel('Value')
        plt.ylabel('Density')
        plt.legend()
        plt.show()

    return chi_squared, reduced_chi_squared, dof, p_value

def get_signal(t, y, window=100, plot=True):
    '''
    Smooth the signal with a moving average and compute the noise
    '''
    # A convolution with a constant kernel is a moving average.
    ma = lambda x, w: np.convolve(x, np.ones(w) / w, mode='same')
    smooth = ma(y, window)
    # ignore the edges where the moving average is not well-defined
    m = window // 2  
    residual = (y - smooth)[m:-m] / np.sqrt(1 - 1 / window)
    t_r, smooth_r = t[m:-m], smooth[m:-m]
    slope = np.gradient(smooth_r, t_r)

    if plot:
        plt.figure(figsize=(10, 6))
        plt.subplot(2, 1, 1)
        plt.plot(t, y, label='Original Signal', alpha=0.5)
        plt.plot(t_r, smooth_r, label='Smoothed Signal', color='orange')
        plt.title('Signal and Smoothed Signal')
        plt.xlabel('Time (ns)')
        plt.ylabel('Signal')
        plt.legend()

        plt.subplot(2, 1, 2)
        plt.plot(t_r, residual, label='Residual (Noise)', color='green')
        plt.title('Residual (Noise) after Smoothing')
        plt.xlabel('Time (ns)')
        plt.ylabel('Residual')
        plt.legend()
        plt.tight_layout()
        plt.show()

    return t_r, smooth_r, slope, residual

def noise_during_signal(t, y, window=20, n_groups=6):
    '''
    get the residual noise during the signal by subtracting the smoothed signal from the raw data, then check if its gaussian
    '''
    t_r, smooth_r, slope, residual = get_signal(t, y, window)

    # mask for the signal region
    sig = t_r > 0
    chi_squared, reduced_chi_squared, dof, p_value = check_background_gaussianity(residual[sig], False)
    print(f"Residual noise inside signal region: chi2 = {chi_squared:.2f}, chi2/dof = {reduced_chi_squared:.2f}, dof = {dof}, p-value = {p_value:.3g}")
    # print noise stats
    mean, std_dev = background_stats(residual[sig])
    print(f"Residual noise inside signal region: mean = {mean:.3f}, std = {std_dev:.3f}, n = {len(residual[sig])}")
    # plot the histogram of the residual noise inside the signal region
    plt.figure(figsize=(8, 5))
    plt.hist(residual[sig], bins=30, density=True, alpha=0.5, label='Residual Noise Histogram')
    mean, std_dev = background_stats(residual[sig])
    x = np.linspace(mean - 4 * std_dev, mean + 4 * std_dev, 100)
    plt.plot(x, norm.pdf(x, mean, std_dev), 'r-', label='Fitted Gaussian')
    plt.title('Residual Noise Histogram and Fitted Gaussian (Signal Region)')
    plt.xlabel('Residual Value')
    plt.ylabel('Density')
    plt.legend()
    plt.show()


def extract_signal(t, y, bin_width=0.5, plot=True):
    '''
    put signal and its errors into bins, subtract the baseline noise, and compute the error on the mean for each bin.
    '''
    baseline = y[t <= -10]
    # baseline noise stats
    b_mean, b_std = background_stats(baseline)
    b_mean_err = b_std / np.sqrt(len(baseline))
    edges = np.arange(t.min(), t.max() + bin_width, bin_width)
    t_bin, signal, sigma = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (t >= lo) & (t < hi)
        n = m.sum()
        t_bin.append(t[m].mean())
        signal.append(y[m].mean() - b_mean)
        n_err = np.std(y[m], ddof=1) / np.sqrt(n)
        sigma.append(np.sqrt(n_err ** 2 + b_mean_err ** 2))
    t_bin, signal, sigma = np.array(t_bin), np.array(signal), np.array(sigma)
    print(f"binned signal: {len(t_bin)} bins of {bin_width} ns")

    if plot:
        plt.figure(figsize=(10, 5))
        plt.plot(t, y - b_mean, '.', ms=2, alpha=0.3, label='Raw data (background subtracted)')
        plt.errorbar(t_bin, signal, yerr=sigma, fmt='.', ms=4, color='#c5493d', elinewidth=1, label='Extracted signal S = N - B')
        plt.title('Extracted signal with errors')
        plt.xlabel('Time (ns)')
        plt.ylabel('Amplitude')
        plt.legend()
        plt.show()

    return t_bin, signal, sigma

def part_b(file_path):
    worksheet_name = 'מדידה אמיתית'
    data = pd.read_excel(file_path, sheet_name=worksheet_name)
    t = data.iloc[:, 0].values * 1e9  # type: ignore 
    y = data.iloc[:, 1].values
    baseline = y[t <= -10]  # the signal only starts around t = 0
    background_mean, background_std = background_stats(baseline)
    background_mean_err = background_std / np.sqrt(len(baseline))
    print(f"Baseline ({len(baseline)} samples): mean B = {background_mean:.3f} +- {background_mean_err:.3f}, std = {background_std:.3f}")
    print(f"Background Gaussianity Test:")
    chi_squared, reduced_chi_squared, dof, p_value = check_background_gaussianity(baseline, True)
    print(f"chi2 = {chi_squared:.2f}, chi2/dof = {reduced_chi_squared:.2f}, dof = {dof}, p-value = {p_value:.3g}")
    noise_during_signal(t, y)
    # estimate the signal and its errors (the dataset part c fits against)
    t_bin, signal, sigma = extract_signal(t, y)
    base, sig_region = t_bin <= -10, t_bin > 0
    i, j = np.argmin(signal), np.argmax(signal)
    print(f"signal minimum: S = {signal[i]:.0f} +- {sigma[i]:.0f} at t = {t_bin[i]:.2f} ns, significance S/sigma = {abs(signal[i]) / sigma[i]:.0f} sigma")
    print(f"signal maximum: S = {signal[j]:.0f} +- {sigma[j]:.0f} at t = {t_bin[j]:.2f} ns, significance S/sigma = {abs(signal[j]) / sigma[j]:.0f} sigma")
    print(f"median bin error: baseline {np.median(sigma[base]):.1f}, signal region {np.median(sigma[sig_region]):.1f}")
    return t_bin, signal, sigma

if __name__ == "__main__":
    part_b('Exam1.xlsx')
