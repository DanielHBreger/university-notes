import pandas as pd
import numpy as np
import prettytable

from q4_common import y1, y2, y3, fit_model, check_fit_quality
from q4_part_b import extract_signal

def load_sheet(file_path, worksheet_name):
    data = pd.read_excel(file_path, sheet_name=worksheet_name)
    t = data.iloc[:, 0].values * 1e9  # type: ignore  # convert s to ns
    y = data.iloc[:, 1].values
    return t, y

def align_simulation(t_sim, y_sim, t_data, s_data, sigma_data, bin_width=0.5):
    '''
    The times in the simulation don't match the measurement times, so the simulation needs to be shifted to align. The correct shift is found by minimizing the chi-squared between the simulation and the measurement.
    '''
    # sort the simulation by time
    order = np.argsort(t_sim)
    t_sim, y_sim = t_sim[order], y_sim[order]
    # Since the simulation data has a different resolution compared to the measurement data, trying to fit one against the other will probably fail. because of that, the simulation is averaged in the same way as the measurement data
    sub = t_data[:, None] + np.linspace(-bin_width / 2, bin_width / 2, 11)[None, :]
    best = None
    for shift in np.arange(0.0, 40.0, 0.05):
        # find which measurement bins are fully covered by the simulation after the shift
        covered = ((t_data - bin_width / 2 >= t_sim.min() + shift) & (t_data + bin_width / 2 <= t_sim.max() + shift))
        # compute the mean simulation value
        prediction = np.interp(sub[covered] - shift, t_sim, y_sim).mean(axis=1)
        # chi-squared between the binned simulation and the measurement
        c = np.sum(((s_data[covered] - prediction) / sigma_data[covered]) ** 2)
        if best is None or c / covered.sum() < best[0]: 
            best = (c / covered.sum(), shift, covered, prediction)
    if best is None:
        raise RuntimeError("No offset gives enough overlap between simulation and measurement")
    _, shift, covered, prediction = best
    return shift, covered, prediction

def part_c(file_path):
    # real measurement data
    t_m, y_m = load_sheet(file_path, 'מדידה אמיתית')
    t_bin, s, sigma = extract_signal(t_m, y_m, plot=False)
    m = t_bin >= 0
    t_bin, s, sigma = t_bin[m], s[m], sigma[m]

    table = prettytable.PrettyTable()
    table.field_names = ["Hypothesis", "chi2", "chi2/dof", "p-value"]
    fits = {}
    for model, name in ((y1, "Y1"), (y2, "Y2"), (y3, "Y3")):
        params, chi_squared = fit_model(model, t_bin, s, sigma, num_points=200, steps=6)
        fits[name] = (model, params)
        red, p = check_fit_quality(chi_squared, len(t_bin), len(params))
        print(f"{name} best fit: params = {np.round(params, 4)}")
        table.add_row([name, f"{chi_squared:.0f}", f"{red:.2f}", f"{p:.11g}"])

    t_sim, y_sim = load_sheet(file_path, 'סימולציה')
    shift, covered, prediction = align_simulation(t_sim, y_sim, t_bin, s, sigma)
    chi_squared = np.sum(((s[covered] - prediction) / sigma[covered]) ** 2)
    red, p = check_fit_quality(chi_squared, covered.sum(), 1)  # 1 fitted offset
    print(f"simulation time offset: {shift:.2f} ns ({covered.sum()} bins covered)")
    table.add_row(["Simulation", f"{chi_squared:.0f}", f"{red:.2f}", f"{p:.11g}"])
    print(table)

    import matplotlib.pyplot as plt
    plt.figure(figsize=(11, 6))
    plt.errorbar(t_bin, s, yerr=sigma, fmt='.', ms=3, color='gray', elinewidth=0.8, label='Binned measurement')
    t_dense = np.linspace(t_bin.min(), t_bin.max(), 1000)
    for name, (model, params) in fits.items(): plt.plot(t_dense, model(t_dense, *params), lw=1.5, label=f'{name} fit')
    plt.plot(t_bin[covered], prediction, lw=1.5, ls='--', label=f'Simulation (shifted {shift:.1f} ns)')
    plt.xlabel('Time (ns)')
    plt.ylabel('Amplitude')
    plt.title('Hypotheses and simulation vs the measured signal')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    part_c("Exam1.xlsx")
