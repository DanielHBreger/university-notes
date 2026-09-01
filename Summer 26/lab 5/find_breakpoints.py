import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution
from prettytable import PrettyTable

# Read the trace data
data = pd.read_csv('trace1.csv')
ch1 = data['CH1(V)'].values
ch2 = data['CH2(V)'].values

resistance = 216  # ohms
current_full = (ch1 - ch2) / resistance

# Discard points outside this range (edge/saturation effects)
x_min, x_max = -9.33, 8.28
keep = (ch1 >= x_min) & (ch1 <= x_max)
ch1 = ch1[keep]
current = current_full[keep]

N_SEGMENTS = 5
N_BREAKPOINTS = N_SEGMENTS - 1  # interior breakpoints; outer edges are fixed
MIN_POINTS_PER_SEGMENT = 50


def build_segments(breakpoints):
    edges = np.concatenate(([x_min], np.sort(breakpoints), [x_max]))
    return list(zip(edges[:-1], edges[1:]))


def segment_fit(lo, hi):
    mask = (ch1 >= lo) & (ch1 <= hi)
    x_seg, y_seg = ch1[mask], current[mask]
    slope, intercept = np.polyfit(x_seg, y_seg, 1)
    y_pred = slope * x_seg + intercept
    ss_res = np.sum((y_seg - y_pred) ** 2)
    return x_seg, y_seg, slope, intercept, ss_res


MIN_SEGMENT_WIDTH = 0.3  # volts


def combined_r2(breakpoints):
    """R^2 of the full piecewise-linear model against the global variance."""
    ss_res_total = 0.0
    for lo, hi in build_segments(breakpoints):
        if hi - lo < MIN_SEGMENT_WIDTH:
            return None
        mask = (ch1 >= lo) & (ch1 <= hi)
        if mask.sum() < MIN_POINTS_PER_SEGMENT:
            return None
        try:
            _, _, _, _, ss_res = segment_fit(lo, hi)
        except np.linalg.LinAlgError:
            return None
        ss_res_total += ss_res
    ss_tot_total = np.sum((current - current.mean()) ** 2)
    return 1 - ss_res_total / ss_tot_total


def objective(breakpoints):
    r2 = combined_r2(breakpoints)
    return 1.0 if r2 is None else -r2  # differential_evolution minimizes


bounds = [(x_min + 0.5, x_max - 0.5)] * N_BREAKPOINTS

result = differential_evolution(
    objective, bounds, seed=0, maxiter=500, popsize=30,
    tol=1e-10, polish=True,
)

best_breakpoints = np.sort(result.x)
best_r2 = -result.fun
segments = build_segments(best_breakpoints)

print(f"Best combined R^2: {best_r2:.6f}")
print(f"Interior breakpoints (V): {np.round(best_breakpoints, 3).tolist()}")
print()
print("Paste into plot_current.py:")
print("segments = [")
for lo, hi in segments:
    print(f"    ({lo:.3f}, {hi:.3f}),")
print("]")
print()

# Per-segment statistics table
table = PrettyTable()
table.field_names = ["Segment (V)", "N", "Slope (S)", "Intercept (A)", "R²"]
seg_fits = []
for lo, hi in segments:
    x_seg, y_seg, slope, intercept, ss_res = segment_fit(lo, hi)
    ss_tot = np.sum((y_seg - y_seg.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot
    seg_fits.append(slope)
    table.add_row([f"[{lo:.3f}, {hi:.3f}]", len(x_seg), f"{slope:.4e}",
                    f"{intercept:.4e}", f"{r2:.4f}"])
print(table)

# --- Chua diode design condition: 1/|Ga| < R_total < 1/|Gb| ---
# Ga = inner (steepest, near v=0) slope; Gb = the two flanking "shoulder"
# slopes averaged (the outer, pre-saturation segments), per the standard
# Chua diode piecewise-linear model: i = Gb*v + 0.5*(Ga-Gb)*(|v+Bp|-|v-Bp|)
Ga = seg_fits[2]
Gb = 0.5 * (seg_fits[1] + seg_fits[3])

R0 = 990  # ohms

r_low, r_high = sorted([1 / abs(Ga), 1 / abs(Gb)])

print()
print(f"Ga (inner slope)  = {Ga:.6e} S")
print(f"Gb (outer slope)  = {Gb:.6e} S")
print(f"Required range: 1/|Ga| < R_total < 1/|Gb|  =>  {r_low:.2f} ohm < R_total < {r_high:.2f} ohm")

# Plot the optimized piecewise-linear fit over the raw data
plt.figure(figsize=(10, 6))
plt.plot(ch1, current, 'b-', linewidth=0.5, alpha=0.25, label='Raw')
for lo, hi in segments:
    x_seg, y_seg, slope, intercept, _ = segment_fit(lo, hi)
    x_fit = np.array([lo, hi])
    plt.plot(x_fit, slope * x_fit + intercept, color='crimson', linewidth=2.5)
for bp in best_breakpoints:
    plt.axvline(bp, color='gray', linestyle='--', linewidth=0.8)

plt.xlabel('CH1 Voltage (V)', fontsize=12)
plt.ylabel('Current (A)', fontsize=12)
plt.title(f'Optimized Breakpoints (combined R² = {best_r2:.4f})', fontsize=14)
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('optimized_breakpoints.png', dpi=150)
plt.show()
