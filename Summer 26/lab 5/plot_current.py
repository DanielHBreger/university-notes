import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from prettytable import PrettyTable

# Read the trace data
data = pd.read_csv('trace1.csv')

# Extract columns
time = data['Time(s)'].values
ch1 = data['CH1(V)'].values
ch2 = data['CH2(V)'].values


def block_average(arr, window):
    """Average consecutive non-overlapping windows of `window` samples."""
    n = (len(arr) // window) * window
    return arr[:n].reshape(-1, window).mean(axis=1)


window = 50  # number of samples per averaging window
ch1_avg = block_average(ch1, window)
ch2_avg = block_average(ch2, window)

# Current through the resistor, from the voltage drop across it
voltage_diff = ch1_avg - ch2_avg
resistance = 216  # ohms
current = voltage_diff / resistance

segments = [
    (-9.33, -6),
    (-6, -0.7),
    (-0.7, 0.5),
    (0.5, 5.3),
    (5.3, 8.28)
]

current_raw = (ch1 - ch2) / resistance

# Linear fit for each segment
fits = []
table = PrettyTable()
table.field_names = ["Segment (V)", "N", "Slope (S)", "Slope err (S)",
                      "Intercept (A)", "Intercept err (A)", "R²"]

for x_min, x_max in segments:
    mask = (ch1 >= x_min) & (ch1 <= x_max)
    x_seg, y_seg = ch1[mask], current_raw[mask]

    (slope, intercept), cov = np.polyfit(x_seg, y_seg, 1, cov=True)
    slope_err, intercept_err = np.sqrt(np.diag(cov))

    y_pred = slope * x_seg + intercept
    ss_res = np.sum((y_seg - y_pred) ** 2)
    ss_tot = np.sum((y_seg - y_seg.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot

    fits.append((x_min, x_max, slope, intercept))
    table.add_row([
        f"[{x_min:.2f}, {x_max:.2f}]",
        len(x_seg),
        f"{slope:.4e}",
        f"{slope_err:.2e}",
        f"{intercept:.4e}",
        f"{intercept_err:.2e}",
        f"{r2:.4f}",
    ])

print(table)

# Plot current against CH1 (voltage across the nonlinear resistor / Chua diode)
plt.figure(figsize=(10, 6))
plt.plot(ch1, current_raw, 'b-', linewidth=0.5, alpha=0.25, label='Raw')
plt.plot(ch1_avg, current, color='darkblue', linewidth=1.0, label=f'{window}-point average')

for x_min, x_max, slope, intercept in fits:
    x_fit = np.array([x_min, x_max])
    plt.plot(x_fit, slope * x_fit + intercept, color='crimson', linewidth=2.5)

plt.xlabel('CH1 Voltage (V)', fontsize=12)
plt.ylabel('Current (A)', fontsize=12)
plt.title('Current vs CH1 Voltage (Chua Diode Characteristic)', fontsize=14)
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()

# Save and show the plot
plt.savefig('current_vs_voltage.png', dpi=150)
plt.show()

# Print overall statistics
overall = PrettyTable()
overall.field_names = ["Quantity", "Min", "Max"]
overall.add_row(["Voltage difference (V)", f"{voltage_diff.min():.6f}", f"{voltage_diff.max():.6f}"])
overall.add_row(["Current (A)", f"{current.min():.6e}", f"{current.max():.6e}"])
print(overall)
