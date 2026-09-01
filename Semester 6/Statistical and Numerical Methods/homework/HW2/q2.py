import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scienceplots as sp
import scipy.stats as stats

plt.style.use(['science', 'no-latex'])

def part_a():
    # 10000 measurements from one voltmeter
    number_of_samples = 10000
    measurements = np.random.normal(loc=0, scale=0.01, size=number_of_samples)
    weights = np.ones(number_of_samples) / number_of_samples * 100
    plt.figure(figsize=(10, 6))
    plt.hist(measurements, bins=25, weights=weights)
    plt.xlabel("Measurement Value", fontsize=16)
    plt.ylabel("Percentage of Measurements", fontsize=16)
    plt.title("Histogram of Measurements", fontsize=20)
    plt.xlim(-0.05, 0.05)
    plt.show()

def part_b():
    # 50 measurements from 50 voltmeters 
    voltmeters = 50
    measurements_per_voltmeter = 1
    measurements = np.random.normal(loc=0, scale=0.01, size=(voltmeters, measurements_per_voltmeter))
    flat_measurements = measurements.flatten()
    plt.figure(figsize=(10, 6))
    plt.hist(flat_measurements, bins=10, weights=np.ones(flat_measurements.size) / flat_measurements.size * 100)
    plt.xlabel("Mean Measurement Value", fontsize=16)
    plt.ylabel("Percentage of Voltmeters", fontsize=16)
    plt.title("Histogram of Measurements from 50 Voltmeters", fontsize=20)
    plt.xlim(-0.05, 0.05)
    plt.show()

def part_c():
    deviations = [0.1, 0.05, 0.01, 0.02]
    measurements_per_device = [20, 10, 5, 15]
    measurements = []
    for deviation, num_measurements in zip(deviations, measurements_per_device):
        device_measurements = np.random.normal(loc=0, scale=deviation, size=num_measurements)
        measurements.append(device_measurements)
    flat_measurements = np.concatenate(measurements)
    # per-measurement sigma: each measurement inherits its device's sigma
    sigmas = np.concatenate([np.full(n, s) for s, n in zip(deviations, measurements_per_device)])
    inv_var = 1 / sigmas**2
    weighted_mean = np.sum(inv_var * flat_measurements) / np.sum(inv_var)
    sigma_mean = 1 / np.sqrt(np.sum(inv_var))
    print(f"Weighted Mean: {weighted_mean:.6f}")
    print(f"Uncertainty on Mean: {sigma_mean:.6f}")

    plt.figure(figsize=(10, 6))
    plt.hist(flat_measurements, bins=50, weights=np.ones(flat_measurements.size) / flat_measurements.size * 100)
    plt.xlabel("Measurement Value", fontsize=16)
    plt.ylabel("Percentage of Measurements", fontsize=16)
    plt.title("Histogram of Measurements from Multiple Devices", fontsize=20)
    plt.xlim(-0.3, 0.3)
    plt.show()
           
def part_d():
    voltages = [0.5, 1.0, 3.9, 7.8, 10, 13.4, 15.2, 17.5, 21.1, 25]
    measurements = []
    for voltage in voltages:
        device_measurements = np.random.normal(loc=voltage, scale=0.01, size=5)
        measurements.append(device_measurements)
    # linear fit and check how good the fit is
    all_measurements = np.array(measurements).flatten()
    all_voltages = np.repeat(voltages, 5)
    slope, intercept, r_value, p_value, std_err = stats.linregress(all_voltages, all_measurements)
    print(f"Linear Fit: y = {slope:.4f}x + {intercept:.4f}")
    print(f"slope uncertainty: {std_err:.4f}")
    print(f"intercept uncertainty: {std_err:.4f}")
    print(f"R-squared: {r_value**2:.10f}")
    # scatter plot of all the measurements vs true voltage grouped by measurement number at each voltage
    plt.figure(figsize=(10, 6))
    for i in range(5):
        plt.scatter(voltages, [measurements[j][i] for j in range(len(voltages))], label=f'Measurement {i+1}')
    plt.xlabel("True Voltage (V)", fontsize=16)
    plt.ylabel("Measured Voltage (V)", fontsize=16)
    plt.title("Measured Voltage vs True Voltage", fontsize=20)
    plt.legend()
    plt.grid()
    plt.show()

def main():
    part_a()
    part_b()
    part_c()
    part_d()

if __name__ == "__main__":
    main()