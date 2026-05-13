import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erf

def plot_profile_at_times(times, x, T):
    plt.figure(figsize=(10, 6))
    for t in times:
        plt.plot(x, T(x, t), label=f't={t:.2f}')
    plt.xlabel('Position (x)')
    plt.ylabel('Temperature (T)')
    plt.title('Temperature Profile at Different Times')
    plt.legend()
    plt.grid()
    plt.show()

def main():
    # temperature and coefficient are both assumed to be 1 for simplicity
    # the position also isn't really the real position, scaling will be off but I imagine it's not very important for this.
    held_temp = 1
    xi = lambda x,t: x/np.sqrt(t)
    solution = lambda x,t: held_temp*(1-erf(xi(x,t)/2)) 
    times = np.linspace(0, 10, 10)
    x = np.linspace(0, 10, 1000)
    plot_profile_at_times(times, x, solution)

if __name__ == "__main__":
    main()