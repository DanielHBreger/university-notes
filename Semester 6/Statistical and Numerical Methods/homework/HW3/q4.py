import numpy as np
import scipy.stats as stats

def main():
    # country: (population, scholarship winners per 100k, error on number of winners)
    data = {
        "Israel": (9e6,3.6, 0.2),
        "Germany": (82.8e6,4.0,None),
        "Denmark": (5.8e6,3.5,None),
        "USA": (325.7e6,10.9,None),
        "Sri Lanka": (105e3,1.9,None),
        "Ivory Coast": (33e3,3.2,None),
        "Puerto Rico": (22e3,4.8,4.8),
    }
    lambda_calc = lambda pop, winners_per_100k: winners_per_100k * (pop/1e5)
    for country, (population, winners_per_100k, error) in data.items():
        lam = lambda_calc(population, winners_per_100k)
        standard_error = np.sqrt(lam)/(population/1e5) if error is None else error
        print(f"{country}: lambda = {lam:.2f}, standard error = {standard_error:.2f}")

if __name__ == "__main__":
    main()