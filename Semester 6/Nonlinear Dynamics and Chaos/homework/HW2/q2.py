import numpy as np

def main():
    x0 = 1
    next = lambda x: x/2 + 2/x
    for i in range(6):
        print(x0)
        x0 = next(x0)

if __name__ == "__main__":
    main()