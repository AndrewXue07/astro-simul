import numpy as np

def main(n, d):
    n_circle = 0
    n_total = 0
    dist = 0 

    for i in range(n):
        
        for j in range(d):
            dist += np.random.uniform(-1, 1) **2 # distance formula in 10 dimensions

        if dist <= 1:
            n_circle += 1

        n_total += 1
        dist = 0

    print((n_circle/n_total) * 2 **d) # cube in ten dimension volume = 2^10

if __name__ == "__main__":
    main(100000, 6)
