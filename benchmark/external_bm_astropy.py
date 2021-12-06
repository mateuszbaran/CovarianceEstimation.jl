from astropy.stats import biweight_midcovariance
import numpy as np
from time import time

p = [20, 40, 200, 400]

C_20  = np.genfromtxt('tmpmat/C_20.csv')
C_40  = np.genfromtxt('tmpmat/C_40.csv')
C_200 = np.genfromtxt('tmpmat/C_200.csv')
C_400 = np.genfromtxt('tmpmat/C_400.csv')

C = [C_20, C_40, C_200, C_400]

X_40x20   = np.genfromtxt('tmpmat/X_40x20.csv')
X_20x40   = np.genfromtxt('tmpmat/X_20x40.csv')
X_400x200 = np.genfromtxt('tmpmat/X_400x200.csv')
X_200x400 = np.genfromtxt('tmpmat/X_200x400.csv')

X = [X_40x20, X_20x40, X_400x200, X_200x400]

# Ensure that numpy has fully initialised before timing anything.
C_bw = biweight_midcovariance(X[1].T)

times = np.zeros(len(p))
res = np.zeros(len(p))
for i in range(len(p)):
    Xi = X[i].T
    Ci = C[i]
    start = time()
    C_bw = biweight_midcovariance(Xi)
    times[i] = time() - start
    res[i] = np.linalg.norm(C_bw - Ci)

print("biweight_midcovariance results")
print(res)
print("")
print(times)