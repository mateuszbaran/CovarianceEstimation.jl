from sklearn.covariance import oas, ledoit_wolf
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

times = np.zeros(len(p))
res = np.zeros(len(p))
for i in range(len(p)):
    Xi = X[i]
    Ci = C[i]
    start = time()
    C_oas, _ = oas(Xi)
    times[i] = time() - start
    res[i] = np.linalg.norm(C_oas - Ci)

print("OAS results")
print(res)
print("")
print(times)

times = np.zeros(len(p))
res = np.zeros(len(p))
for i in range(len(p)):
    Xi = X[i]
    Ci = C[i]
    start = time()
    C_lw, _ = ledoit_wolf(Xi)
    times[i] = time() - start
    res[i] = np.linalg.norm(C_lw - Ci)

print("LW results")
print(res)
print("")
print(times)
