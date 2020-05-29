#L2.py

# imports space seperated lists from "analytic.txt" and "numerical.txt"
# computes L2 norm, then plots on the same axis

#!!!! currently only works with 1 timestep

import numpy as np
import matplotlib.pyplot as plt

print ("reading in files....")

anaylitic = np.loadtxt("analytic.txt")

numerical = np.loadtxt("numerical.txt")

if len(anaylitic) != len(numerical):
    print("ERROR: arrays must be same length")
else:
    print("success")

print("computing L2 norm...")
#the differnce between anaylitic and numerical
residual = numerical-anaylitic #order doesn't matter becasue we are taking the square

#L2 norm - sqrt of the sum of the sqaures of the residual
l2 = np.sqrt(np.dot(residual,residual.T))/len(residual)
print("success")

print("\nL2 norm:")
print(l2)

#plotting
plt.plot(anaylitic, label="Analytic")
plt.plot(numerical, label="Numerical")
plt.legend()
plt.show()
