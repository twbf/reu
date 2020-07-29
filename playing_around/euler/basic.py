import math
import numpy as np
import matplotlib.pyplot as plt

def f(x,y):
    return math.sqrt(y+1)

x0 = 0
y0 = 0
xend = 10
step = 0.01

numSteps = (xend - x0)/step

y = y0
x = x0

i=0

print("steps:")
print(numSteps)

y_historical = np.zeros([1001])

while i < numSteps:
    x += step
    y += step * f(x,y)
    i += 1
    y_historical[i] = y

print("value:")
print(y)

plt.plot([1, 2])
#plt.show()
