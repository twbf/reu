import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

a = 1.0
beta = 2.0
v = 2.0
u = 1.0

bigA = v*u**(3/2)*(1+1j)/(2*np.sqrt(a))

b = -beta*u/np.sqrt(a)+2j*np.sqrt(a)

def mheight(z,t):
    tmp1 = bigA*(t+1j*b)
    tmp2 = (z-((t+1j*b)**2)/4)**(3/2)
    return tmp1/tmp2

def mspeed(z,t):
    tmp1 = bigA
    tmp2 = (z-((t+1j*b)**2)/4)**(3/2)
    return 2*tmp1/tmp2

t = 20
x=500

height = np.zeros((x,t))
speed = np.zeros((x,t))

for j in range(t):
    for i in range(x):
        height[i][j] = mheight(i/5,j)
        speed[i][j] = mspeed(i/5,j)


plt.plot(height)
plt.show()
plt.plot(speed)
plt.show()
