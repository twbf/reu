import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

a = 1.0
beta = 2.0
v = 2.0
u = 1.0

bigA = v*u**(3/2)*(1+1j)/(2*np.sqrt(a))

b = -beta*u/np.sqrt(a)+2j*np.sqrt(a)

xtt = np.arange(0, 500, .1)
fig, ax = plt.subplots()
line, = ax.plot(xtt, np.sin(xtt))

def mheight(z,t):
    tmp1 = bigA*(t+1j*b)
    tmp2 = (z-((t+1j*b)**2)/4)**(3/2)
    return tmp1/tmp2


def heightgen(t):
    x=5000
    height = np.zeros((x))
    for i in range(x):
        height[i] = mheight(i/5,t/10)
    line.set_ydata(height)
    return line,

def init():
    line.set_ydata([np.nan] * len(xtt))
    return line,

def mspeed(z,t):
    tmp1 = bigA
    tmp2 = (z-((t+1j*b)**2)/4)**(3/2)
    return 2*tmp1/tmp2


ani = animation.FuncAnimation(
    fig, heightgen, init_func=init, interval=1, blit=True)
plt.show()
