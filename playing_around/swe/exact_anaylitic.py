import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

a = 0.1
beta = 1
v = 5.0
u = 10.0

bigA = v*u**(3/2)*(1+1j)/(2*np.sqrt(a))

b = -beta*u/np.sqrt(a)+2j*np.sqrt(a)

def bigN(z,t):
    tmp1 = bigA*(t+1j*b)
    tmp2 = (z-((t+1j*b)**2)/4)**(3/2)
    return tmp1/tmp2

def bigU(z,t):
    tmp1 = bigA
    tmp2 = (z-((t+1j*b)**2)/4)**(3/2)
    return 2*tmp1/tmp2

def mheight(x,t):  #this might be wrong - change of variable
    return bigN(x,t) -(bigU(x,t)**2)/2

def mspeed(x,t): #this might be wrong - change of variable
    return bigU(x,t)

def transform(z,tao):
    return bigN(z,tao)-(bigU(z, tao)**2)/2 , tao + bigU(z, tao)


def heightgen(tao):
    steps=500
    height = np.zeros((steps))
    height1 = np.zeros((steps))
    for z in range(steps):
        x,t = transform(z-110,tao-50)
        height[z] = bigN(z-100,tao-50)
        height1[z] = x
    line.set_ydata(height)
    line2.set_ydata(height1)
    return line, line2,

def init():
    line2.set_ydata([np.nan] * len(xtt))
    line.set_ydata([np.nan] * len(xtt))
    return line, line2,

xtt = np.arange(-20, 480, 1)
fig, ax = plt.subplots()
line, = ax.plot(xtt, 5*np.sin(xtt))
line2, = ax.plot(xtt, np.sin(xtt))

Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)


ani = animation.FuncAnimation(
    fig, heightgen, init_func=init, interval=1, blit=True)
ani.save('im.mp4', writer=writer)
plt.show()
