import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.integrate as integrate
import scipy.special as special

#notes:
#assumption of 0 initial velocity
#phi(velocity) isn't implemeented yet
#uses numerical integration to infinity which is slow

m = 1000
beta = np.sqrt(m/(m+1))


def initial_eta(x):
    return np.e**(-(x-5)**2)

def initial_eta_prime(x):
    return -2*(x-5)*np.e**(-(x-5)**2)

def initial_u(x):
    return 0

def b(k):
    return 0

def j(x,k):
    return special.jv(1/m,2*k*np.sqrt(x+initial_eta(x)))

def a(k):
    result = integrate.quad(lambda x: initial_eta(x)*(x+initial_eta(x))*j(x,k)*(1+initial_eta_prime(x)), 0, 20)
    return 2*k*result[0]

def wave_plus(k,lambd): #combination of Sin's and Cos's
    return a(k)*np.cos(beta*k*lambd)+b(k)*np.sin(beta*k*lambd)

def wave_minus(k,lambd): #combination of Sin's and Cos's
    return a(k)*np.cos(beta*k*lambd)-b(k)*np.sin(beta*k*lambd)

def psi(s,lambd):  #height
    result = integrate.quad(lambda k: wave_plus(k,lambd)*special.jv(1/m,2*k*np.sqrt(s)), 0, 20) #cuts at 20
    return s**(-1/(2*m))*result[0]

def phi(s,lambd):  #velocity
    result = integrate.quad(lambda k: wave_minus(k,lambd)*special.jv(1/m+1,2*k*np.sqrt(s)), 0, 20) #cuts at 20
    return 1/beta*s**(-1/(2*m)-1/2)*result[0]

def greenspan(s,lambd):
    u = phi(s,lambd)
    eta = psi(s,lambd)+u**2/2
    #x = s-eta
    #t = lambd+u
    return eta



ndx = 2
dx = 1
ndt = 2
dt = 0.1


#x = np.zeros((ndt, ndx))
x1 = np.zeros((ndt, ndx))
for c in range(ndt):
    for i in range(ndx):
        print(c,i)
        #x[c][i] = greenspan(i*dx + 1,c*dt+0.6)
        x1[c][i] = psi(i*dx + 1,c*dt+ 0.6)  # +0.001 is becasue solution breaks down at 0
    #plt.plot(x[c])
    plt.plot(x1[c])

#print(x)
print(x1)
#plt.set_ylabel('psi')
#plt.set_xlabel('x')
plt.show()
