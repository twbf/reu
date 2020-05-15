
# to start venv in Github put source ./venv/bin/activate

import numpy as np
import matplotlib.pyplot as plt
import Math

print("simple 1-d swe simulation")
print("initializing......")

num_t = 10
num_x = 40

t = 0
g = 9.81
endt = 0.01
endx = 1

t_steps = endt/num_t
x_steps = endx/num_x


h_array = np.zeros((num_t,num_x))
u_array = np.zeros((num_t,num_x))

h_initial = np.array([1,1,1,1,1,1,1,1,1,1,1,1,1.1,1.2,1.3,1.5,1.75,1.9,1.9,1.9,1.75,1.5,1.25,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
u_initial = np.array([0,0,0,0,0,0,0,0,0,0,0,0,-0.5,-1,-1,-1,-1,-0.5,0,0.5,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

print("initial h:")
print(h_initial)
print("initial u:")
print(u_initial)

h_array[0] = h_initial
u_array[0] = u_initial

print("finite difference method ......")

for i in range(num_t-1):

    h = h_array[i]
    u = u_array[i]

    h_new = h_array[i+1]
    u_new = u_array[i+1]

    for j in range(i+1, num_x - i-1): # -i because of boundry finite derivitives

        h_new[j] = h[j] + (t_steps/x_steps)*(h[j-1]*u[j-1] - h[j+1]*u[j+1])
        u_new[j] = u[j] - (t_steps)*(u[j]*(u[j+1]-u[j-1])/x_steps + g*(h[j+1]-h[j-1])/x_steps)

    h_array[i+1] = h_new
    u_array[i+1] = u_new


print("h array:")
print(h_array)

print(" ")
print(" ")

print("U array:")
print(u_array)

for i in range(num_t):
    plt.plot(h_array[i])
    plt.show()


print("Done")
