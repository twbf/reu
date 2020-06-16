
# to start venv in Github put source ./venv/bin/activate

import numpy as np
import matplotlib.pyplot as plt

print("simple 1-d swe simulation")
print("initializing......")

num_t = 60
num_x = 200

t = 0
g = 9.81
endt = 40
endx = 2

t_steps = endt/num_t
x_steps = endx/num_x

def initial_height(x):
    return np.exp(x+np.cos((4*x)**2))

def initial_speed(x):
    return 1


h_array = np.zeros((num_t,num_x))
u_array = np.zeros((num_t,num_x))

h_initial = np.zeros((num_x))
u_initial = np.zeros((num_x))

for i in range(1, num_x):
    h_initial[i] = initial_height(x_steps*i)
    u_initial[i] = initial_speed(x_steps*i)


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

        u_new[j] = u[j] - (t_steps)*(u_new[j-1]*(u[j+1]-u[j-1])/2*x_steps + g*(h[j+1]-h[j-1])/2*x_steps)
        h_new[j] = h[j] + (t_steps)*(u_new[j]*(h[j+1]-h[j-1])/2*x_steps + h_new[j-1]*(u[j+1]-u[j-1])/2*x_steps)


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
