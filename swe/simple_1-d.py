
# to start venv in Github put source ./venv/bin/activate

import numpy as np

print("simple 1-d swe simulation")
print("initializing......")

t_steps = 10
x_steps = 10
t = 0
g = -9.81
endt = 1

h_array = np.zeros((t_steps,x_steps))
u_array = np.zeros((t_steps,x_steps))

h_initial = np.array([1,0,1,1,1,0,0,0,0,0])
u_initial = np.array([1,0,1,1,1,0,0,0,0,0])

print("initial h:")
print(h_initial)
print("initial u:")
print(u_initial)

h_array[0] = h_initial
u_array[0] = u_initial

print("finite difference method ......")

for i in range(t_steps-1):

    h = h_array[i]
    u = u_array[i]

    h_new = h_array[i+1]
    u_new = u_array[i+1]

    for j in range(x_steps - i-1): # -i because of boundry finite derivitives

        h_new[j] = h[j] + (t_steps/x_steps)*(h[j]*u[j] - h[j+1]*u[j+1])
        u_new[j] = u[j] - (t_steps)*(u[j]*(u[j+1]-u[j])/x_steps + g*(h[j+1]-h[j])/x_steps)

    h_array[i+1] = h_new
    u_array[i+1] = u_new


print("h array:")
print(h_array)
