
import numpy as np

t_steps = 5
x_steps = 5
t = 0
endt = 1

h_array = np.zeros((t_steps,x_steps))
u_array = np.zeros((t_steps,x_steps))

h_initial = np.array([1,0,1,1,1])
u_initial = np.array([1,0,1,1,1])

set initials

for i in range(tsteps):

    for j in range(xsteps - i): # -i because of boundry conditions

        h_t1 = (tstep/xstep)(h_x*u_x - hx1*ux1) + h_t
