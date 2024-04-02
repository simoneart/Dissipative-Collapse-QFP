import numpy as np
import matplotlib.pyplot as plt
from functions import *
from matplotlib.ticker import MaxNLocator
from matplotlib import cm
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation

'''
Full QFP eq. in 1-D due to computational cost. We restrict to the CSL model,
the DP should be thermodynamically equivalent (rescale the paramenters using the
3-D counterpart)
'''

no_points = 200
L, R = -2., 2.

#for higher dimensions use a regular grid, no need to define a line for every direction
q_grid = np.linspace(L, R, no_points, dtype=np.float64)
p_grid = np.linspace(L, R, no_points, dtype=np.float64)
Grid = [q_grid, p_grid]
dx = np.float64((R - L)/no_points) #the grid is evenly spaced


#IC
var0 = 0.1
d = np.array([0,0])
V = np.array([[var0,0],[0,var0]], dtype=np.float64)
w0 = np.array([[gauss_Wf([q,p], d, V) for p in p_grid] for q in q_grid], dtype=np.float64)

w0 /= np.sum(w0 * dx * dx)

#time grid
dt = np.float64(0.001)
Num_steps = 1000
times = np.array([k*dt for k in range(Num_steps)], dtype=np.float64)

print('Starting the dynamics')

#solution of the dynamics
wt = dyn(w0,Grid,times,dx,dt)

print('End of the simulation')

'''
#target state of the dynamics
Vtarget = np.array([[0.2,0],[0,0.2]])
Wtarget = np.array([[gauss_Wf([q,p], d, Vtarget) for p in p_grid] for q in q_grid])

ent_prod = entropy_production(wt, Wtarget, dt, Num_steps, [q_grid,p_grid])

at = [times[i] for i in range(len(times)-1)]

plt.figure(1)
plt.grid()
plt.plot(at, ent_prod)
plt.xlabel('t')
plt.ylabel('Entropy production rate')
'''

def plotheatmap(u_k, k):
    # Clear the current plot figure
    plt.clf()

    plt.title(f"Wigner function at t = {k*dt:.3f} unit time")
    plt.xlabel("q")
    plt.ylabel("p")

    # This is to plot u_k (u at time-step k)
    plt.pcolormesh(q_grid, p_grid, u_k, cmap=plt.cm.jet)
    plt.colorbar()

    return plt


def animate(k):
    plotheatmap(wt[k], k)


print('Preparing the animation')

anim = animation.FuncAnimation(plt.figure(), animate, interval=1, frames=Num_steps, repeat=False)
anim.save("1D_1p_diff.gif")

