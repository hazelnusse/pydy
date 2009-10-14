#!/usr/bin/env python
import pendulum_lib as p
from scipy.integrate import odeint
from numpy import array, arange, zeros, pi
import matplotlib.pyplot as plt

def plot_energy(t, x):
    # Plot the kinetic energy, potential energy, and total energy
    ke = zeros((n,1))
    pe = zeros((n,1))
    te = zeros((n,1))
    for i in range(n):
        ke[i], pe[i] = p.energy(x[i,:], params)
        te[i] = ke[i] + pe[i]

    plt.figure()
    plt.plot(t, ke, label='KE')
    plt.plot(t, pe, label='PE')
    plt.plot(t, te, label='TE')
    plt.legend()
    plt.title('Energy of pendulum during integration\n m = %0.2f, g = %0.2f, l = %0.2f, b = %0.2f'%(m, g, l, b))
    plt.xlabel('Time [s]')
    plt.ylabel('Energy [kg * m ** 2 / s**2]')
    plt.show()
    ############################

m = 1.
g = 9.8
l = 1.
b = 0.
params = [m, g, l, b]
x0 = [pi/4., 0.0]

# Integration time
ti = 0.0
ts = 0.01
tf = 100.0
t = arange(ti, tf+ts, ts)
n = len(t)
# Integrate the differential equations
x = odeint(p.eoms, x0, t, args = (params,))
plot_energy(t, x)
stop

# Animate using Visual-Python
from visual import display, rate, arrow, cylinder, sphere
# Animation playback speed multiplier (1 == realtime)
k = 1.0

# Set up the window
scene = display(title='Point mass pendulum animation @ %0.2f realtime'%k,
        background=(0,0,0), forward=(1,0,0), up=(0,0,1),
        width=800, height=800)
black = (0,0,0)
red = (1, 0, 0)
green = (0, 1, 0)
blue = (0, 0, 1)
P = zeros((n,3))

# Call the animate function to get the A[3] unit vector
for i, state in enumerate(x[:,0]):
    P[i] = p.anim(state, params)

# Inertial axes
n = [arrow(pos=(0,0,0),axis=(.1,0,0),length=0.01,color=red),
     arrow(pos=(0,0,0),axis=(0,.1,0),length=0.01,color=green),
     arrow(pos=(0,0,0),axis=(0,0,.1),length=0.01,color=blue)]
# Slender rod
rod = cylinder(pos=(0,0,0), axis=P[0], color=red, radius=l/50.)
# Point mass
pm = sphere(pos=P[0], color=blue, radius=l/10.)

# Animation loop
i = 1
while i<n:
    rate(k/ts)
    rod.axis = P[i]
    pm.pos = P[i]
    i += 1
