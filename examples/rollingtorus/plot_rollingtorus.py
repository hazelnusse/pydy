#!/usr/bin/env python
import rollingtorus_lib as rt
from scipy.integrate import odeint
from numpy import array, arange, zeros, roots
import matplotlib.pyplot as plt

# Dimensions of a penny
m = .100   # A penny has a mass of 3.1g
g = 9.81        # Gravitational acceleration
r1 = .594       # Radius of a penny
r2 = 0.025
I = (r1**2/2 + 5*r2**2/8)*m  # Central moment of inertia about any diameter
J = (r1**2 + 3*r2**2/4)*m    # Central moment of inertia about normal axis
params = [m, g, r1, r2, I, J]

def plot_energy(t, x):
    # Plot the kinetic energy, potential energy, and total energy
    ke = zeros((n,1))
    pe = zeros((n,1))
    te = zeros((n,1))
    for i in range(n):
        ke[i], pe[i] = rt.energy(x[i,:], params)
        te[i] = ke[i] + pe[i]

    plt.figure()
    plt.plot(t, ke, label='KE')
    plt.plot(t, pe, label='PE')
    plt.plot(t, te, label='TE')
    plt.legend()
    plt.title(('Energy of rolling torus\nm = %0.2f, g = %0.2f, r1 = %0.2f, '+\
              'r2 = %0.2f, I = %0.2f, J = %0.2f')%(m, g, r1, r2, I, J))
    plt.xlabel('Time [s]')
    plt.ylabel('Energy [kg * m ** 2 / s**2]')
    plt.show()
    ############################


# states = [q1, q2, q3, q4, q5, u1, u2, u3]
# q1, q2, q3 are Body Fixed (Euler) 3-1-2 angles
# q1:  Yaw / heading angle
# q2:  Lean angle
# q3:  Spin angle
# q4, q5 are N[1], N[2] Inertial positions
# u1, u2, u3 are the generalized speeds.
# Gravity is in the positive N[3] direction

# Specify the initial conditions of the coordinates and the generalized speeds
qi = [0., -0.3, 0.0, .001, 0.001]

# Alternatively, specify any other intial generalized speeds
ui = [-.25,2.0,1.0]

# Inital states
xi = qi + ui

# Integration time
ti = 0.0
ts = 0.01
tf = 40.0
t = arange(ti, tf+ts, ts)
n = len(t)
# Integrate the differential equations
x = odeint(rt.eoms, xi, t, args = (params,))

# Plot energy
plot_energy(t, x)

# Animate using Visual-Python
CO = zeros((n,3))
B2 = zeros((n,3))
C1 = zeros((n, 3))
C3 = zeros((n, 3))

# Animation playback speed multiplier (1 == realtime)
k = 0.5

for i, state in enumerate(x[:,:5]):
    CO[i], B2[i], C1[i], C3[i] = rt.anim(state, params)
    # Make the out of plane axis shorter since this is what control the height
    # of the cone
    C1[i] *= r1
    C3[i] *= r1

from visual import display, rate, arrow, curve, cone, box, ring
import visual
black = (0,0,0)
red = (1, 0, 0)
green = (0, 1, 0)
blue = (0, 0, 1)
white = (1, 1, 1)
NO = (0,0,0)
N1 = (1, 0, 0)
N2 = (0, 1, 0)
N3 = (0, 0, 1)

scene = display(title='Rolling torus @ %0.2f realtime'%k, up=(0,0,-1),\
        uniform=1, background=black, forward=N1)

# Inertial reference frame arrows
N = [arrow(pos=NO,axis=N1,color=red),
     arrow(pos=NO,axis=N2,color=green),
     arrow(pos=NO,axis=N3,color=blue)]

torus = ring(pos=CO[0], axis=B2[0], radius=r1, thickness=r2, color=blue)
# Arrows for body fixed coordinates in plane of torus
c1 = arrow(pos=CO[0], axis=C1[0], up=C3[0], color=red)
c3 = arrow(pos=CO[0], axis=C3[0], up=C1[0], color=green)
# Ground contact path
trail = curve()
trail.append(pos=(x[0,3], x[0,4], 0.0), color=white)

i = 1
while i<n:
    torus.pos = CO[i]
    torus.axis = B2[i]
    c1.pos = CO[i]
    c3.pos = CO[i]
    c1.axis = C1[i]
    c3.axis = C3[i]
    c1.up = C3[i]
    c3.up = C1[i]
    trail.append(pos=(x[i,3], x[i,4], 0.0), color=white)
    i += 1
    rate(k/ts)
