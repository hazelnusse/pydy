#!/usr/bin/env python
from rxnw_eoms import f, animate
from scipy.integrate import odeint
from numpy import arange, zeros, pi, sin
import matplotlib.pyplot as plt

# Define the parameters
I1 = 2575.0
I2 = 0.0625
J1 = 5000.0
J2 = 0.125
l = 1
m1 = 100.0
m2 = 1.0
M = 10.
# Assemble parameters into a list (must be correct order)
params = [m1, m2, I1, I2, J1, J2, l, M]

# Initial Conditions
q0 = [0.0, 0.0, 0.0, pi/4, -pi/4, 0.0, 0.0, 1.0]
u0 = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
x0 = q0 + u0

# Integration time
ti = 0.0
ts = 0.01
tf = 10.0
t = arange(ti, tf+ts, ts)
n = len(t)

# Integrate the differential equations
x = odeint(f, x0, t, args = (params,))
u1 = x[:, 8]
u2 = x[:, 9]
u3 = x[:, 10]
"""
# Plot the results and save a .pdf
plt.figure(1)
plt.subplot(311)
plt.ylabel('$u_1$')
plt.plot(t, u1)
plt.subplot(312)
plt.ylabel('$u_2$')
plt.plot(t, u2)
plt.subplot(313)
plt.ylabel('$u_3$')
plt.plot(t, u3)
plt.xlabel('$t$')
plt.savefig('rxnw_plot.pdf')
plt.figure(2)
plt.plot(t, x[:,2])
"""
# Generate a nice animation
AO = zeros((n, 3))
BO = zeros((n, 3))
CO = zeros((n, 3))
A1 = zeros((n, 4))
A2 = zeros((n, 4))
A3 = zeros((n, 4))
B1 = zeros((n, 4))
B3 = zeros((n, 4))
C1 = zeros((n, 4))
C3 = zeros((n, 4))
# Animation rate multiplier
k = 1.0

# Following parameters are only for visualization purposes, i.e, the don't
# change the dynamics
A_l = 1.0    # Length of body A along A[3] axis
A_wh = 0.5   # Width and height of body A along A[1] and A[2] axes
cyl_l = 0.05     # Length of reaction wheels
r = 1.0         # Radius of reaction wheels
for i, state in enumerate(x[:,:8]):
    AO[i], BO[i], CO[i], A1[i], A2[i], A3[i], B1[i], B3[i], C1[i], C3[i]  = animate(state, params)
    A3[i] *= A_l
    A2[i] *= cyl_l
    B1[i] *= r
    B3[i] *= r
    C1[i] *= r
    C3[i] *= r

from visual import box, display, rate, arrow, cylinder
scene = display(title='Rigid body with reaction wheel @ %0.2f realtime'%k,
        width=800, height=800, uniform=1, background=(1,1,1), up=(0,0,1), forward=(1,0,0))
black = (1,1,1)
red = (1, 0, 0)
green = (0, 1, 0)
blue = (0, 0, 1)
grey = (0.5,0.5,0.5)
N = [arrow(pos=(0,0,0),axis=(.5,0,0),color=red),
     arrow(pos=(0,0,0),axis=(0,.5,0),color=green),
     arrow(pos=(0,0,0),axis=(0,0,.5),color=blue)]
A = box(pos=AO[0, :], axis=A3[0, :3], height=A_wh, width=A_wh, color=red,
        up=A1[0,:3])

B = cylinder(pos=BO[0,:], axis=A2[0, :3], radius=r, color=green)
C = cylinder(pos=CO[0,:], axis=-A2[0, :3], radius=r, color=blue)
B1a = arrow(pos=BO[0,:]+cyl_l*A2[0,:3], axis=B1[0,:3], color=blue)
B3a = arrow(pos=BO[0,:]+cyl_l*A2[0,:3], axis=B3[0,:3], color=blue)
C1a = arrow(pos=CO[0,:]-cyl_l*A2[0,:3], axis=C1[0,:3], color=green)
C3a = arrow(pos=CO[0,:]-cyl_l*A2[0,:3], axis=C3[0,:3], color=green)


i = 1
while i<n:
    rate(k/ts)
    A.pos = AO[i, :]
    A.axis = A3[i, :3]
    A.up = A1[i,:3]
    B.pos = BO[i, :]
    B.axis = A2[i, :3]
    C.pos = CO[i, :]
    C.axis = -A2[i, :3]
    B1a.pos = BO[i,:]+cyl_l*A2[i,:3]
    B1a.axis = B1[i,:3]
    B3a.pos = BO[i,:]+cyl_l*A2[i,:3]
    B3a.axis = B3[i,:3]
    C1a.pos = CO[i,:]-cyl_l*A2[i,:3]
    C1a.axis = C1[i,:3]
    C3a.pos = CO[i,:]-cyl_l*A2[i,:3]
    C3a.axis = C3[i,:3]
    i += 1
