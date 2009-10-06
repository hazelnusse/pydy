#!/usr/bin/env python
import doublependulum_lib as dp
from scipy.integrate import odeint
from numpy import array, arange, zeros, pi

m1 = .8
m2 = 1.1
g = 9.8
l1 = 1.1
l2 = 2.1
b1 = 0.5
b2 = 0.5
params = [m1, m2, g, l1, l2, b1, b2]
x0 = [0, pi, 0, -0.1]

# Integration time
ti = 0.0
ts = 0.01
tf = 40.0
t = arange(ti, tf+ts, ts)
n = len(t)
# Integrate the differential equations
x = odeint(dp.eoms, x0, t, args = (params,))
kepe = zeros((n,2))
te = zeros((n,1))
for i in range(n):
    kepe[i] = dp.energy(x[i], params)
    te[i] = kepe[i,0] + kepe[i,1]

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
P1 = zeros((n,3))
P2 = zeros((n,3))
B3 = zeros((n,3))

# Call the animate function to get the A[3] unit vector
for i, state in enumerate(x[:,:2]):
    P1[i], P2[i], B3[i] = dp.anim(state, params)
    B3[i] *= -l2

# Inertial axes
N = [arrow(pos=(0,0,0),axis=(.1,0,0),length=0.01,color=red),
     arrow(pos=(0,0,0),axis=(0,.1,0),length=0.01,color=green),
     arrow(pos=(0,0,0),axis=(0,0,.1),length=0.01,color=blue)]

# Slender rod
rod1 = cylinder(pos=(0,0,0), axis=P1[0], color=red, radius=l1/50.)
rod2 = cylinder(pos=P1[0], axis=P2[0], color=red, radius=l2/50.)
# Point mass
pm1 = sphere(pos=P1[0], color=blue, radius=m1/10.)
pm2 = sphere(pos=P2[0], color=blue, radius=m2/10.)


# Animation loop
i = 1
while i<n:
    rate(k/ts)
    rod1.axis = P1[i]
    rod2.pos = P1[i]
    rod2.axis = B3[i]
    pm1.pos = P1[i]
    pm2.pos = P2[i]
    i += 1
