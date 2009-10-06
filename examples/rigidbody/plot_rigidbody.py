#!/usr/bin/env python
import rigidbody_lib as rb
from scipy.integrate import odeint
from numpy import array, arange, zeros

# Dimensions of rigid body in the three body fixed directions
# Following are the dimensions of an iPhone 3G taken from apple.com
h = 0.1155      # meters  in the 1 direction
w = 0.0621      # meters  in the 2 direction
d = 0.0123      # meters  in the 3 direction
m = 0.135       # kilograms
g = 0.0081        # meters / sec**2
I11 = m*(w**2 + d**2)/12.
I22 = m*(h**2 + d**2)/12.
I33 = m*(h**2 + w**2)/12.

params = [m, 0, I11, I22, I33]

# states = [q1, q2, q3, q4, q5, q6, u1, u2, u3, u4, u5, u6]
# q1, q2, q3 are Body Fixed (Euler) 3-1-2 angles
# q4, q5, q6 are x, y, z Inertial positions
# u1, ..., u6 are the generalized speeds.
# Gravity is in the positive z direction, defined to be downwards

# Specify the initial conditions of the coordinates and the generalized speeds
q0 = [0.0, 0.0, 0.0, .05, 0., 0.]
# Intermediate inertia axis is the body-2 axis, exhibits instability
u0 = [0.0, 2.0, 0.15, 0., 0., 0.0]
x0 = q0 + u0

# Integration time
ti = 0.0
ts = 0.01
tf = 40.0
t = arange(ti, tf+ts, ts)
n = len(t)
# Integrate the differential equations
x = odeint(rb.eoms, x0, t, args = (params,))

# Animate using Visual-Python
AO = zeros((n,3))
A1 = zeros((n,3))
A3 = zeros((n,3))

# Animation playback speed multiplier (1 == realtime)
k = 1.0

for i, state in enumerate(x[:,:6]):
    AO[i], A1[i], A3[i] = rb.anim(state, params)
    A1[i] *= h

from visual import box, display, rate, arrow
black = (0,0,0)
red = (1, 0, 0)
green = (0, 1, 0)
blue = (0, 0, 1)
scene = display(title='Rigid body animation @ %0.2f realtime'%k, width=800, height=800, up=(0,0,-1),\
        uniform=1, background=black, forward=(1,0,0))
N = [arrow(pos=(0,0,0),axis=(.1,0,0),length=0.01,color=red),
     arrow(pos=(0,0,0),axis=(0,.1,0),length=0.01,color=green),
     arrow(pos=(0,0,0),axis=(0,0,.1),length=0.01,color=blue)]

body = box(pos=AO[0], axis=A1[0], up=A3[0],\
        height=d, width=w, color=red)
i = 1
while i<n:
    body.pos = AO[i]
    body.axis = A1[i]
    body.up = A3[i]
    i += 1
    rate(k/ts)
