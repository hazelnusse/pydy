from rigidbody_eoms import *

from scipy.integrate import odeint
from numpy import array, arange, zeros
import matplotlib.pyplot as plt

# Dimensions of rigid body in the three body fixed directions
# Following are the dimensions of an iPhone 3G taken from apple.com
h = 0.1155      # meters  in the 1 direction
w = 0.0621      # meters  in the 2 direction
d = 0.0123      # meters  in the 3 direction
m = 0.135       # kilograms
g = 9.81        # meters / sec**2
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
q0 = [0.1, 0.1, 0.1, 0., 0., 0.]
u0 = [0.0, 0.0, 1.0, 0., 0., 0.]
x0 = q0 + u0

# Integration time
ti = 0.0
ts = 0.01
tf = 40.0
t = arange(ti, tf+ts, ts)
n = len(t)
# Integrate the differential equations
x = odeint(f, x0, t, args = (params,))

# Animate using Visual-Python
CO_pos = zeros((n,3))
C_axis_angle = zeros((n,4))

# Animation playback speed multiplier (1 == realtime)
k = 0.5

for i, state in enumerate(x[:,:6]):
    CO_pos[i], C_axis_angle[i] = animate(state, params)
    C_axis_angle[i] *= w

from visual import box, display, rate, arrow
scene = display(title='Rigid body animation @ %0.2f realtime', width=800, height=800, up=(0,0,-1),\
        uniform=1, background=(1,1,1), forward=(1,0,0), exit=0)
black = (1,1,1)
red = (1, 0, 0)
green = (0, 1, 0)
blue = (0, 0, 1)
n = [arrow(pos=(0,0,0),axis=(.1,0,0),length=0.01,color=red),
     arrow(pos=(0,0,0),axis=(0,.1,0),length=0.01,color=green),
     arrow(pos=(0,0,0),axis=(0,0,.1),length=0.01,color=blue)]
body = box(pos=CO_pos[0, :], axis=C_axis_angle[0, :3], angle=C_axis_angle[0, 3],\
        length=w, height=h, width=d, color=(0.5,0.5,0.5))
i = 1
while i<n:
    rate(k/ts)
    body.pos = CO_pos[i, :]
    body.axis = C_axis_angle[i, :3]
    body.angle = C_axis_angle[i, 3]
    i += 1
