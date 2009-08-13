from rollingdisc_eoms import *

from scipy.integrate import odeint
from numpy import array, arange, zeros
import matplotlib.pyplot as plt

# Dimensions of disc
m = 1.0
g = 1.0
r = 1.0
params = [m, g, r]

# states = [q1, q2, q3, q4, q5, u1, u2, u3]
# q1, q2, q3 are Body Fixed (Euler) 3-1-2 angles
# q4, q5 are N[1], N[2] Inertial positions
# u1, u2, u3 are the generalized speeds.
# Gravity is in the positive N[3] direction

# Specify the initial conditions of the coordinates and the generalized speeds
q0 = [0., 0., 0., 0., 0.]
u0 = qdot2u(q0, [0.1, -0.6, 3.0], params)
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
C1_axis = zeros((n, 4))
C3_axis = zeros((n, 4))

for i, state in enumerate(x[:,:5]):
    CO_pos[i], C_axis_angle[i], C1_axis[i], C3_axis[i] = animate(state, params)

from visual import ring, display, rate, arrow, curve
scene = display(title='Rolling disc animation', width=800, height=800, up=(0,0,-1),\
        uniform=1, background=(1,1,1), forward=(1,0,0), exit=0,\
        autocenter=True, autoscale=True)
black = (0,0,0)
red = (1, 0, 0)
green = (0, 1, 0)
blue = (0, 0, 1)
n = [arrow(pos=(0,0,0),axis=(.1,0,0),length=0.01,color=red),
     arrow(pos=(0,0,0),axis=(0,.1,0),length=0.01,color=green),
     arrow(pos=(0,0,0),axis=(0,0,.1),length=0.01,color=blue)]
body = ring(pos=CO_pos[0, :], axis=C_axis_angle[0, :3],\
        angle=C_axis_angle[0, 3], radius=r, thickness=0.01, color=(0.5,0.5,0.5))
c1c3 = [arrow(pos=CO_pos[0,:],axis=0.3*C1_axis[0,:3],length=.01,color=red),
            arrow(pos=CO_pos[0,:],axis=0.3*C3_axis[0,:3],length=.01,color=blue)]
body.trail = curve()
body.trail.append(pos=(x[0,3], x[0,4], 0.0), color=black)

i = 1
while i<n:
    rate(1./ts)
    body.pos = CO_pos[i, :]
    body.axis = C_axis_angle[i, :3]
    c1c3[0].pos = body.pos
    c1c3[1].pos = body.pos
    c1c3[0].axis = 0.3*C1_axis[i,:3]
    c1c3[1].axis = 0.3*C3_axis[i,:3]
    body.trail.append(pos=(x[i,3], x[i,4], 0.0), color=black)
    i += 1
