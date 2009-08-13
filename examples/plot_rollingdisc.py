from rollingdisc_eoms import *

from scipy.integrate import odeint
from numpy import array, arange, zeros, roots

# Dimensions of a penny
m = 3.1/1000.   # A penny has a mass of 3.1g
g = 9.81        # Gravitational acceleration
r = 0.019/2.    # Radius of a penny
# Dimensions of a quarter 
m = 5.67/1000.  # A quarter has a mass of 5.67g
r = 0.02426/2.    # Radius of a penny
params = [m, g, r]

# states = [q1, q2, q3, q4, q5, u1, u2, u3]
# q1, q2, q3 are Body Fixed (Euler) 3-1-2 angles
# q1:  Yaw / heading angle
# q2:  Lean angle
# q3:  Spin angle
# q4, q5 are N[1], N[2] Inertial positions
# u1, u2, u3 are the generalized speeds.
# Gravity is in the positive N[3] direction

# Specify the initial conditions of the coordinates and the generalized speeds
qi = [0., 1.1, 0., .001, 0.001]
# Steady turning conditions require that u1 = 0 and that
# u3**2 - 2*c2*u2/s2*u3 - 4*g*c2/(5*r) = 0
# Given lean angle q2 and gen. speed u2, u3 must be a root of the above
# polynomial
q2i = qi[1]
u2i = .0
coefs = [1., -2*u2i/tan(q2i), -4*g*cos(q2i)/(5*r)]
u3i = roots(coefs)
print u3i
ui = [0.0,u2i,u3i[0]]
xi = qi + ui

# Integration time
ti = 0.0
ts = 0.001
tf = 40.0
t = arange(ti, tf+ts, ts)
n = len(t)
# Integrate the differential equations
x = odeint(f, xi, t, args = (params,))

# Animate using Visual-Python
CO_pos = zeros((n,3))
C_axis_angle = zeros((n,4))
C1_axis = zeros((n, 4))
C3_axis = zeros((n, 4))

# Animation playback speed multiplier (1 == realtime)
k = 0.2

for i, state in enumerate(x[:,:5]):
    CO_pos[i], C_axis_angle[i], C1_axis[i], C3_axis[i] = animate(state, params)
    # Make the out of plane axis shorter since this is what control the height
    # of the cone
    C_axis_angle[i] *= 0.001
    C1_axis[i] *= r
    C3_axis[i] *= r

from visual import display, rate, arrow, curve, cone, box
black = (0,0,0)
red = (1, 0, 0)
green = (0, 1, 0)
blue = (0, 0, 1)
white = (1, 1, 1)
NO = (0,0,0)
scene = display(title='Rolling disc @ %0.2f realtime'%k, up=(0,0,-1),\
        uniform=1, background=black, forward=(1,0,0), exit=0,\
        autocenter=True, autoscale=True)
# Inertial reference frame arrows
n = [arrow(pos=NO,axis=(.001,0,0),color=red),
     arrow(pos=NO,axis=(0,.001,0),color=green),
     arrow(pos=NO,axis=(0,0,.001),color=blue)]
# Two cones are used to look like a thin disc
body1 = cone(pos=CO_pos[0, :], axis=C_axis_angle[0, :3],\
        angle=C_axis_angle[0, 3], radius=r, color=blue)
body2 = cone(pos=CO_pos[0, :], axis=-C_axis_angle[0, :3],\
        angle=C_axis_angle[0, 3], radius=r, color=blue)
# Body fixed coordinates in plane of disc, can't really be seen through cones
c1c3 = [arrow(pos=CO_pos[0,:],axis=C1_axis[0,:3],length=.01,color=red),
            arrow(pos=CO_pos[0,:],axis=C3_axis[0,:3],length=.01,color=green)]
trail = curve()
trail.append(pos=(x[0,3], x[0,4], 0.0), color=white)

i = 1
while i<n:
    rate(k/ts)
    body1.pos = CO_pos[i, :]
    body1.axis = C_axis_angle[i, :3]
    body2.pos = CO_pos[i, :]
    body2.axis = -C_axis_angle[i, :3]
    c1c3[0].pos = body1.pos
    c1c3[1].pos = body1.pos
    c1c3[0].axis = C1_axis[i,:3]
    c1c3[1].axis = C3_axis[i,:3]
    trail.append(pos=(x[i,3], x[i,4], 0.0), color=white)
    i += 1
