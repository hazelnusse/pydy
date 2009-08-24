from rollingtorus_eoms import *
from scipy.integrate import odeint
from numpy import array, arange, zeros, roots


# Dimensions of a penny
m = .100   # A penny has a mass of 3.1g
g = 9.81        # Gravitational acceleration
r1 = .594       # Radius of a penny
r2 = 0.025
I = (r1**2/2 + 5*r2**2/8)*m  # Central moment of inertia about any diameter
J = (r1**2 + 3*r2**2/4)*m    # Central moment of inertia about normal axis
params = [m, g, r1, r2, I, J]

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
x = odeint(f, xi, t, args = (params,))

# Animate using Visual-Python
CO_pos = zeros((n,3))
C_axis_angle = zeros((n,4))
C1_axis = zeros((n, 4))
C3_axis = zeros((n, 4))

# Animation playback speed multiplier (1 == realtime)
k = 0.5

for i, state in enumerate(x[:,:5]):
    CO_pos[i], C_axis_angle[i], C1_axis[i], C3_axis[i] = animate(state, params)
    # Make the out of plane axis shorter since this is what control the height
    # of the cone
    C_axis_angle[i] *= 0.001
    C1_axis[i] *= r1
    C3_axis[i] *= r1

from visual import display, rate, arrow, curve, cone, box, ring
import visual
black = (0,0,0)
red = (1, 0, 0)
green = (0, 1, 0)
blue = (0, 0, 1)
white = (1, 1, 1)
NO = (0,0,0)
#scene = display(title='Rolling torus @ %0.2f realtime'%k, up=(0,0,-1),\
#        uniform=1, background=black, forward=(1,0,0), exit=0)

scene = display(title='Rigid body animation @ %0.2f realtime'%k, width=800,
        height=800, up=(0,0,-1), uniform=1, background=(0,0,0), forward=(1,0,0), exit=0)
# Inertial reference frame arrows
N = [arrow(pos=NO,axis=(.001,0,0),color=red),
     arrow(pos=NO,axis=(0,.001,0),color=green),
     arrow(pos=NO,axis=(0,0,.001),color=blue)]
# Two cones are used to look like a thin disc
torus = ring(pos=CO_pos[0, :], axis=C_axis_angle[0, :3],\
        angle=C_axis_angle[0, 3], radius=r1, thickness=r2, color=blue)
# Body fixed coordinates in plane of disc, can't really be seen through cones
c1c3 = [arrow(pos=CO_pos[0,:],axis=C1_axis[0,:3],length=.01,color=red),
            arrow(pos=CO_pos[0,:],axis=C3_axis[0,:3],length=.01,color=green)]
trail = curve()
trail.append(pos=(x[0,3], x[0,4], 0.0), color=white)

i = 1
while i<n:
    torus.pos = CO_pos[i, :]
    torus.axis = C_axis_angle[i, :3]
    c1c3[0].pos = CO_pos[i, :]
    c1c3[1].pos = CO_pos[i, :]
    c1c3[0].axis = C1_axis[i,:3]
    c1c3[1].axis = C3_axis[i,:3]
    trail.append(pos=(x[i,3], x[i,4], 0.0), color=white)
    i += 1
    rate(k/ts)
