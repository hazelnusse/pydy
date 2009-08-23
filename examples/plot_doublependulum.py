from doublependulum_eoms import f, animate
from scipy.integrate import odeint
from numpy import array, arange, zeros, pi

m1 = .1
m2 = .1
g = 9.8
l1 = 0.1
l2 = 0.1
b1 = .1
b2 = 0.001
params = [m1, m2, g, l1, l2, b1, b2]
x0 = [pi/4, 0., 0, -0.1]

# Integration time
ti = 0.0
ts = 0.01
tf = 40.0
t = arange(ti, tf+ts, ts)
n = len(t)
# Integrate the differential equations
x = odeint(f, x0, t, args = (params,))

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
P1 = zeros((n,4))
P2 = zeros((n,4))

# Call the animate function to get the A[3] unit vector
for i, state in enumerate(x[:,:2]):
    P1[i], P2[i] = animate(state, params)
    P1[i] *= -l1
    P2[i] *= -l2
# Inertial axes
N = [arrow(pos=(0,0,0),axis=(.1,0,0),length=0.01,color=red),
     arrow(pos=(0,0,0),axis=(0,.1,0),length=0.01,color=green),
     arrow(pos=(0,0,0),axis=(0,0,.1),length=0.01,color=blue)]

# Slender rod
rod1 = cylinder(pos=(0,0,0), axis=P1[0, :3], color=red, radius=l1/50.)
rod2 = cylinder(pos=P1[0, :3], axis=P2[0, :3], color=red, radius=l2/50.)
# Point mass
pm1 = sphere(pos=P1[0, :3], color=blue, radius=m1/10.)
pm2 = sphere(pos=P1[0, :3]+P2[0, :3], color=blue, radius=m2/10.)


# Animation loop
i = 1
while i<n:
    rod1.axis = P1[i,:3]
    rod2.pos = P1[i,:3]
    rod2.axis = P2[i, :3]
    pm1.pos = P1[i, :3]
    pm2.pos = P1[i, :3] + P2[i, :3]
    i += 1
    rate(k/ts)
