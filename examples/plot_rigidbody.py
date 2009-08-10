from rigidbody_eoms import *

from scipy.integrate import odeint
from numpy import array, arange
import matplotlib.pyplot as plt

# params = [m, g, I11, I22, I33]
params = array([1., 1., 10., 20., 30.])
# states = [q1, q2, q3, q4, q5, q6, u1, u2, u3, u4, u5, u6]
# q1, q2, q3 are Body Fixed (Euler) 3-1-2 angles
# q4, q5, q6 are x, y, z Inertial positions
# u1, ..., u6 are the generalized speeds.
# Gravity is in the positive z direction, defined to be downwards

# Specify the initial conditions of the coordinates and coordinate rates
q0 = [0., 0., 0., 0., 0., 0.]
qdot0 = [.00, 1.0, 0.001, 0., 0., 0.]

# Form the initial generalized speeds based upon the inital coordinates and
# their rates
u0 = qdot2u(q0, qdot0, params)
x0 = array(q0 + u0)

# Integration time 
ti = 0.0
ts = 0.01
tf = 40.0
t = arange(ti, tf+ts, ts)

# Integrate the differential equations
x = odeint(f, x0, t, args = (params,))

# Plot some results
plt.plot(t, x[:,0:3])
plt.legend(('q1 (3-axis)', 'q2 (1-axis)', 'q3 (2-axis)'))
plt.show()

# Animate using Visual-Python

