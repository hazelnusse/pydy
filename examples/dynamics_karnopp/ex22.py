from pydy import *

N = NewtonianReferenceFrame('N')
q, qd = N.declare_coords('q', 3)

# Use Python list unpacking to access each individual variable
q1, q2, q3 = q             # r, theta, z
q1d, q2d, q3d = qd         # dr/dt, dtheta/dt, dz/dt

# Create the U reference frame by a rotation about the 3 axis
U = N.rotate('U', 3, q2)

# Locate P relative to N.O
P = N.O.locate('P', q1*U[1] + q3*U[3])

# Caclulate velocity and acceleration of P in N
v_p_n = dt(P.rel(N.O), N)
a_p_n = dt(v_p_n, N)

print 'Position', P.rel(N.O)
print 'Position (expressed in N coordinates)', P.rel(N.O).express(N)
print 'Velocity', v_p_n
print 'Velocity (expressed in N coordinates)', v_p_n.express(N)
print 'Acceleration', a_p_n
print 'Acceleration (expressed in N coordinates)', a_p_n.express(N)
