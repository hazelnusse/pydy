from sympy import solve
from pydy import *

# Create a Newtonian reference frame
N = NewtonianReferenceFrame('N')

# Declare parameters, coordinates, speeds
m, g, I11, I22, I33 = N.declare_parameters('m g I11 I22 I33')
(q1, q2, q3, q4, q5, q6), q_list, qdot_list = N.declare_coords('q', 6)
(u1, u2, u3, u4, u5, u6), u_list, udot_list = N.declare_speeds('u', 6)

# Frame fixed to the rigid body
C = N.rotate("C", 'BODY312', (q1, q2, q3), I=(I11, I22, I33, 0, 0, 0))

# Locate the mass center
CO = N.O.locate('CO', q4*N[1] + q5*N[2] + q6*N[3], mass=m)

# Define the generalized speeds
u_defs = N.define_speeds(
        [Eq(u_list[i-1], dot(C.ang_vel(), C[i])) for i in (1, 2, 3)] + \
        [Eq(u_list[i-1], dot(CO.vel(), N[i-3])) for i in (4, 5, 6)])

# Form transformation matrix of u := T * q', return T, inv(T), and the kindiffs
T, Tinv, kindiffs = N.form_kindiffs(u_defs, qdot_list, method='ADJ')

print 'Kinematic differential equations'
for qd in qdot_list:
    print qd, '=', kindiffs[qd]

# Set the angular velocity and velocity to only involve generalized speeds,
# must be done *exactly* as they were defined when .define_speeds was called in
# order for the correct equations of motion to result.
C._wrel = Vector(u1*C[1] + u2*C[2] + u3*C[3])
CO._vrel = Vector(u4*N[1] + u5*N[2] + u6*N[3])

# Sets the kinematic differential equations
N.setkindiffs(kindiffs)

# Apply gravity
N.gravity(g*N[3])

# Form Kane's equations and solve them for the udots
kanes_eqns = N.form_kanes_equations()

dyndiffs = solve(kanes_eqns, udot_list)
print 'Dynamic differential equations'
for r, ud in enumerate(udot_list):
    print ud, ' = ', dyndiffs[ud]

N.setdyndiffs(dyndiffs)
N.output_eoms('rigidbody_eoms.py', (CO, N.O), (C[2], q3))
