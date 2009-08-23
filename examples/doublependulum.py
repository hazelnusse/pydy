from sympy import simplify
from pydy import *

# Create a Newtonian reference frame
N = NewtonianReferenceFrame('N')

# Declare parameters, coordinates, speeds
m1, m2, g, l1, l2, b1, b2 = N.declare_parameters('m1, m2 g l1 l2 b1 b2')
(q1, q2), q_list, qdot_list = N.declare_coords('q', 2)
(u1, u2), u_list, udot_list = N.declare_speeds('u', 2)

# Orient the A reference frame
A = N.rotate('A', 1, q1)
B = A.rotate('B', 1, q2)
# Define the position from the hinge
P1 = N.O.locate('P1', -l1*A[3], mass=m1)
P2 = P1.locate('P2', -l2*B[3], mass=m2)

# Define the generalized speed
u_defs = N.define_speeds(
        [Eq(u1, dot(A.ang_vel(), A[1])),
         Eq(u2, dot(B.ang_vel(), A[1]))])

# Compute the kinematic differential equations
T, Tinv, kindiffs = N.form_kindiffs(u_defs, qdot_list)

print 'Kinematic differential equations'
for u_def, qd in zip(u_defs, qdot_list):
    print qd, '=', kindiffs[qd]

# Set the angular velocity and velocities to only involve expressions involving
# the generalized speeds
A._wrel = Vector(u1*A[1])
B._wrel = Vector((u2-u1)*A[1])
P1._vrel = cross(A.ang_vel(), P1.rel(N.O))
P2._vrel = cross(B.ang_vel(), P2.rel(P1))

# Set the kinematic differential equations
N.setkindiffs(kindiffs)

# Apply gravity
N.gravity(-g*N[3])
# Apply a friction torque at hinges
A.apply_torque(-b1*u1*A[1], N)
B.apply_torque(-b2*(u2-u1)*A[1], A)

# Form Kane's equations and solve them for the udots
kanes_eqns = N.form_kanes_equations()

dyndiffs = N.solve_kanes_equations()
print 'Dynamic differential equations'
for ud in udot_list:
    print ud, '=', dyndiffs[ud]

N.setdyndiffs(dyndiffs)
N.output_eoms('doublependulum_eoms.py', (A[3], 0), (B[3], 0))
