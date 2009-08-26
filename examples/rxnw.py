#!/usr/bin/env python
from pydy import *
from sympy import solve

# Declare a Newtonian Reference Frame
N = NewtonianReferenceFrame('N')

# Declare parameters, coordinates, and speeds
m1, m2, I1, I2, J1, J2, l, M = N.declare_parameters('m1 m2 I1 I2 J1 J2 l M')
(q1, q2, q3, q4, q5, q6, q7, q8), q_list, qdot_list = N.declare_coords('q', 8)
(u1, u2, u3, u4, u5, u6, u7, u8), u_list, udot_list = N.declare_speeds('u', 8)

# Orient the three rigid bodies, assign inertias
A = N.rotate('A', 'BODY123', (q1, q2, q3), I=(I1, I1, J1, 0, 0, 0))
B = A.rotate('B', 2, q4, I=(I2, J2, I2, 0, 0, 0), I_frame=A)
C = A.rotate('C', 2, q5, I=(I2, J2, I2, 0, 0, 0), I_frame=A)

# Locate the mass centers, assign masses
AO = N.O.locate('AO', q6*N[1] + q7*N[2] + q8*N[3], mass=m1)
BO = AO.locate('BO', -l*A[2], mass=m2)
CO = AO.locate('CO', l*A[2], mass=m2)

# Define the generalized speeds
u_defs = N.define_speeds(
        [Eq(u_list[i-1], dot(A.ang_vel(), A[i])) for i in (1, 2, 3)] +\
        [Eq(u_list[3], dot(B.ang_vel(A), A[2])),\
         Eq(u_list[4], dot(C.ang_vel(A), A[2]))] +\
        [Eq(u_list[i+4], dot(AO.vel(), A[i])) for i in (1, 2, 3)])

# Form the kinematic differential equations and the transformation matrices
# between the qdot's and u's
T, Tinv, kindiffs = N.form_kindiffs(u_defs, qdot_list)
print "Kinematic differential equations"
for qd in qdot_list:
    print qd, ' = ', kindiffs[qd]

# Manually setting all velocities to only involve speeds
# Need to improve the API here, but this works.
A._wrel = Vector(u1*A[1] + u2*A[2] + u3*A[3])
B._wrel = Vector(u4*A[2])
C._wrel = Vector(u5*A[2])
AO._vrel = Vector(u6*A[1] + u7*A[2] + u8*A[3])
BO._vrel = cross(B.ang_vel(N), -l*A[2])
CO._vrel = cross(C.ang_vel(N), l*A[2])
N.setkindiffs(kindiffs)

# Apply the torque between each rigid body
A.apply_torque(-M*A[2], B)
A.apply_torque(-M*A[2], C)

# Form Kane's equations
kanes_eqns = N.form_kanes_equations()
# Solve Kane's equations for the udot's
dyn_diffs = N.solve_kanes_equations()
print "Dynamic differential equations:"
for u in udot_list:
    print u, '=', dyn_diffs[u]

N.setdyndiffs(dyn_diffs)

N.output_eoms('rxnw_eoms.py', (AO, N.O), (BO, N.O), (CO, N.O), (A[1], 0),
        (A[2], 0), (A[3], 0), (B[1], 0), (B[3], 0), (C[1], 0), (C[3], 0))
