#!/usr/bin/env python
from pydy import *

# Create a Newtonian reference frame
N = NewtonianReferenceFrame('N')

# Declare parameters, coordinates, speeds
m, g, l, b = N.declare_parameters('m g l b')
(q1,), q_list, qdot_list = N.declare_coords('q', 1)
(u1,), u_list, udot_list = N.declare_speeds('u', 1)

# Orient the A reference frame
A = N.rotate('A', 1, q1)
# Define the position from the hinge
P = N.O.locate('P', -l*A[3], mass=m)

# Define the generalized speed
u_defs = N.define_speeds([Eq(u1, dot(A.ang_vel(), A[1]))])

# Compute the kinematic differential equations
T, Tinv, kindiffs = N.form_kindiffs(u_defs, qdot_list)

print 'Kinematic differential equations'
for qd in qdot_list:
    print qd, '=', kindiffs[qd]

# Set the angular velocity and velocities to only involve expressions involving
# the generalized speeds
A._wrel = Vector(u1*A[1])
P._vrel = cross(A.ang_vel(), P.rel(N.O))

# Set the kinematic differential equations
N.setkindiffs(kindiffs)

# Apply gravity
N.gravity(-g*N[3])
# Apply a friction torque at hinge
A.apply_torque(-b*u1*A[1], N)

# Form Kane's equations and solve them for the udots
kanes_eqns = N.form_kanes_equations()
dyndiffs = N.solve_kanes_equations()
print 'Dynamic differential equations'
for ud in udot_list:
    dyndiffs[ud] = dyndiffs[ud].expand()
    print ud, '=', dyndiffs[ud]

N.setdyndiffs(dyndiffs)
N.output_eoms('pendulum_eoms.py', (A[3], 0))
