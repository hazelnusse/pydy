#!/usr/bin/env python
from sympy import simplify, solve
from pydy import *

# Create a Newtonian reference frame
N = NewtonianReferenceFrame('N')

# Declare parameters, coordinates, speeds
params = m1, m2, g, l1, l2, b1, b2 = N.declare_parameters('m1, m2 g l1 l2 b1 b2')
(q1, q2), (q1d, q2d) = q, qd = N.declare_coords('q', 2)
(u1, u2), (u1d, u2d) = u, ud = N.declare_speeds('u', 2)

# Orient the A reference frame
A = N.rotate('A', 1, q1)
B = A.rotate('B', 1, q2)
# Define the position from the hinge
P1 = N.O.locate('P1', -l1*A[3], mass=m1)
P2 = P1.locate('P2', -l2*B[3], mass=m2)

# Define the generalized speed
u_defs = [Eq(u1, dot(A.ang_vel(), A[1])),
          Eq(u2, dot(B.ang_vel(), A[1]))]

kd = solve(u_defs, qd)
kd_eqs = [Eq(lhs, rhs) for lhs, rhs in zip(qd, [kd[q1d], kd[q2d]])]

# Set the angular velocity and velocities to only involve expressions involving
# the generalized speeds
A.abs_ang_vel = Vector(u1*A[1])
B.abs_ang_vel = Vector(u2*A[1])
P1.abs_vel = cross(A.ang_vel(N), P1.rel(N.O))
P2.abs_vel = P1.abs_vel + cross(B.ang_vel(N), P2.rel(P1))

A.abs_ang_acc = dt(A.abs_ang_vel, N)
B.abs_ang_acc = dt(B.abs_ang_vel, N)
P1.abs_acc = dt(P1.abs_vel, N)
P2.abs_acc = dt(P2.abs_vel, N)

# Apply gravity
N.gravity(-g*N[3])
# Apply a friction torque at hinges
A.apply_torque(-b1*u1*A[1], N)
B.apply_torque(-b2*(u2-u1)*A[1], A)

# Form Kane's equations and solve them for the udots
kanes_eqns = N.form_kanes_equations()
dd, subs_dict = N.solve_kanes_equations(dummy_vars=True)

# Energy calculations
ke = (m1*P1.abs_vel.mag_sqr + m2*P2.abs_vel.mag_sqr)/S(2)
pe = g*(m1*P1.rel(N.O).dot(N[3]) + m2*P2.rel(N.O).dot(N[3]))
energy_eqs = [Eq(Symbol('ke'), ke), Eq(Symbol('pe'), pe)]

# Animation equations
anim_eqs = animate(N, ('P1', P1.rel(N.O)), ('P2', P2.rel(N.O)), ('B3', B[3]))

print "Kinematic differential equations"
for e in kd_eqs:
    print e
print "Mass matrix and right hand side terms of dynamic equations"
for k, v in subs_dict.items():
    print k, "==", v
print 'Dynamic differential equations'
for i in (0,1):
    print dd[i]

output_string = "from __future__ import division\n"
output_string += "from math import sin, cos\n\n"

ds = """\
Double Pendulum equations of motion.

    _x is an array/list in the following order:
        q1:  Angle of first pendulum link relative to vertical (0 downwards)
        q2:  Angle of second pendulum link relative to first (0 is fully extended)
        u1:  A[1] measure number of the inertial angular velocity of the first link.
        u2:  A[1] measure number of the inertial angular velocity of the second link.

    _params is an array/list in the following order:
        m1:  Mass of first pendulum point mass.
        m2:  Mass of second pendulum point mass.
        l1:  Length of first pendulum link.
        l2:  Length of second pendulum link.
         g:  Gravitational constant.
"""
output_string += generate_function("eoms", kd_eqs+dd, q+u, params,
        nested_terms=[subs_dict], docstring=ds, time=True)

ds = """\
Kinetic and Potential Energy of double pendulum.

    _x is an array/list in the following order:
        q1:  Angle of first pendulum link relative to vertical (0 downwards)
        q2:  Angle of second pendulum link relative to first (0 is fully extended)
        u1:  A[1] measure number of the inertial angular velocity of the first link.
        u2:  A[1] measure number of the inertial angular velocity of the second link.

    _params is an array/list in the following order:
        m1:  Mass of first pendulum point mass.
        m2:  Mass of second pendulum point mass.
        l1:  Length of first pendulum link.
        l2:  Length of second pendulum link.
         g:  Gravitational constant.
    Returns a list/array of kinetic energy and potential energy, respectively.
"""
output_string += generate_function("energy", energy_eqs, q+u, params, docstring=ds)

ds = """\
Calculate configuration of double pendulum for purposes of animation.

    _x is an array/list of the configuration coordinates of the disc:
        q1:  Angle of first pendulum link relative to vertical (0 downwards)
        q2:  Angle of second pendulum link relative to first (0 is fully extended)
        u1:  A[1] measure number of the inertial angular velocity of the first link.
        u2:  A[1] measure number of the inertial angular velocity of the second link.

    _params is the radius of the disc.
        l1:  Length of first pendulum link.
        l2:  Length of second pendulum link.

    Output is four 3-tuples in the following order:
          P1:  Position of first point mass.
          P2:  Position of second point mass.
        B[3]:  Axis along second link of pendulum.
"""
output_string += generate_function("anim", anim_eqs, q, params, triples=True,
        docstring=ds)
file = open('doublependulum_lib.py', 'w')
file.write(output_string)
file.close()
