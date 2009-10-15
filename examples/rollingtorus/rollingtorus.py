#!/usr/bin/env python
from sympy import solve, simplify
from pydy import *

# Create a Newtonian reference frame
N = NewtonianReferenceFrame('N')

# Constants
m, g, r1, r2, I, J = params = N.declare_parameters("m g r1 r2 I J")

# Assume uniform density torus
I = m*r1**2/2 + 5*m*r2**2/8
J = m*r1**2 + 3*m*r2**2/4

# Declare generalized coordinates and generalized speeds
q, qd = N.declare_coords('q', 5)
u, ud = N.declare_speeds('u', 3)
q1, q2, q3, q4, q5 = q
q1d, q2d, q3d, q4d, q5d = qd
u1, u2, u3 = u
u1d, u2d, u3d = ud

# Dictionary used for substituting sin(q1)/cos(q1) with tan(q1)
tan_lean = {sin(q2)/cos(q2): tan(q2)}

# Intermediate reference frames
A = N.rotate("A", 3, q1)
B = A.rotate("B", 1, q2)

# Frame fixed to the torus
C = B.rotate("C", 2, q3, I=(I, J, I, 0, 0, 0), I_frame=B)

# Locate the mass center of torus
CO = N.O.locate('CO', -r2*A[3] - r1*B[3], frame=C, mass=m)

# Fixed inertial reference point
N1 = CO.locate('N1', r1*B[3] + r2*A[3] - q4*N[1] - q5*N[2])

# Define the generalized speeds to be the B frame measure numbers of the angular
# Must be of the form:  A*q' == u
u_rhs = [dot(C.ang_vel(), B[i]) for i in (1, 2, 3)]

# Form the list of equations mapping qdots to generalized speeds
qd_to_u_eqs = [Eq(ui, ui_rhs) for ui, ui_rhs in zip(u, u_rhs)]

# Form the matrix that maps qdot's to u's
qd_to_u = coefficient_matrix(u_rhs, qd[:-2])

# Invert the matrix
u_to_qd = qd_to_u.inv()

# Form the right hand sides of the kindiffs
qd_rhs = u_to_qd * Matrix(u)
# Create a list of Equations:
kd = []
for qdot, eqn in zip(qd[:-2], qd_rhs):
    kd.append(Eq(qdot, eqn))

kd[2] = Eq(kd[2].lhs, kd[2].rhs.subs(tan_lean))

kd_dict = eqn_list_to_dict(kd)

# Form the velocity of the rear wheel center in two different ways in order to
# bring the non holonomic constraints into evidence.
vcon1 = dt(CO.rel(N1), N)
vcon2 = cross(C.ang_vel(N), CO.rel(N.O))
eq1 = dot(vcon1 - vcon2, N[1])
eq2 = dot(vcon1 - vcon2, N[2])

# eq1 and eq2 are linear in the time derivatives of the coordinates, determine
# the matrix of coefficients here.
T = coefficient_matrix([eq1, eq2], qd)
Td_inv, Ti = transform_matrix(T, qd, qd[3:])

contact_rates = Td_inv*Ti*Matrix(qd[:3])

# Append the contact point kinematic differential equations to the list,
# keeping them in implicit form.
kd.append(Eq(q4d, contact_rates[0]))
kd.append(Eq(q5d, contact_rates[1]))
kd_dict = eqn_list_to_dict(kd)

# Set velocities and angular velocities using only generalized speeds
B.abs_ang_vel = Vector(u1*B[1] + u3*B[3])
C.abs_ang_vel = Vector(u1*B[1] + u2*B[2] + u3*B[3])
CO.abs_vel = cross(C.abs_ang_vel, CO.rel(N.O))

# Set accelerations and angular accelerations
C.abs_ang_acc = dt(C.abs_ang_vel, N)
CO.abs_acc = dt(CO.abs_vel, N)

# Apply gravity
N.gravity(g*A[3])

# Form Kane's equations and solve them for the udots
kanes_eqns = N.form_kanes_equations(subs_dict=kd_dict)
N.mass_matrix = N.mass_matrix.expand()

# Solve for the du/dt terms
dyndiffs, dyn_dict = N.solve_kanes_equations(dummy_vars=True)

# Energy calculations
ke = ((m*CO.abs_vel.mag_sqr + dot(C.abs_ang_vel, dot(C.inertia, C.abs_ang_vel)))).expand()/2.
# Take zero potential to be when the disc is upright
pe = -m*g*CO.rel(N.O).dot(N[3]) - m*g*(r1+r2)
energy_eqs= [Eq(Symbol('ke'), ke),
             Eq(Symbol('pe'), pe)]

# Begin code for generation of output functions
output_string = "from __future__ import division\n"
output_string += "from numpy import sin, cos, tan\n\n"
ds = """\
Equations of motion for rolling torus.

    _x is an array/list with the following ordering:
        q1:  Yaw angle (0 is aligned with N[1])
        q2:  Lean angle (0 is upright)
        q3:  Spin angle
        q4:  N[1] measure number of contact point position relative to origin
        q5:  N[2] measure number of contact point position relative to origin
        u1:  B[1] measure number of angular velocity of C in N
        u2:  B[2] measure number of angular velocity of C in N
        u3:  B[3] measure number of angular velocity of C in N

    _params is a array/list with the following ordering:
         m:  Mass of disc
         g:  Gravitational constant
         r1: Major radius of torus
         r2: Minor radius of torus
"""

output_string += generate_function("eoms", kd+dyndiffs, q+u,\
        params, docstring=ds, time=True, nested_terms=[dyn_dict])

# Animation equations
anim_eqs = animate(N, ('CO', CO.rel(N1)), ('B2', B[2]), ('C1', C[1]),\
        ('C3', C[3]))
ds = """\
Calculate position and orientation of torus for purposes of animation.

    _x is an array/list of the configuration coordinates of the disc:
        q1:  Yaw angle (0 is aligned with N[1])
        q2:  Lean angle (0 is upright)
        q3:  Spin angle
        q4:  N[1] measure number of contact point position relative to origin
        q5:  N[2] measure number of contact point position relative to origin

    _params is the radius of the torus.

    Output is four 3-tuples in the following order:
          CO:  Postion of the center of the disc in inertial coordinates
        B[2]:  Vector normal to the plane of the disc (axis of rotation)
        C[1]:  1st body fixed unit vector in the plane of the disc
        C[3]:  3rd body fixed unit vector in the plane of the disc
"""
output_string += generate_function("anim", anim_eqs, q, params,\
        triples=True, docstring=ds)

ds = """\
Kinetic and Potential Energy of rolling torus.

    _x is an array/list in the following order:
        q2:  Lean angle
        u1:  B[1] measure number of angular velocity (roll rate)
        u2:  B[2] measure number of angular velocity (spin rate)
        u3:  B[3] measure number of angular velocity (yaw-like rate)

    _params is an array/list in the following order:
         m:  Disc mass
         g:  Gravitational constant
         r1: Major radius of torus
         r2: Minor radius of torus

    Returns a list/array of kinetic energy and potential energy, respectively.
"""
output_string += generate_function("energy", energy_eqs, q+u, params, docstring=ds)

file = open('rollingtorus_lib.py', 'w')
file.write(output_string)
file.close()
