#!/usr/bin/env python
from pydy import *

# Create a Newtonian reference frame
N = NewtonianReferenceFrame('N')

# Declare parameters, coordinates, speeds
params = N.declare_parameters('m g I11 I22 I33')
q, qd = N.declare_coords('q', 6)
u, ud = N.declare_speeds('u', 6)
# Unpack the lists
m, g, I11, I22, I33 = params
q1, q2, q3, q4, q5, q6 = q
q1d, q2d, q3d, q4d, q5d, q6d = qd
u1, u2, u3, u4, u5, u6 = u
u1d, u2d, u3d, u4d, u5d, u6d = ud

# Frame fixed to the rigid body
A = N.rotate("A", 'BODY312', (q1, q2, q3), I=(I11, I22, I33, 0, 0, 0))

# Locate the mass center
AO = N.O.locate('AO', q4*N[1] + q5*N[2] + q6*N[3], mass=m)

# Define the generalized speeds
u_rhs = [dot(A.ang_vel(), A[i]) for i in (1, 2, 3)] + \
        [dot(AO.vel(), N[i-3]) for i in (4, 5, 6)]

# Form the list of equations mapping qdots to generalized speeds
qd_to_u_eqs = [Eq(ui, ui_rhs) for ui, ui_rhs in zip(u, u_rhs)]
# Form the matrix that maps qdot's to u's
qd_to_u = coefficient_matrix(u_rhs, qd)

u_to_qd = qd_to_u.inv(method='ADJ').subs({cos(q3)**2: 1-sin(q3)**2}).expand()
print u_to_qd
qd_rhs = u_to_qd * Matrix(u)

# Create a list of Equations:
kd_eqs = []
for qdot, eqn in zip(qd, qd_rhs):
    kd_eqs.append(Eq(qdot, eqn.subs({sin(q2)/cos(q2):tan(q2)})))

# Set velocities and angular velocities using only generalized speeds
A.abs_ang_vel = Vector(u1*A[1] + u2*A[2] + u3*A[3])
AO.abs_vel = Vector(u4*N[1] + u5*N[2] + u6*N[3])

# Set accelerations and angular accelerations
A.abs_ang_acc = dt(A.abs_ang_vel, N)
AO.abs_acc = dt(AO.abs_vel, N)

# Apply gravity
N.gravity(g*N[3])

# Form Kane's equations and solve them for the udots
kanes_eqns = N.form_kanes_equations()

rhs = zeros((6,1))
for i in range(6):
    rhs[i] = kanes_eqns[i].rhs
dd = N.mass_matrix.inv() * rhs

dd_eqs = [Eq(lhs, rhs) for lhs, rhs in zip(ud, dd)]

# Energy calculations
ke = ((m*AO.abs_vel.mag_sqr + dot(A.abs_ang_vel, dot(A.inertia,
    A.abs_ang_vel))).expand())/2
pe = -m*g*AO.rel(N.O).dot(N[3])
energy_eqs= [Eq(Symbol('ke'), ke),
             Eq(Symbol('pe'), pe)]

# Animation equations
anim_eqs = animate(N, ('AO', AO.rel(N.O)), ('A1', A[1]), ('A3', A[3]))

output_string = "from __future__ import division\n"
output_string += "from math import sin, cos, tan\n\n"

ds = """\
Rigidy body equations of motion.

    _x is an array/list in the following order:
        q1:  Yaw   \
        q2:  Lean   |-(Euler 3-1-2 angles used to orient A
        q3:  Pitch /
        q4:  N[1] displacement of mass center.
        q5:  N[2] displacement of mass center.
        q6:  N[3] displacement of mass center.
        u1:  A[1] measure number of angular velocity
        u2:  A[2] measure number of angular velocity
        u3:  A[3] measure number of angular velocity
        u4:  N[1] velocity of mass center.
        u5:  N[2] velocity of mass center.
        u6:  N[3] velocity of mass center.

    _params is an array/list in the following order:
        m:  Mass of first pendulum point mass.
        g:  Gravitational constant.
      I11:  Principal moment of inertia about A[1]
      I22:  Principal moment of inertia about A[2]
      I33:  Principal moment of inertia about A[3]
"""

output_string += generate_function("eoms", kd_eqs+dd_eqs, q+u, params, docstring=ds,
        time=True)

ds = """\
Kinetic and Potential Energy of rigid body.

    _x is an array/list in the following order:
        q1:  Yaw   \
        q2:  Lean   |-(Euler 3-1-2 angles used to orient A
        q3:  Pitch /
        q4:  N[1] displacement of mass center.
        q5:  N[2] displacement of mass center.
        q6:  N[3] displacement of mass center.
        u1:  A[1] measure number of angular velocity
        u2:  A[2] measure number of angular velocity
        u3:  A[3] measure number of angular velocity
        u4:  N[1] velocity of mass center.
        u5:  N[2] velocity of mass center.
        u6:  N[3] velocity of mass center.

    _params is an array/list of:
        m:  Mass of first pendulum point mass.
        g:  Gravitational constant.
      I11:  Principal moment of inertia about A[1]
      I22:  Principal moment of inertia about A[2]
      I33:  Principal moment of inertia about A[3]

    Returns a list/array of kinetic energy and potential energy, respectively.
"""
output_string += generate_function("energy", energy_eqs, q+u, params, docstring=ds)

ds = """\
Calculate configuration of pendulum for purposes of animation.

    _x is an array/list of the configuration coordinates of the disc:

    _params is an array/list of:
        m:  Mass of first pendulum point mass.
        g:  Gravitational constant.
      I11:  Principal moment of inertia about A[1]
      I22:  Principal moment of inertia about A[2]
      I33:  Principal moment of inertia about A[3]

    Output is:
          AO:  Position of mass center.
        A[1]:  Body fixed A[1] unit vector
        A[3]:  Body fixed A[3] unit vector
"""
output_string += generate_function("anim", anim_eqs, q, params, docstring=ds,
        triples=True)

file = open('rigidbody_lib.py', 'w')
file.write(output_string)
file.close()
