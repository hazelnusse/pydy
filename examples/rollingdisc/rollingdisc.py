#!/usr/bin/env python
from sympy import solve, simplify
from pydy import *

# Create a Newtonian reference frame
N = NewtonianReferenceFrame('N')

# Constants
params = N.declare_parameters("m g r")

# Declare generalized coordinates and generalized speeds
q, qd = N.declare_coords('q', 5)
u, ud = N.declare_speeds('u', 3)

# Unpack the lists
q1, q2, q3, q4, q5 = q
q1d, q2d, q3d, q4d, q5d = qd
u1, u2, u3 = u
u1d, u2d, u3d = ud
m, g, r = params

# Dictionary used for substituting sin(q1)/cos(q1) with tan(q1)
tan_lean = {sin(q2)/cos(q2): tan(q2)}

# Define the inertia of the rolling disc
I = m*r**2/4  # Central moment of inertia about any diameter
J = m*r**2/2  # Central moment of inertia about normal axis

# Intermediate reference frames
A = N.rotate("A", 3, q1)
B = A.rotate("B", 1, q2)

# Frame fixed to the disc.
C = B.rotate("C", 2, q3, I=(I, J, I, 0, 0, 0), I_frame=B)

# Locate the mass center of disc
CO = N.O.locate('CO', -r*B[3], frame=C, mass=m)

# Fixed inertial reference point
N1 = CO.locate('N1', r*B[3] - q4*N[1] - q5*N[2])

# Define the generalized speeds to be the B frame measure numbers of the
# angular velocity of C relative to the inertial frame.
u_rhs = [dot(C.ang_vel(N), B[1]),
         dot(C.ang_vel(N), B[2]),
         dot(C.ang_vel(N), B[3])]

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

# Form the velocity of the rear wheel center in two different ways in order to
# bring the non holonomic constraints into evidence.
vcon1 = dt(CO.rel(N1), N)
vcon2 = cross(C.ang_vel(N), CO.rel(N.O))
eq1 = dot(vcon1 - vcon2, N[1])
eq2 = dot(vcon1 - vcon2, N[2])

contact_rates = solve([eq1, eq2], [q4d, q5d])
# Append the contact point kinematic differential equations to the list,
# keeping them in implicit form.
kd.append(Eq(q4d, contact_rates[q4d]))
kd.append(Eq(q5d, contact_rates[q5d]))

# Generate a dictionary to be used for replacing qdots in accerlations with
# expressions involving the generalized speeds.
kd_dict = eqn_list_to_dict(kd)

# Set velocities and angular velocities using only generalized speeds
C.abs_ang_vel = Vector(u1*B[1] + u2*B[2] + u3*B[3])
CO.abs_vel = C.abs_ang_vel.cross(CO.rel(N.O))

# Set accelerations and angular accelerations
C.abs_ang_acc = dt(C.abs_ang_vel, N).express(B).subs(kd_dict).subs(tan_lean)
CO.abs_acc = dt(CO.abs_vel, N).express(B).subs(kd_dict).subs(tan_lean)

# Apply gravity
N.gravity(g*N[3])

# Form Kane's equations and solve them for the udots
kanes_eqns = N.form_kanes_equations()
dyndiffs = N.solve_kanes_equations()
# End derivation of equations of motion

# Upright stability analysis -- spin rate as a parameter
##############################################################################
min_eqs = [kd[1]] + dyndiffs
q2s = N.q_list_s[1]
u1s, u2s, u3s = N.u_list_s
# Convert the right hand side only use Sympy Symbols instead of Functions
for i, eq in enumerate(min_eqs):
    lhs = min_eqs[i].lhs.subs(N.symbol_dict)
    rhs = min_eqs[i].rhs.subs(N.symbol_dict)
    min_eqs[i] = Eq(lhs, rhs)

x = [N.q_list_s[1]] + N.u_list_s
# Point of linearization:
lin_diff_eqs = []
J = zeros((4,4))
for i in range(4):
    for j in range(4):
        J[i, j] = min_eqs[i].rhs.diff(x[j])
# Steady upright conditions
J_s = J.subs({q2s: 0, u1s: 0, u3s: 0})
evals = J_s.eigenvals()
eval_eqs = []
for i, ev in enumerate(evals.keys()):
    eval_eqs.append(Eq(Symbol("ev%d"%(i+1)), ev))
critical_speed = [Eq(Symbol('cs'), solve(eval_eqs[0].rhs**2, [u2s])[0])]

##############################################################################


# Energy calculations
ke = ((m*CO.abs_vel.mag_sqr + dot(C.abs_ang_vel, dot(C.inertia, C.abs_ang_vel)))).expand()/2.
pe = -m*g*CO.rel(N.O).dot(N[3]) - m*g*r  # Take 0 to be when the disc is upright
energy_eqs= [Eq(Symbol('ke'), ke),
             Eq(Symbol('pe'), pe)]

# Begin code for generation of output functions
output_string = "from __future__ import division\n"
output_string += "from numpy import sin, cos, tan\n\n"
ds = """\
Equations of motion for rolling disc.

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
         r:  Radius of disc
"""

output_string += generate_function("eoms", kd+dyndiffs, q+u,\
        params, docstring=ds, time=True)


# Animation equations
anim_eqs = animate(N, ('CO', CO.rel(N1)), ('B2', B[2]), ('C1', C[1]),\
        ('C3', C[3]))
ds = """\
Calculate position and orientation of disc for purposes of animation.

    _x is an array/list of the configuration coordinates of the disc:
        q1:  Yaw angle (0 is aligned with N[1])
        q2:  Lean angle (0 is upright)
        q3:  Spin angle
        q4:  N[1] measure number of contact point position relative to origin
        q5:  N[2] measure number of contact point position relative to origin

    _params is the radius of the disc.

    Output is four 3-tuples in the following order:
          CO:  Postion of the center of the disc in inertial coordinates
        B[2]:  Vector normal to the plane of the disc (axis of rotation)
        C[1]:  1st body fixed unit vector in the plane of the disc
        C[3]:  3rd body fixed unit vector in the plane of the disc
"""
output_string += generate_function("anim", anim_eqs, q, [r],\
        triples=True, docstring=ds)
# Mapping between qdot and u, in case it is desired to specify initial
# coordinate rates rather than initial body fixed angular velocity measure
# numbers
ds = """\
Mapping from time derivates of coordinates to generalized speeds.

    _x is an array/list in the following order:
         q2:  Lean angle
        q1d:  Yaw rate
        q2d:  Lean rate
        q3d:  Spin rate
"""
output_string += generate_function("qdot_to_u", qd_to_u_eqs, [q2, q1d, q2d,
    q3d], docstring=ds)

# Output the eigenvalue equations, taken to be a function of spin rate and
# parameters
ds = """\
Eigenvalues of the linearized equations of motion about the upright
configuration with spin rate treated as a parameter.
"""
output_string += generate_function("evals", eval_eqs, [u2], [g, r],
        docstring=ds)

ds = """\
Critical speed of rolling disc.

    _x is an array/list in the following order:
        g:  Gravitational constant
        r:  Radius of disc
"""
output_string += generate_function("critical_speed", critical_speed, [g,r], docstring=ds)

ds = """\
Kinetic and Potential Energy of rolling disc.

    _x is an array/list in the following order:
        q2:  Lean angle
        u1:  B[1] measure number of angular velocity (roll rate)
        u2:  B[2] measure number of angular velocity (spin rate)
        u3:  B[3] measure number of angular velocity (yaw-like rate)

    _params is an array/list in the following order:
         m:  Disc mass
         g:  Gravitational constant
         r:  Radius of disc

    Returns a list/array of kinetic energy and potential energy, respectively.
"""
output_string += generate_function("energy", energy_eqs, q+u, params, docstring=ds)

file = open('rollingdisc_lib.py', 'w')
file.write(output_string)
file.close()
