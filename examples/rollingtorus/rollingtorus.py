#!/usr/bin/env python
from sympy import solve, simplify
from pydy import *

# Create a Newtonian reference frame
N = NewtonianReferenceFrame('N')

# Constants
m, g, r1, r2, I, J = N.declare_parameters("m g r1 r2 I J")

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

# Form the velocity of the rear wheel center in two different ways in order to
# bring the non holonomic constraints into evidence.
vcon1 = dt(CO.rel(N1), N)
vcon2 = cross(C.ang_vel(N), CO.rel(N.O))
eq1 = dot(vcon1 - vcon2, N[1])
eq2 = dot(vcon1 - vcon2, N[2])

T = coefficient_matrix([eq1, eq2], qd)
Td_inv, Ti = transform_matrix(T, qd, qd[3:])

contact_rates = Td_inv*Ti*Matrix(qd[:3])

# Append the contact point kinematic differential equations to the list,
# keeping them in implicit form.
kd.append(Eq(q4d, contact_rates[0]))
kd.append(Eq(q5d, contact_rates[1]))

# Generate a dictionary to be used for replacing qdots in accerlations with
# expressions involving the generalized speeds.
kd_dict = eqn_list_to_dict(kd)

# Set velocities and angular velocities using only generalized speeds
B.abs_ang_vel = Vector(u1*B[1] + u3*B[3])
C.abs_ang_vel = Vector(u1*B[1] + u2*B[2] + u3*B[3])
CO.abs_vel = C.abs_ang_vel.cross(CO.rel(N.O))

# Set accelerations and angular accelerations
C.abs_ang_acc = dt(C.abs_ang_vel, N)
CO.abs_acc = dt(CO.abs_vel, N).express(B)

# Apply gravity
N.gravity(g*A[3])

# Form Kane's equations and solve them for the udots
kanes_eqns = N.form_kanes_equations()
N.mass_matrix[0,0] = N.mass_matrix[0,0].expand().subs(N.csqrd_dict).expand()
print kanes_eqns[0].rhs
print kanes_eqns[1].rhs
print kanes_eqns[2].rhs

print N.mass_matrix
stop

dyndiffs = N.solve_kanes_equations()

print 'Dynamic differential equations'
for ud in udot_list[:3]:
    dyndiffs[ud] = dyndiffs[ud].expand()
    print ud, '=', dyndiffs[ud]
