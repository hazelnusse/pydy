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
q0, q1, q2, q3, q4, q5 = q
q0d, q1d, q2d, q3d, q4d, q5d = qd
u0, u1, u2, u3, u4, u5 = u
u0d, u1d, u2d, u3d, u4d, u5d = ud

# Frame fixed to the rigid body
C = N.rotate("C", 'BODY312', (q0, q1, q2), I=(I11, I22, I33, 0, 0, 0))

# Locate the mass center
CO = N.O.locate('CO', q3*N[1] + q4*N[2] + q5*N[3], mass=m)

# Define the generalized speeds
u_rhs = [dot(C.ang_vel(), C[i]) for i in (1, 2, 3)] + \
        [dot(CO.vel(), C[i-3]) for i in (4, 5, 6)]

# Form the list of equations mapping qdots to generalized speeds
qd_to_u_eqs = [Eq(ui, ui_rhs) for ui, ui_rhs in zip(u, u_rhs)]
# Form the matrix that maps qdot's to u's
qd_to_u = coefficient_matrix(u_rhs, qd)

# Invert the matrix
u_to_qd_upper = qd_to_u[0:3,0:3].inverse_ADJ().subs(N.csqrd_dict).expand()
u_to_qd_lower_adj = qd_to_u[3:,3:].adjugate().expand().subs(N.csqrd_dict).expand()
u_to_qd_lower_det = qd_to_u[3:,3:].det(method="berkowitz").expand().subs(N.csqrd_dict).expand()
# Note: the above determinant should be 1
u_to_qd_lower = u_to_qd_lower_adj / u_to_qd_lower_det
u_to_qd = block_diag([u_to_qd_upper, u_to_qd_lower])

qd_rhs = u_to_qd * Matrix(u)

# Create a list of Equations:
u_to_qd_eqs = []
for qdot, eqn in zip(qd, qd_rhs):
    u_to_qd_eqs.append(Eq(qdot, eqn))

# Set velocities and angular velocities using only generalized speeds
C.abs_ang_vel = Vector(u0*C[1] + u1*C[2] + u2*C[3])
CO.abs_vel = Vector(u3*C[1] + u4*C[2] + u5*C[3])

# Set accelerations and angular accelerations
C.abs_ang_acc = dt(C.abs_ang_vel, C)
CO.abs_acc = dt(CO.abs_vel, C) + cross(C.abs_ang_vel, CO.abs_vel)

# Apply gravity
N.gravity(g*N[3])

# Form Kane's equations and solve them for the udots
kanes_eqns = N.form_kanes_equations()
dyndiffs = N.solve_kanes_equations()

print kanes_eqns
for i, d in enumerate(dyndiffs[3:]):
    rhs = d.rhs.expand()
    dyndiffs[i+3] = Eq(d.lhs, rhs)

print dyndiffs

