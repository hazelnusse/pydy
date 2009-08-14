from sympy import (symbols, Function, S, solve, simplify,
        collect, Matrix, lambdify, trigsimp, expand, Eq, pretty_print,
        solve_linear_system, var)

from pydy import *
from scipy.integrate import odeint
from numpy import array, arange

# Create a Newtonian reference frame
N = NewtonianReferenceFrame('N')

m, g, I11, I22, I33 = N.declare_parameters('m g I11 I22 I33')

# Declare generalized coordinates and generalized speeds
(q1, q2, q3, q4, q5, q6), q_list, qdot_list = N.declare_coords('q', 6)
(u1, u2, u3, u4, u5, u6), u_list, udot_list = N.declare_speeds('u', 6)

# Intermediate reference frames
A = N.rotate("A", 3, q1)
B = A.rotate("B", 1, q2)

# Frame fixed to the torus rigid body
C = B.rotate("C", 2, q3, I=(I11, I22, I33, 0, 0, 0))

# Locate the mass center
CO = N.O.locate('CO', q4*N[1] + q5*N[2] + q6*N[3], mass=m)

# Define the generalized speeds
u_defs = [Eq(dot(C.ang_vel(), C[i]), u_list[i-1]) for i in (1, 2, 3)]\
        + [Eq(dot(CO.vel(), N[i-3]), u_list[i-1]) for i in (4, 5, 6)]

# Form transformation matrix of u := T * q'
T, Tinv, kindiffs = N.form_transform_matrix(u_defs, qdot_list, method='ADJ')

print 'Kinematic differential equations'
for r, qd in enumerate(qdot_list):
    kindiffs[qd] = trigsimp(expand(kindiffs[qd]))
    print qd, ' = ', kindiffs[qd]

N.setkindiffs(kindiffs)

# Apply gravity
N.gravity(g*A[3])

# Form Kane's equations and solve them for the udots
kanes_eqns = N.form_kanes_equations()
dyndiffs = solve(kanes_eqns, udot_list)
print 'Dynamic differential equations'
for r, ud in enumerate(udot_list):
    dyndiffs[ud] = simplify(expand(dyndiffs[ud].subs(N.csqrd_dict)))
    print ud, ' = ', dyndiffs[ud]

N.setdyndiffs(dyndiffs)

N.output_eoms('rigidbody_eoms.py', (CO, N.O), (C[2], q3))
