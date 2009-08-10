from sympy import (symbols, Function, S, solve, simplify,
        collect, Matrix, lambdify, trigsimp, expand, Eq, pretty_print,
        solve_linear_system, var)

from pydy import *
from scipy.integrate import odeint
from numpy import array, arange
import matplotlib.pyplot as plt

# Create a Newtonian reference frame
N = NewtonianReferenceFrame('N')

m, g, I11, I22, I33 = N.declare_parameters('m g I11 I22 I33')

# Declare generalized coordinates and generalized speeds
(q1, q2, q3, q4, q5, q6), q_list, qdot_list = N.declare_coords('q', 6)
(u1, u2, u3, u4, u5, u6), u_list, udot_list = N.declare_speeds('u', 6)

# Intermediate reference frames
A = N.rotate("A", 3, q1)
B = A.rotate("B", 1, q2)

# Frame fixed to the torus rigid body.
C = B.rotate("C", 2, q3, I=(I11, I22, I33, 0, 0, 0))

# Locate the mass center of torus
CO = N.O.locate('CO', q4*N[1] + q5*N[2] + q6*N[3], mass=m)

# Define the generalized speeds
u_rhs = [dot(C.ang_vel(), C[i]) for i in (1, 2, 3)]\
        + [dot(CO.vel(), N[i]) for i in (1, 2, 3)]

# Form transformation matrix of u = T * q'
T = N.form_transform_matrix(u_rhs, qdot_list)

# Need to implement a method to automate the solving for the qdots
system = T.col_insert(6,Matrix(u_list))
d = {}
system2 = zeros(system.shape)
for i in range(system.shape[0]):
    for j in range(system.shape[1]):
        if system[i, j] != 0:
            s = Symbol("a", dummy=True)
            d[s] = system[i, j]
            system2[i, j] = s

kindiffs = solve_linear_system(system2, *qdot_list)

for qd in qdot_list:
    kindiffs[qd] = trigsimp(kindiffs[qd].subs(d)).expand()

N.setkindiffs(kindiffs)

# Apply gravity
N.gravity(g*A[3])

# Form Kane's equations and solve them for the udots
kanes_eqns = N.form_kanes_equations()
dyndiffs = solve(kanes_eqns, udot_list)

N.setdyndiffs(dyndiffs)

for r, ud in enumerate(udot_list):
    dyndiffs[ud] = simplify(expand(dyndiffs[ud].subs(N.csqrd_dict)))

N.output_eoms('rigidbody_eoms.py')
