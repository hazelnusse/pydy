from sympy import symbols, Function, S, solve, simplify, \
        collect, Matrix, lambdify, trigsimp, expand, Eq, pretty_print

from pydy import *

# Create a Newtonian reference frame
N = NewtonianReferenceFrame('N')

# Constants
m, g, r = N.declare_parameters("m g r")

I = m*r**2/4  # Central moment of inertia about any diameter
J = m*r**2/2  # Central moment of inertia about normal axis

# Declare generalized coordinates and generalized speeds
(q1, q2, q3, q4, q5), q_list, qdot_list = N.declare_coords('q', 5, list=True)
(u1, u2, u3), u_list, udot_list = N.declare_speeds('u', 3, list=True)

# Intermediate reference frames
A = N.rotate("A", 3, q1)
B = A.rotate("B", 1, q2)

# Frame fixed to the torus rigid body.
C = B.rotate("C", 2, q3, I=(I, J, I, 0, 0, 0), I_frame=B)

# Locate the mass center of torus
CO = N.O.locate('CO', -r*B[3], frame=C, mass=m)

# Fixed inertial reference point
N1 = CO.locate('N1', r*B[3] - q4*N[1] - q5*N[2])

# Define the generalized speeds to be the B frame measure numbers of the angular
# Must be of the form:  A*q' == u
u_defs = [Eq(dot(C.ang_vel(), B[i]), u_list[i-1]) for i in (1, 2, 3)]
# Must be of the form:  B*q' == 0
constrainteqs = [Eq(dot(N1.vel(), N[1]), 0), Eq(dot(N1.vel(), N[2]), 0)]
implicit_rates = solve(constrainteqs, q4.diff(t), q5.diff(t))

T, Tinv, kindiffs = N.form_transform_matrix(u_defs, qdot_list[:3], method='GE')
print 'Kinematic differential equations'
for qd in qdot_list[:3]:
    print qd, '=', kindiffs[qd]

print 'Implicitly defined rates:'
for qd in implicit_rates:
    print qd, '=', implicit_rates[qd]

A._wrel = express(A._wrel, B).subs(kindiffs)
B._wrel = express(B._wrel, B).subs(kindiffs)
C._wrel = express(C._wrel, B).subs(kindiffs)
#print C.ang_vel()
CO._vrel = cross(C.ang_vel(), CO.rel(N.O))
#print CO.vel()
#stop
# Substitute the kinematic differential equations into velocity expressions,
# form partial angular velocities and partial velocites, form angular
# accelerations and accelerations
N.setkindiffs(kindiffs, implicit_rates)
#N.setkindiffs({},{})

# Apply gravity
N.gravity(g*A[3])

# Form Kane's equations and solve them for the udots
kanes_eqns = N.form_kanes_equations()
dyndiffs = solve(kanes_eqns, udot_list)

print 'Dynamic differential equations'
for ud in udot_list:
    dyndiffs[ud] = dyndiffs[ud].expand()
    print ud, '=', dyndiffs[ud]
N.setdyndiffs(dyndiffs)
N.output_eoms('rollingdisc_eoms.py', (CO, N1), (B[2], q3), (C[1], 0), (C[3], 0))
