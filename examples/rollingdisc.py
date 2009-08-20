from sympy import solve, simplify
from pydy import *

# Create a Newtonian reference frame
N = NewtonianReferenceFrame('N')

# Constants
m, g, r = N.declare_parameters("m g r")

I = m*r**2/4  # Central moment of inertia about any diameter
J = m*r**2/2  # Central moment of inertia about normal axis

# Declare generalized coordinates and generalized speeds
(q1, q2, q3, q4, q5), q_list, qdot_list = N.declare_coords('q', 5, list=True)
(u1, u2, u3, u4, u5), u_list, udot_list = N.declare_speeds('u', 5, list=True)

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
u_defs = N.define_speeds(
        [Eq(u_list[i-1], dot(C.ang_vel(), B[i])) for i in (1, 2, 3)] + \
        [Eq(u4, q4.diff(t)), Eq(u5, q5.diff(t))])

T, Tinv, kindiffs = N.form_kindiffs(u_defs, qdot_list, method='ADJ')

# Set angular velocity and velocity expressions to only involve generalized
# speeds
A._wrel = A._wrel.express(B).subs(kindiffs)
B._wrel = B._wrel.express(B).subs(kindiffs)
C._wrel = C._wrel.express(B).subs(kindiffs)
CO._vrel = cross(C.ang_vel(), CO.rel(N.O))
# A little tricky, but basically, the velocity of N1 relative to the disc
# center CO, as viewed by and observer fixed in N
N1._vrel = Vector(-u4*N[1] - u5*N[2]) + dt(N.O.rel(CO), N)

# Must be of the form:  B*u == 0
constrainteqs = [Eq(dot(N1.vel(), N[1]), 0), Eq(dot(N1.vel(), N[2]), 0)]
B, T, dependent_speeds = N.impose_constraints(constrainteqs, dependent=[u4,u5])

print 'Kinematic differential equations'
for qd in qdot_list:
    print qd, '=', kindiffs[qd]

print 'Dependent speeds:'
for u in u_list:
    if u in dependent_speeds:
        print u, '=', dependent_speeds[u]

# Set the kindiffs and dependent_speeds
N.setkindiffs(kindiffs, dependent_speeds)

# Apply gravity
N.gravity(g*A[3])

# Form Kane's equations and solve them for the udots
kanes_eqns = N.form_kanes_equations()
dyndiffs = solve(kanes_eqns, udot_list)

print 'Dynamic differential equations'
for ud in udot_list[:3]:
    dyndiffs[ud] = dyndiffs[ud].expand()
    print ud, '=', dyndiffs[ud]

N.setdyndiffs(dyndiffs)
N.output_eoms('rollingdisc_eoms.py', (CO, N1), (B[2], q3), (C[1], 0), (C[3], 0))
