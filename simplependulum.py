from math import sin, cos, pi

from numpy import array, arange
from sympy import symbols, Function, S, solve, simplify, \
        collect, Matrix, lambdify

from pydy import ReferenceFrame, cross, dot, dt, express, expression2vector, \
    coeff


# Declare parameters
m, g, l, t = symbols("m, g, l, t")
# Declare generalized coordinate speed as implicit functions of time
q1 = Function("q1")(t)
u1 = Function("u1")(t)
# Declare a symbols for:
# u1p:  Time derivative of generalized speed (u1p = d(u1)/dt)
# au1:  Auxilliary generalized speed (For determination of constraint force)
# cf1:  Constraint Force the time derivative of the generalized speed
u1p, au1, cf1 = symbols("u1p, au1, cf1")


# Declare an inertial ("Newtonian") reference frame
N = ReferenceFrame("N")
# N[3] is aligned with the local gravitational field, i.e., downward as drawn below
# N[2] is parallel to the axis of rotation of the hinge and comes out of the page as drawn below
# N[1] is to the right if N[2] is pointed out of the page as you would draw it:
#        NO
#    -----o--> N[1] -----
#         |\
#         | \
#         |  \
#         |   O    <-- PO (point mass)

P = N.rotate("P", 2, q1)
# Define the position from the hinge (taken to be the inertial origin) to the point
r_NO_PO = l * P[3]
print "R_NO_PO> = ", r_NO_PO

# Create a dictionary whose keys are symbols to be replaced with the values
q_u_dots = [q1.diff(t), au1, u1.diff(t)]
u_up = [u1, au1, u1p]
gen_speeds = dict(zip(q_u_dots, u_up))


P.set_omega(P.get_omega(N).subs(gen_speeds), N, force=True)
print "W_P_N> = ", P.get_omega(N)

# Velocity of PO relative to N, with the au1*P[3] term representing
# a component of velocity which it cannot actually posess.  This is for
# computation of the constraint force acting along the length of the rod.
V_PO_N = au1*P[3] + dt(r_NO_PO, N, t)
print "V_PO_N> = ", V_PO_N


#print r_NO_PO
#print v_PO_N
#print P.get_omega(N)

# Compute the partial velocities of the points NO and PO
VPO = coeff(V_PO_N, u_up[:-1])
print "Partial velocities:\n", VPO

# Now set the auxilliary speed to zero
V_PO_N = V_PO_N.subs({au1: 0})
print "V_PO_N> = ", V_PO_N

# Determine the acceleration of PO in N
A_PO_N = dt(V_PO_N, N, t)
A_PO_N = A_PO_N.subs(gen_speeds)
print "A_PO_N> = ", A_PO_N

# Forces acting at point PO
F_PO = cf1*P[3] + m*g*N[3]
print "F_PO> = ", F_PO

# Inertia Force acting at PO
RSTAR_PO = -m*A_PO_N
print "RSTAR_PO>", RSTAR_PO

FR = [dot(ViPO, F_PO) for ViPO in VPO]
FRSTAR = [dot(ViPO, RSTAR_PO) for ViPO in VPO]
print "FR = \n", FR
print "FRSTAR = \n", FRSTAR

# Kane's dynamic equations:
zero = Matrix(FR) + Matrix(FRSTAR)

# Solve the equations for the time derivative of 
# the independent generalized speed and the constraint force
eqns = [x for x in zero[:2]]
r = solve(eqns, [u1p, cf1])

print "_"*80
print "Equations of motion in first order form: "
print "q1p = ", gen_speeds[q1.diff(t)]
print "u1p = ", r[u1p]

print "_"*80
print "\nConstraint force acting in the P3> direction: "
print "cf1 = ", r[cf1]

# Numerical integration of the equations of motion


