from math import sin, cos, pi

from numpy import array, arange
from sympy import symbols, Function, S, solve, simplify, \
        collect, Matrix, lambdify

from pydy import * 

# Declare parameters
m, g, l, t = symbols("m, g, l, t")
# Declare generalized coordinate speed as implicit functions of time
q1 = Function("q1")(t)
u1 = Function("u1")(t)
# Declare a symbols for:
# u1p:  Time derivative of generalized speed (u1p = d(u1)/dt)
# au1:  Auxilliary generalized speed (For determination of constraint force)
# cf1:  Constraint Force the time derivative of the generalized speed
au1, cf1 = symbols("au1, cf1")


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

A = N.rotate("A", 2, q1)
# Define the position from the hinge (taken to be the inertial origin) to the point
P = N.O.locate('P', l*A[3])

# Create a dictionary whose keys are symbols to be replaced with the values

kindiffs = {q1.diff(t): u1}

P.vel = P.vel.subs(kindiffs) + au1*A[3]
A.W[N] = A.W[N].subs(kindiffs)

print "W_A_N> = ", A.W[N]
print "V_P_N> = ", P.vel

partial_v = P.vel.partials([u1, au1])
print partial_v

# Velocity of PO relative to N, with the au1*P[3] term representing
# a component of velocity which it cannot actually posess.  This is for
# computation of the constraint force acting along the length of the rod.
# Now set the auxilliary speed to zero
P.vel = P.vel.subs({au1: 0})

# Determine the acceleration of PO in N
P.acc = dt(P.vel, N).subs(kindiffs)
A.alf = dt(A.W[N], N)
print "A_P_N> = ", P.acc
print "W_A_N> = ", A.alf

# Forces acting at point PO
P.force = Vector(cf1*A[3] + m*g*N[3])
print "F_PO> = ", P.force

# Inertia Force acting at PO
R_star = InertiaForce(m, P.acc)
print "RSTAR_PO>", R_star

GAF = [dot(vp, P.force) for vp in partial_v]
GIF = [dot(vp, R_star) for vp in partial_v]
print "Generalized Active Forces:\n", GAF
print "Generalized Inertia Forces:\n", GIF

EOMS = [AF + IF for AF, IF in zip(GAF, GIF)]
dyn_eqns = solve(EOMS, u1.diff(t), cf1)
print 'Kinematic equations:\n', kindiffs
print 'Dynamic equations:\n', dyn_eqns


def f(x, t0, args):
    sp_dict = {q1: x[0], u1: x[1], m: args[0], g: args[1], l: args[2]}
    return [kindiffs[q1.diff(t)].subs(sp_dict),
            dyn_eqns[u1.diff(t)].subs(sp_dict)]

print f([1,0], 0, (1, 1, 1))
stop
from scipy.integrate import odeint
from numpy import arange
t = arange(0,1.1,0.1)
x0 = [0.1, 0.0]
x = odeint(f, x0, t)
print x

