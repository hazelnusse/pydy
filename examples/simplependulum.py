from math import sin, cos, pi

from numpy import array, arange
from sympy import symbols, Function, S, solve, simplify, \
        collect, Matrix, lambdify, ccode, pretty

from pydy import *

# Declare parameters
# au1:  Auxilliary generalized speed (For determination of constraint force)
# cf1:  Constraint Force measure number in A[3] direction acting on particle P
m, g, l, t = symbols('m, g, l, t')
au1, cf1 = symbols('au1, cf1')

# Declare generalized coordinate and speed as implicit functions of time
q1 = GeneralizedCoordinate('q1')
u1 = GeneralizedCoordinate('u1')

u_list = [u1, au1]

# Declare an inertial ("Newtonian") reference frame
N = NewtonianReferenceFrame("N")

# N[3] is aligned with the local gravitational field, i.e., downward as drawn below
# N[2] is parallel to the axis of rotation of the hinge and comes out of the page as drawn below
# N[1] is to the right if N[2] is pointed out of the page as you would draw it:
#        NO
#    -----o--> N[1] -----
#         |\
#         | \
#         |  \
#         |   O    <-- P (point mass)

A = N.rotate('A', 2, q1)
# Define the position from the hinge (taken to be the inertial origin) to the point
P = N.O.locate('P', l*A[3], mass=m, force=m*g*N[3])

# Create a dictionary whose keys are symbols to be replaced with the values

kode = {q1.diff(t): u1}

N.subkindiffs(kode)
N.setgenspeeds(u_list)

print 'W_A_N> =', A.ang_vel(N)
print 'V_P_N> =', P.vel()
stop

partial_v = P.vel[N].partials([u1, au1])
print partial_v

# Velocity of PO relative to N, with the au1*P[3] term representing
# a component of velocity which it cannot actually posess.  This is for
# computation of the constraint force acting along the length of the rod.
# Now set the auxilliary speed to zero
P.vel[N] = P.vel[N].subs({au1: 0})

# Determine the acceleration of PO in N
P.acc = dt(P.vel[N], N).subs(kindiffs)
A.alf = dt(A.W[N], N)
print 'A_P_N> =', P.acc
print 'W_A_N> =', A.alf

# Forces acting at point PO
P.force = Vector(cf1*A[3] + m*g*N[3])
print 'F_PO> =', P.force

# Inertia Force acting at PO
R_star = InertiaForce(m, P.acc)
print 'RSTAR_PO> =', R_star

GAF = [dot(vp, P.force) for vp in partial_v]
GIF = [dot(vp, R_star) for vp in partial_v]
print 'Generalized Active Forces:\n', GAF[0], ' ', GAF[1]
print 'Generalized Inertia Forces:\n', GIF[0], ' ',  GAF[1]

EOMS = [AF + IF for AF, IF in zip(GAF, GIF)]
dyn_eqns = solve(EOMS, u1.diff(t), cf1)
print 'Kinematic equations:\n', kindiffs
print 'Dynamic equations:\n', dyn_eqns
stop

def f(x, t0, args):
    print args
    sp_dict = {q1: x[0], u1: x[1], m: args[0], g: args[1], l: args[2]}
    return [kindiffs[q1.diff(t)].subs(sp_dict).n(),
            dyn_eqns[u1.diff(t)].subs(sp_dict).n()]

print f([1,0], 0, (1, 1, 1))
#stop

from scipy.integrate import odeint
from numpy import arange
t = arange(0,1.1,0.1)
x0 = [0.1, 0.0]
x = odeint(f, x0, t, (1, 1, 1))
print x

