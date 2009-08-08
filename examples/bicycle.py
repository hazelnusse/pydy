from sympy import *
from pydy import *

rrt, rft, rr, rf, lr, ls, lf = symbols('rrt rft rr rf lr ls lf')
l1, l2, l3, l4 = symbols('l1 l2 l3 l4')
(q1, q2, q3, q4, q5, q6), q_list, qdot_list = gcs('q', 6, list=True)
(u1, u2, u3), u_list, udot_list = gcs('u', 3, list=True)

N = NewtonianReferenceFrame('N')
N.setcoords(q_list, qdot_list, u_list, udot_list)

# Reference Frames
# Yaw frame
A = N.rotate('A', 3, q1)
# Lean frame
B = A.rotate('B', 1, q2)
# Rear wheel spin frame
C = B.rotate('C', 2, q3)
# Bicycle frame pitch frame
D = B.rotate('D', 2, q4)
# Steer frame
E = D.rotate('E', 3, q5)
# Front wheel spin frame
F = E.rotate('F', 2, q6)

# In last reference frame, use E[2] instead of F[2] for the angular velocity,
# this prevents the ignorable coordinate q8 from appearing in the nonholonomic
# constraint equations.
F._wrel = Vector(q6.diff(t)*E[2])
E._wrel_children[F] = Vector(-q6.diff(t)*E[2])

# Unit vector in the plane of the front wheel, pointed towards the ground
g = Vector(A[3] - (dot(E[2], A[3]))*E[2]).normalized()

# Locate rear wheel center relative to point fixed in N, coincident with rear
# wheel contact
CO = N.O.locate('CO', -rrt*A[3] - rr*B[3], C)

# Second point fixed in N, serves as an inertial reference point
#N1 = CO.locate('N1', rr*B[3] + rrt*A[3] - q7*N[1] - q8*N[2])

# Locate mass center of ricycle with rigidly attached rider
DO = CO.locate('DO', l1*D[1] + l2*D[3], D)
# Locate top of steer axis
DE = CO.locate('DE', lr*D[1], D)
# Locate mass center of fork/handlebar assembly
EO = DE.locate('EO', l3*E[1] + l4*E[3], E)
# Locate front wheel center
FO = DE.locate('FO', lf*E[1] + ls*E[3], E)
# Locate front wheel contact point (fixed in the front wheel)
FN = FO.locate('FN', rf*g + rft*A[3], F)

# Form the holonomic constraint and its time derivative
kinematic_chain(N.O, FN, vec_list=[A[3]])

# Front wheel nonholonomic constraint equations
nh_front = [dot(FN.vel(), A[i]) for i in (1, 2)]

N.set_nhc_eqns(nh_front)

# Take the dependent qdots to be yaw rate, rear wheel rate, pitch rate, and
# translation rates in the ground plane
dependent_rates = N.solve_constraint_matrix([q1.diff(t), q3.diff(t),
    q4.diff(t)])

print dependent_rates
