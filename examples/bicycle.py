from sympy import *
from pydy import *

rrt, rft, rr, rf, lr, ls, lf = symbols('rrt rft rr rf lr ls lf')
l1, l2, l3, l4 = symbols('l1 l2 l3 l4')
(q1, q2, q3, q4, q5, q6, q7, q8), q_list, qdot_list = gcs('q', 8, list=True)

N = NewtonianReferenceFrame('N')
N.q_list = q_list
N.qdot_list = qdot_list

# Reference Frames
A = N.rotate('A', 3, q3)
B = A.rotate('B', 1, q4)
C = B.rotate('C', 2, q5)
D = B.rotate('D', 2, q6)
E = D.rotate('E', 3, q7)
F = E.rotate('F', 2, q8)

# Unit vector in the plane of the front wheel, pointed towards the ground
g = Vector(A[3] - (dot(E[2], A[3]))*E[2]).normalized()

# Locate rear wheel center relative to point fixed in N, coincident with rear
# wheel contact
CO = N.O.locate('CO', -rrt*A[3] - rr*B[3], C)

# Locate point moving in N, coincident with rear wheel contact
NC = CO.locate('NC', rrt*A[3] + rr*B[3])
# Second point fixed in N, serves as an inertial reference point
N2 = NC.locate('N2', -q1*N[1] - q2*N[2])

# Rear wheel nonholonomic constraint equations
nh_rear = [dot(N.O.vel(N2, N), A[1]), dot(N.O.vel(N2, N), A[2])]
q1q2_kindiffs = solve(nh_rear, q1.diff(t), q2.diff(t))
print q1q2_kindiffs
stop

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
stop


loop = Vector(-rrt*A[3] - rr*B[3] + lr*D[1] + ls*D[3] + lf*E[1] + rf*g +
        rft*A[3])
print 'loop =', loop
print loop.dict
