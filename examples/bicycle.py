from sympy import collect
from pydy import *

N = NewtonianReferenceFrame('N')

#rrt, rft, rr, rf, lr, ls, lf, l1, l2, l3, l4 = N.declare_parameters('rrt rft\
#        rr rf lr ls lf l1 l2 l3 l4')
rr, rf, lr, ls, lf, l1, l2, l3, l4 = N.declare_parameters('rr rf lr ls lf l1 l2 l3 l4')
(q1, q2, q3, q4, q5, q6), q_list, qdot_list = N.declare_coords('q', 6)
(u1, u2, u3, u4, u5, u6), u_list, udot_list = N.declare_speeds('u', 6)

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
C._wrel = Vector(q3.diff(t)*D[2])
F._wrel = Vector(q6.diff(t)*E[2])
#E._wrel_children[F] = Vector(-q6.diff(t)*E[2])

# Unit vector in the plane of the front wheel, pointed towards the ground
g = Vector(A[3] - (dot(E[2], A[3]))*E[2]).normalized

# Locate rear wheel center relative to point fixed in N, coincident with rear
# wheel contact
#CO = N.O.locate('CO', -rrt*A[3] - rr*B[3], C)
CO = N.O.locate('CO', - rr*B[3], C)
# Locate mass center of ricycle with rigidly attached rider
DO = CO.locate('DO', l1*D[1] + l2*D[3], D)
# Locate top of steer axis
DE = CO.locate('DE', lr*D[1], D)
# Locate mass center of fork/handlebar assembly
EO = DE.locate('EO', l3*E[1] + l4*E[3], E)
# Locate front wheel center
FO = DE.locate('FO', lf*E[1] + ls*E[3], E)
# Locate front wheel contact point (fixed in the front wheel)
#FN = FO.locate('FN', rf*g + rft*A[3], F)
FN = FO.locate('FN', rf*g, F)

u_defs = N.define_speeds(
        [Eq(u_list[i-1], dot(C.ang_vel(), B[i])) for i in (1, 2, 3)] + \
        [Eq(u_list[3], dot(D.ang_vel(B), B[2])),
         Eq(u_list[4], dot(E.ang_vel(D), D[3])),
         Eq(u_list[5], dot(F.ang_vel(E), E[2]))])

T, Tinv, kindiffs = N.form_kindiffs(u_defs, qdot_list)

print 'Kinematic differential equations'
for qd in qdot_list:
    print qd, '=', kindiffs[qd]

A._wrel = A._wrel.express(B).subs(kindiffs)
B._wrel = B._wrel.express(B).subs(kindiffs)
C._wrel = C._wrel.express(B).subs(kindiffs)
D._wrel = D._wrel.express(B).subs(kindiffs)
E._wrel = E._wrel.express(D).subs(kindiffs)
F._wrel = F._wrel.express(E).subs(kindiffs)

CO._vrel = cross(C.ang_vel(), CO.rel(N.O))
DO._vrel = cross(D.ang_vel(), DO.rel(CO))
DE._vrel = cross(D.ang_vel(), DE.rel(CO))
EO._vrel = cross(E.ang_vel(), EO.rel(DE))
FO._vrel = cross(E.ang_vel(), FO.rel(DE))
FN._vrel = cross(F.ang_vel(), FN.rel(FO))
"""
print A._wrel
print B._wrel
print C._wrel
print D._wrel
print E._wrel
print F._wrel

print C.ang_vel()
print D.ang_vel()
print E.ang_vel().express(E)
print F.ang_vel().express(E)
print CO._vrel
print DO._vrel
print DE._vrel
print EO._vrel
print FO._vrel
print FN._vrel
"""

constraint_eqs = [Eq(dot(FN.vel(), A[i]), 0) for i in (1,2,3)]

for eq in constraint_eqs:
    print eq
    for u in u_list:
        print collect(eq.lhs, u_list)
        print eq.lhs.coeff(u)
stop

B_constraints, T, dependent_speeds = N.impose_constraints(constraint_eqs,
        dependent=[u1,u3,u4])

print 'Dependent speeds:'
for u in u_list:
    if u in dependent_speeds:
        print u, '=', dependent_speeds[u]
