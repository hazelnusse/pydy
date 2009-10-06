
from sympy import solve
from pydy import *


N = NewtonianReferenceFrame('N')
rrt, rft, rr, rf, lr, ls, lf, l1, l2, l3, l4 = N.declare_parameters('rrt rft rr rf lr ls lf l1 l2 l3 l4')

# Declare generalized coordinates and generalized speeds
(q1, q2, q3, q4, q5, q6, q7, q8), q_list, qdot_list = N.declare_coords('q', 8)
(u1, u2, u3), u_list, udot_list = N.declare_speeds('u', 3, list=True)

q = [Symbol(str(qi.func)) for qi in q_list]
qp = [Symbol(str(qi.func)+'p') for qi in q_list]
q_dif_dict = dict(zip(q, qp))
q_sym_dict = dict(zip(q_list, q))


# Independent qdots
qdot_list_i = [q2.diff(t), q5.diff(t), q6.diff(t)]
# Take the dependent qdots to be yaw, rear wheel, pitch, and x, y rates
qdot_list_d = [q1.diff(t), q3.diff(t), q4.diff(t), q7.diff(t), q8.diff(t)]
# Dependent qdot's which are implicitly dependent upon the other three.
qdot_list_d_i = [q7.diff(t), q8.diff(t)]

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
g = Vector(A[3] - (dot(E[2], A[3]))*E[2]).normalized

v = Vector((rft-rrt)*A[3] - rr*B[3] + lr*D[1] + ls*D[3] + lf*E[1] + rf*g)
print v.dot(A[3])
stop

print 'v =',v
for uv, co in v.dict.items():
    print dot(A[3], uv), 'coeff: ', co
stop

print 'v[E[2]]*E[2] dot A[3]', dot(v.dict[E[2]]*E[2], A[3])
print 'v in A dot A[3]', express(v, A).dot(A[3])
print 'v dot A[3]', dot(v, A[3])
stop


# Locate rear wheel center relative to point fixed in N, coincident with rear
# wheel contact
CO = N.O.locate('CO', -rrt*A[3] - rr*B[3], C)

# Locate point fixed in N, to serve as inertial origin
N1 = CO.locate('NC', rrt*A[3] + rr*B[3] - q7*N[1] - q8*N[2])
# Locate mass center of bicycle with rigidly attached rider
DO = CO.locate('DO', l1*D[1] + l2*D[3], D)
# Locate top of steer axis
DE = CO.locate('DE', lr*D[1], D)
# Locate mass center of fork/handlebar assembly, fixed in E
EO = DE.locate('EO', l3*E[1] + l4*E[3], E)
# Locate front wheel center, fixed in E
FO = DE.locate('FO', lf*E[1] + ls*E[3], E)
# Locate front wheel contact point, fixed in F
FN = FO.locate('FN', rf*g + rft*A[3], F)

print dot(FN.rel(FO), A[3])
stop

hc = dot(FN.rel(N.O), A[3]).subs(q_sym_dict)
for a in hc.args:
    print a

stop

print hc.atoms()
print hc.atoms() & set(q)
def difft(expr):
    diff_set = hc.atoms() & set(q)
    d = S(0)
    for coord in diff_set:
        d += expr.diff(coord)*q_dif_dict[coord]
    return d

dhc = difft(hc)
print dhc
