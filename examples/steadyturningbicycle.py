from sympy import solve
from pydy import *

N = NewtonianReferenceFrame('N')
rrt, rft, rr, rf, lr, ls, lf, l1, l2, l3, l4 = N.declare_parameters('rrt rft rr rf lr ls lf l1 l2 l3 l4')

# Declare generalized coordinates and generalized speeds
(q1, q2, q3, q4, q5, q6, q7, q8), q_list, qdot_list = N.declare_coords('q', 8, list=True)
(u1, u2, u3), u_list, udot_list = N.declare_speeds('u', 3, list=True)

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

# Locate rear wheel center relative to point fixed in N, coincident with rear
# wheel contact
CO = N.O.locate('CO', -rrt*A[3] - rr*B[3], C)

# Locate point fixed in N, to serve as inertial origin
N1 = CO.locate('NC', rrt*A[3] + rr*B[3] - q7*N[1] - q8*N[2])
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
hcs, dhcs = kinematic_chain(N.O, FN, vec_list=[A[3]])
stop
print hcs, dhcs
print hcs[0].diff(t)
stop
# Rear wheel nonholonomic constraint equations
nh_rear = [dot(N1.vel(), N[i]) for i in (1, 2)]
soln = solve(nh_rear, qdot_list_d[-2:])

# Front wheel nonholonomic constraint equations
nh_front = [Eq(0, dot(FN.vel(), A[i])) for i in (1, 2)]

u_rhs = [Eq(u1, dot(D.ang_vel(), B[1])),
         Eq(u2, dot(D.ang_vel(), B[2])),
         Eq(u3, dot(E.ang_vel(), D[3]))]
T, Tinv, kindiffs = N.form_transform_matrix(nh_front+u_rhs, qdot_list[:-2], expand=False)
print T
stop


N.set_nhc_eqns(nh_front)#, nh_rear)

# There are three generalized speeds but 6 qdots involved in nh_front and the
# differentiated kinematic chain

# Solve the matrix for the dependent rates
dependent_rates = N.solve_constraint_matrix(qdot_list_d[:-2])

print 'Dependent rates'
for dr in dependent_rates:
    print dr, ' = ', dependent_rates[dr]
print 'Implicitly defined dependent rates'
for dr in soln:
    print dr, ' = ', soln[dr]


u_rhs = [dot(D.ang_vel(), B[i]) for i in (1, 3)] + [dot(E.ang_vel(), D[3])]
print 'u_rhs = ', u_rhs
stop


# Create the equations that define the generalized speeds, then solve them for
# the time derivatives of the generalized coordinates
u_definitions = [u - u_r for u, u_r in zip(u_list, u_rhs)]
kindiffs = solve(u_definitions, qdot_list)
print 'Kinematic differential equations'

for qd in qdot_list[:3]:
    kindiffs[qd] = expand(kindiffs[qd])
    print qd, '=', kindiffs[qd]
    eoms.append(kindiffs[qd])

# Form the expressions for q1' and q2', taken to be dependent speeds
# Completely optional, these have no influence on the dynamics.
nh = [dot(N1.vel(), N[1]), dot(N1.vel(), N[2])]
dependent_rates = solve(nh, q4.diff(t), q5.diff(t))
print 'Dependent rates:'
for qd in dependent_rates:
    dependent_rates[qd] = expand(dependent_rates[qd].subs(kindiffs))
    print qd, '=', dependent_rates[qd]
    eoms.append(dependent_rates[qd])

# Substitute the kinematic differential equations into velocity expressions,
# form partial angular velocities and partial velocites, form angular
# accelerations and accelerations



N.setkindiffs(kindiffs)

N.gravity(g*A[3])

kanes_eqns = N.form_kanes_equations()

print kanes_eqns
