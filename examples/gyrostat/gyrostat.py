#!/usr/bin/env python
from pydy import *
from sympy import factor

# Create a Newtonian reference frame
N = NewtonianReferenceFrame('N')

# Declare parameters, coordinates, speeds
params = N.declare_parameters('l1 l2 l3 ma mb g I11 I22 I33 I12 I23 I13 I J K T')
q, qd = N.declare_coords('q', 7)
u, ud = N.declare_speeds('u', 7)
# Unpack the lists
l1, l2, l3, ma, mb, g, I11, I22, I33, I12, I23, I13, I, J, K, T = params
q1, q2, q3, q4, q5, q6, q7 = q
q1d, q2d, q3d, q4d, q5d, q6d, q7d = qd
u1, u2, u3, u4, u5, u6, u7 = u
u1d, u2d, u3d, u4d, u5d, u6d, u7d = ud

# Some extra symbols for convenience
l1a, l3a, l1b, l3b = symbols('l1a l3a l1b l3b')

# Frame fixed to the rigid body
A = N.rotate("A", 'BODY312', (q1, q2, q3), I=(I11, I22, I33, 0, 0, I13))
B = A.rotate("B", 2, q4, I=(I, J, I, 0, 0, 0), I_frame=A)

# Create the point AO
AO = Point('AO')
# Locate AO relative to BO
BO = AO.locate('BO', l1*A[1] + l3*A[3])
# Position from ABO to AO
P_ABO_AO = -mass_center(AO, [(AO, ma), (BO, mb)])
# Position from ABO to BO
P_ABO_BO = -mass_center(BO, [(AO, ma), (BO, mb)])

l_dict_r = {l1a: dot(P_ABO_AO, A[1]),
            l3a: dot(P_ABO_AO, A[3]),
            l1b: dot(P_ABO_BO, A[1]),
            l3b: dot(P_ABO_BO, A[3])}
l_dict = dict([(v, k) for k, v in l_dict_r.items()])

# Locate the mass center of the system
ABO = N.O.locate('ABO', q5*N[1] + q6*N[2] + q7*N[3])
# Overwrite previous definitions of AO and BO
AO = ABO.locate('AO', P_ABO_AO.subs(l_dict), mass=ma)
BO = ABO.locate('AO', P_ABO_BO.subs(l_dict), mass=mb)

# Define the generalized speeds
u_rhs = [dot(A.ang_vel(), A[i]) for i in (1, 2, 3)] + \
        [dot(B.ang_vel(), A[2])] + \
        [dot(ABO.vel(), N[i]) for i in (1, 2, 3)]

# Form the list of equations mapping qdots to generalized speeds
qd_to_u_eqs = [Eq(ui, ui_rhs) for ui, ui_rhs in zip(u, u_rhs)]
# Form the matrix that maps qdot's to u's
qd_to_u = coefficient_matrix(u_rhs, qd)
adj = qd_to_u.adjugate().expand().subs(N.csqrd_dict).expand()
det = qd_to_u.det(method="berkowitz").expand().subs(N.csqrd_dict).expand()
u_to_qd = (adj / det).expand()#.subs({sin(q2)**2:1-cos(q2)**2}).expand()

qd_rhs = u_to_qd * Matrix(u)
# Create a list of kinematic differential equations
u_to_qd_eqs = []
print 'Kinematic differential equations'
for qdot, eqn in zip(qd, qd_rhs):
    u_to_qd_eqs.append(Eq(qdot, eqn.subs({sin(q2)/cos(q2): tan(q2)})))
    print u_to_qd_eqs[-1]

# Set velocities and angular velocities using only generalized speeds
A.abs_ang_vel = Vector(u1*A[1] + u2*A[2] + u3*A[3])
B.abs_ang_vel = Vector(u1*A[1] + u4*A[2] + u3*A[3])
ABO.abs_vel = Vector(u5*N[1] + u6*N[2] + u7*N[3])
AO.abs_vel = ABO.abs_vel + cross(A.ang_vel(N), AO.rel(ABO)).subs(l_dict)
BO.abs_vel = ABO.abs_vel + cross(A.ang_vel(N), BO.rel(ABO)).subs(l_dict)
print AO.abs_vel
print BO.abs_vel
raw_input()
# Set accelerations and angular accelerations
A.abs_ang_acc = dt(A.abs_ang_vel, N)
B.abs_ang_acc = dt(B.abs_ang_vel, N)
AO.abs_acc = dt(AO.abs_vel, N)
BO.abs_acc = dt(BO.abs_vel, N)
print A.abs_ang_acc
raw_input()
print B.abs_ang_acc
raw_input()
print AO.abs_acc
raw_input()
print BO.abs_acc
raw_inpu()

stop
# Apply gravity
N.gravity(g*N[3])
# Apply a torque between the two bodies
B.apply_torque(T*A[2], A)

# Form Kane's equations and solve them for the udots
kanes_eqns = N.form_kanes_equations()

print 'Dynamic differential equations'
for i in range(7):
    print kanes_eqns[i]
#l1a, l2a, l3a, l1b, l2b, l3b = symbols('l1a l2a l3a l1b l2b l3b')
#l_sd = {l1a: dot(P_AO_ABO, A[1]),
#        l2a: dot(P_AO_ABO, A[2]),
#        l3a: dot(P_AO_ABO, A[3]),
#        l1b: dot(P_BO_ABO, A[1]),
#        l2b: dot(P_BO_ABO, A[2]),
#        l3b: dot(P_BO_ABO, A[3])}
#
#P_AO_ABO = Vector(l1a*A[1] + l3a*A[3])
#P_BO_ABO = Vector(l1b*A[1] + l3b*A[3])
#mab = Symbol('mab')
#m_subs_dict = {ma+mb:mab}

I_AO_ABO = inertia_of_point_mass(ma, -P_AO_ABO, A)
I_BO_ABO = inertia_of_point_mass(mb, -P_BO_ABO, A)

# Planar components of B's inertia
I_B_BO_p = Dyad({A[1]*A[1]: I, A[3]*A[3]: I})
# Combined system inertia
I_SYS_ABO = A.inertia + I_AO_ABO + I_B_BO_p + I_BO_ABO
#I_BO_ABO.dict.pop(A[2]*A[2])
#B.inertia.dict.pop(A[2]*A[2])
# Combined inertia of A and B about ABO, except out of plane components of B
I_SYS_ABO = A.inertia + I_AO_ABO + I_BO_ABO + B.inertia
