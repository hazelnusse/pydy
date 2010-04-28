#!/usr/bin/env python
from pydy import *

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
M, Md = symbols('M Md')

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
#AO = ABO.locate('AO', P_ABO_AO.subs(l_dict), mass=ma)
#BO = ABO.locate('AO', P_ABO_BO.subs(l_dict), mass=mb)
# Overwrite previous definitions of AO and BO
AO = ABO.locate('AO', P_ABO_AO, mass=ma)
BO = ABO.locate('AO', P_ABO_BO, mass=mb)

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
#AO.abs_vel = ABO.abs_vel + cross(A.ang_vel(N), AO.rel(ABO)).subs(l_dict)
#BO.abs_vel = ABO.abs_vel + cross(A.ang_vel(N), BO.rel(ABO)).subs(l_dict)
AO.abs_vel = ABO.abs_vel + cross(A.ang_vel(N), AO.rel(ABO))#.subs(l_dict)
BO.abs_vel = ABO.abs_vel + cross(A.ang_vel(N), BO.rel(ABO))#.subs(l_dict)
# Set accelerations and angular accelerations
A.abs_ang_acc = dt(A.abs_ang_vel, N)
B.abs_ang_acc = dt(B.abs_ang_vel, N)
ABO.abs_acc = Vector(u5d*N[1] + u6d*N[2] + u7d*N[3])
AO.abs_acc = dt(AO.abs_vel, N)
BO.abs_acc = dt(BO.abs_vel, N)

# Apply gravity
N.gravity(g*N[3])
# Apply a torque between the two bodies
B.apply_torque(T*A[2], A)

# Form Kane's equations and solve them for the udots
kanes_eqns = N.form_kanes_equations()

print 'Dynamic differential equations'
for i in range(7):
    print kanes_eqns[i]

# Alternative formulation
AO.mass = S(0)
AO.force = Vector(0)
BO.mass = S(0)
BO.force = Vector(0)
ABO.mass = M
ABO.force = Vector(M*g*N[3])
IAB11, IAB22, IAB33, IAB13 = symbols('IAB11 IAB22 IAB33 IAB13')
A.inertia = Inertia(A, (IAB11, IAB22, IAB33, 0, 0, IAB13))
B.inertia = Inertia(A, (0, J, 0, 0, 0, 0))

# Form Kane's equations and solve them for the udots
kanes_eqns2 = N.form_kanes_equations()

print 'Dynamic differential equations combined'
for i in range(7):
    print kanes_eqns2[i]

I_AO_ABO = inertia_of_point_mass(ma, P_ABO_AO, A)
I_BO_ABO = inertia_of_point_mass(mb, P_ABO_BO, A)
I_A_AO = Dyad({A[1]*A[1]: I11, A[2]*A[2]: I22, A[3]*A[3]:I33, A[1]*A[3]:I13,
    A[3]*A[1]:I13})

# Planar components of B's inertia
I_B_BO_p = Dyad({A[1]*A[1]: I, A[3]*A[3]: I})
# Combined system inertia except for B's out of plane moment of inertia
I_SYS_ABO = I_A_AO + I_AO_ABO + I_B_BO_p + I_BO_ABO

print 'I_AB11', dot(A[1], dot(I_SYS_ABO, A[1]))
print 'I_AB22', dot(A[2], dot(I_SYS_ABO, A[2]))
print 'I_AB33', dot(A[3], dot(I_SYS_ABO, A[3]))
print 'I_AB13', dot(A[1], dot(I_SYS_ABO, A[3]))
