from sympy import collect
from pydy import *

# Declare a NewtonianReferenceFrame
N = NewtonianReferenceFrame('N')

# Declare parameters
rr, rf, lr, ls, lf, l1, l2, l3, l4, mcd, mef, IC22, ICD11, ID13, ICD33, ID22, IEF11, IE13, IEF33, IE22, IF22, g = N.declare_parameters('rr rf lr ls lf l1 l2 l3 l4 mcd mef IC22 ICD11 ID13 ICD33 ID22 IEF11 IE13 IEF33 IE22 IF22 g')
# Declare coordinates and their time derivatives
q, qd = N.declare_coords('q', 8)
# Declare speeds and their time derivatives
u, ud = N.declare_speeds('u', 16)
# Unpack the lists
q0, q1, q2, q3, q4, q5, q6, q7 = q
q0d, q1d, q2d, q3d, q4d, q5d, q6d, q7d = qd
u0, u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15 = u
u0d, u1d, u2d, u3d, u4d, u5d, u6d, u7d, u8d, u9d, u10d, u11d, u12d, u13d, u14d, u15d = u

# Reference Frames
# Yaw frame
A = N.rotate('A', 3, q0)
# Lean frame
B = A.rotate('B', 1, q1)
# Bicycle frame with rigidly attached rider
D = N.rotate('D', 'BODY312', (q0, q1, q3), I=(ICD11, ID22, ICD33, 0, 0, ID13))
# Rear wheel
C = N.rotate('C', 'BODY312', (q0, q1, q2), I=(0, IC22, 0, 0, 0, 0), I_frame=D)
# Fork / handle bar assembly
E = D.rotate('E', 3, q4, I=(IEF11, IE22, IEF33, 0, 0, IE13))
# Front wheel
F = E.rotate('F', 2, q5, I=(0, IF22, 0, 0, 0, 0), I_frame=E)

# Unit vector in the plane of the front wheel, pointed towards the ground
fo_fn_uv = Vector(N[3] - dot(E[2], N[3])*E[2]).normalized

# Some manual manipulations to express g in the E frame without expanding the
# coefficients
N3inE = N[3].express(E)
e1c_expr = rf*N3inE.dict[E[1]]*fo_fn_uv.dict[N[3]]
e3c_expr = rf*N3inE.dict[E[3]]*fo_fn_uv.dict[N[3]]
n1, d1 = e1c_expr.as_numer_denom()
n2, d2 = e3c_expr.as_numer_denom()
assert d1 == d2
# Compute time derivatives using quotient rule and prevent unneccessary expansion
e1c_expr_dt = (n1.diff(t)*d1 - n1*d1.diff(t))/d1**2
e3c_expr_dt = (n2.diff(t)*d2 - n2*d2.diff(t))/d2**2
# Create variables for to act as place holders in the vector for FO to FN
e1c = Function('e1c')(t)
e3c = Function('e3c')(t)
fo_fn_subs_dict = {e1c: e1c_expr, e3c: e3c_expr, e1c.diff(t): e1c_expr_dt,
        e3c.diff(t): e3c_expr_dt}
fo_fn = Vector({E[1]: e1c, E[3]: e3c})

# Locate rear wheel center relative to point fixed in N, coincident with rear
# wheel contact
CO = N.O.locate('CO', -rr*B[3], C)
# Locate mass center of ricycle with rigidly attached rider
CDO = CO.locate('CDO', l1*D[1] + l2*D[3], D, mass=mcd)
# Locate top of steer axis
DE = CO.locate('DE', lr*D[1], D)
# Locate front wheel center
FO = DE.locate('FO', lf*E[1] + ls*E[3], E)
# Locate mass center of fork/handlebar assembly
EFO = FO.locate('EFO', l3*E[1] + l4*E[3], E)
# Locate front wheel contact point (fixed in the front wheel)
FN = FO.locate('FN', fo_fn, F)
# Locate another point fixed in N
N1 = CO.locate('N1', rr*B[3] - q6*N[1] - q7*N[2])

# Form the kinematic constraints in terms of the qdots.
constraint_eqns = [dot(FN.vel(), A[1]).expand().subs(N.csqrd_dict),\
                   dot(FN.vel(), A[2]).expand().subs(N.csqrd_dict),\
                   dot(FN.vel(), A[3]).expand().subs(N.csqrd_dict),\
                   dot(N1.vel(), N[1]).expand().subs(N.csqrd_dict),\
                   dot(N1.vel(), N[2]).expand().subs(N.csqrd_dict)]

# Determine the constraint matrix for the qdots
con_matrix_qdots = coefficient_matrix(constraint_eqns, qd)

# Choose which qdots will be taken as independent.  For numerical integration,
# we often like to specify the initial condition in terms of the initial
# coordinate rates, rather than the initial speeds.
# Choose yaw rate, rear wheel rate, pitch rate, and translational rates as the
# independent qdots.  This implies the independent rates are the lean rate,
# steer rate and the front wheel rate.
qd_dep = [q0d, q2d, q3d, q6d, q7d]
adj, det, dep_ci, indep_ci = transform_matrix(con_matrix_qdots, qd, qd_dep)

stop
# Definitions of the generalized speeds in terms of time derivatives of
# coordinates
u_defs = N.define_speeds(
        [Eq(u0, dot(D.ang_vel(), D[1])),
         Eq(u1, dot(D.ang_vel(), D[2])),
         Eq(u2, dot(D.ang_vel(), D[3])),
         Eq(u5, dot(E.ang_vel(), E[3])),
         Eq(u6, dot(C.ang_vel(), D[2])),
         Eq(u7, dot(F.ang_vel(), E[2])),
         Eq(u14, dot(N1.vel(), N[1])),
         Eq(u15, dot(N1.vel(), N[2]))])

for u_def in u_defs:
    print u_def.lhs, ':=', u_def.rhs

# Solve u_defs for the qdots in terms of the u's
T, Tinv, kindiffs = N.form_kindiffs(u_defs, qdot_list)

stop

print 'Resulting kinematic differential equations'
for qd in qdot_list:
    print qd, '=', kindiffs[qd]
stop


# Rigid body angular velocities
# Angular velocity of rear frame with rigidly attached rider
D.abs_ang_vel = Vector(u0*D[1] + u1*D[2] + u2*D[3])
# Angular velocity of fork handlebar assembly
E.abs_ang_vel = Vector(u3*E[1] + u4*E[2] + u5*E[3])
# Angular velocity of rear wheel
C.abs_ang_vel = Vector(u0*D[1] + u6*D[2] + u2*D[3])
# Angular velocity of front wheel
F.abs_ang_vel = Vector(u3*E[1] + u7*E[2] + u5*E[3])

# Mass center velocities
CDO.abs_vel = Vector(u9*D[1] + u10*D[2] + u11*D[3])
EFO.abs_vel = Vector(u12*D[1] + u13*D[2] + u14*D[3])

# Motion constraints
# We have introduced 14 generalized speeds to describe the angular velocity of
# the rigid bodies and the velocity of the mass centers.  Along with these
# generalized speeds, there are 11 motion constraints which reduces the number
# of degrees of freedom of the system to 3, consistent with what we know about
# the bicycle.
# These constraints arise for two reasons:
# 1)  Nonholonomic rolling constraints
# 2)  Introduction of generalized speeds which make the dynamic equations
#     simpler, yet are dependent generalized speeds.

# Rear assembly mass center constraint (incorporates rolling constraints of
# rear wheel)
vcdon = cross(C.ang_vel(), CO.rel(N.O).express(D)) + cross(D.ang_vel(), CDO.rel(CO))
constraint_eqs = [Eq( u9 - dot(vcdon, D[1]), 0),\
                  Eq(u10 - dot(vcdon, D[2]), 0),\
                  Eq(u11 - dot(vcdon, D[3]), 0)]

# Front assembly mass center constraint (incorporates rolling constraints of
# front wheel)
vefon = cross(F.ang_vel(), FO.rel(FN)) + cross(E.ang_vel(), EFO.rel(FO))
constraint_eqs += [Eq(u12 - dot(vefon, E[1]), 0),\
                   Eq(u13 - dot(vefon, E[2]), 0),\
                   Eq(u14 - dot(vefon, E[3]), 0)]

# Motion constraints associated with point DE (top of steer axis)
vden1 = cross(C.ang_vel(), CO.rel(N.O).express(D)) + cross(D.ang_vel(),\
        DE.rel(CO))
vden2 = cross(F.ang_vel(), FO.rel(FN)) + cross(E.ang_vel(), DE.rel(FO))
constraint_eqs += [Eq(dot(vden1 - vden2, D[1]), 0),\
                   Eq(dot(vden1 - vden2, D[2]), 0),\
                   Eq(dot(vden1 - vden2, D[3]), 0)]

# Motion constraints associated with angular velocities of D and E
constraint_eqs += [Eq(u4 - dot(D.ang_vel(), E[1])),\
                   Eq(u5 - dot(D.ang_vel(), E[2]))]

# Form the constraint matrix for B*u = 0
B_con = N.form_constraint_matrix(constraint_eqs, u_list)

# Set the constraint matrix
N.set_constraint_matrix(B_con)

# Form speed transform matrix for ud = (adj/det)*ui = T*ui
adj, det = N.form_speed_transform_matrix(B_con, u_list[3:])

adj = adj.expand().subs(N.csqrd_dict).expand()
det = det.expand().subs(N.csqrd_dict).expand()

# Set the speed transform matrix and determine the dependent speeds
ud = N.set_speed_transform_matrix(adj, det)

# Set the kindiffs and dependent_speeds
N.setkindiffs(kindiffs, ud, acc=False)

# Setting the angular accelerations of the rigid bodies
C.abs_ang_acc = dt(C.ang_vel(), D) + cross(D.ang_vel(), C.ang_vel())
D.abs_ang_acc = dt(D.ang_vel(), D)
E.abs_ang_acc = dt(E.ang_vel(), E)
F.abs_ang_acc = dt(F.ang_vel(), E) + cross(E.ang_vel(), F.ang_vel())

# Setting the accelerations of the mass centers
CDO.abs_acc = dt(CDO.vel(), D) + cross(D.ang_vel(), CDO.vel())
EFO.abs_acc = dt(EFO.vel(), E) + cross(E.ang_vel(), EFO.vel())

# Apply gravity
N.gravity(g*N[3])

# Form Kane's equations and solve them for the udots
kanes_eqns = N.form_kanes_equations()
stop


dyndiffs = N.solve_kanes_equations()
N.setdyndiffs(dyndiffs)
stop
N.output_eoms('bicycle_eoms.py', CO, DO, DE, EO, FO, FN, C[1], C[2], C[3],
        D[1], D[2], D[3], E[1], E[2], E[3], F[1], F[2], F[3])
