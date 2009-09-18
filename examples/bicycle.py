from sympy import collect
from pydy import *

# Declare a NewtonianReferenceFrame
N = NewtonianReferenceFrame('N')

# Declare parameters
params = N.declare_parameters('rr rf lr ls lf l1 l2 l3 l4 mcd mef IC22 ICD11 ID13 ICD33 ID22 IEF11 IE13 IEF33 IE22 IF22 g')
# Declare coordinates and their time derivatives
q, qd = N.declare_coords('q', 8)
# Declare speeds and their time derivatives
u, ud = N.declare_speeds('u', 14)
# Unpack the lists
rr, rf, lr, ls, lf, l1, l2, l3, l4, mcd, mef, IC22, ICD11, ID13, ICD33, ID22, IEF11, IE13, IEF33, IE22, IF22, g = params
q0, q1, q2, q3, q4, q5, q6, q7 = q
q0d, q1d, q2d, q3d, q4d, q5d, q6d, q7d = qd
u0, u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13 = u
u0d, u1d, u2d, u3d, u4d, u5d, u6d, u7d, u8d, u9d, u10d, u11d, u12d, u13d = ud

# Create variables for to act as place holders in the vector from FO to FN
e1c = Function('e1c')(t)
e3c = Function('e3c')(t)
e1c_s, e3c_s = symbols('e1c e3c')

# Some simplifying symbols / trig expressions
s0, s1, s2, s3, s4, s5, c0, c1, c2, c3, c4, c5, e1c_s, e3c_s = symbols('s0 s1 \
        s2 s3 s4 s5 c0 c1 c2 c3 c4 c5 e1c_s e3c_s')

symbol_subs_dict = {sin(q0): s0,
                    cos(q0): c0,
                    sin(q1): s1,
                    cos(q1): c1,
                    cos(q2): c2,
                    sin(q2): s2,
                    cos(q3): c3,
                    sin(q3): s3,
                    cos(q4): c4,
                    sin(q4): s4,
                    cos(q5): c5,
                    sin(q5): s5,
                    e1c    : e1c_s,
                    e3c    : e3c_s}
trig_subs_dict = {c0**2: 1-s0**2,
                  c1**2: 1-s1**2,
                  c2**2: 1-s2**2,
                  c3**2: 1-s3**2,
                  c4**2: 1-s4**2,
                  c5**2: 1-s5**2}

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
kinematic_constraint_eqns = [dot(FN.vel(), A[1]).expand().subs(N.csqrd_dict).expand(),\
                   dot(FN.vel(), A[2]).expand().subs(N.csqrd_dict).expand(),\
                   dot(FN.vel(), A[3]).expand().subs(N.csqrd_dict).expand(),\
                   dot(N1.vel(), N[1]).expand().subs(N.csqrd_dict).expand(),\
                   dot(N1.vel(), N[2]).expand().subs(N.csqrd_dict).expand()]

# Determine the constraint matrix for the qdots
con_matrix_qdots = coefficient_matrix(kinematic_constraint_eqns, qd)

# Choose which qdots will be taken as independent.  For numerical integration,
# we often like to specify the initial condition in terms of the initial
# coordinate rates, rather than the initial generalized speeds.
# Choose yaw rate, rear wheel rate, pitch rate, and translational rates as the
# independent qdots.  This implies the independent rates are the lean rate,
# steer rate and the front wheel rate.

# Docstring for the generated function.
ds = """\
Linear mapping from lean rate, steer rate, front wheel rate to yaw rate, rear
    wheel rate, pitch rate, rear wheel contact point velocity in N[1] and N[2]
    directions.
"""

qd_dep = [q0d, q2d, q3d, q6d, q7d]
qd_indep = [q1d, q4d, q5d]
T_qd_adj, T_qd_det, dep_ci, indep_ci = transform_matrix(con_matrix_qdots, qd, qd_dep)

e1c_s, e3c_s = symbols('e1c e3c')
f_to_s = {e1c: e1c_s, e3c: e3c_s}
fo_fn_symbol_subs_dict = {e1c_s: e1c_expr, e3c_s: e3c_expr}
func_params = params[:5] + (q0, q1, q3, q4)
output_string = linear_transform(T_qd_adj.subs(f_to_s), func_params,\
        "dependent_qdot", T_qd_det.subs(f_to_s), x=qd_indep, y=qd_dep,\
        nested_terms=fo_fn_symbol_subs_dict, docstring=ds)

# Form mapping from u's to qdots
# Definitions of the generalized speeds in terms of time derivatives of
# coordinates
u_lhs = [u0, u1, u2, u5, u6, u7]
u_rhs = [dot(D.ang_vel(), D[1]),
         dot(D.ang_vel(), D[2]),
         dot(D.ang_vel(), D[3]),
         dot(E.ang_vel(), E[3]),
         dot(C.ang_vel(), C[2]),
         collect(dot(F.ang_vel(), E[2]).expand(), qd)]#,

qd_to_u = coefficient_matrix(u_rhs, qd[:-2])
qd_to_u_adj, qd_to_u_det = qd_to_u.adjugate(), qd_to_u.det()

qd_to_u_adj = qd_to_u_adj.expand().subs(N.csqrd_dict).expand()
qd_to_u_det = qd_to_u_det.expand().subs(N.csqrd_dict).expand()

u_to_qd = qd_to_u_adj / qd_to_u_det

# Form velocity of rear wheel center in two distinct but equivalent ways.  This
# allows for the rates of the rear wheel coordinates to be determined.
vco1 = dt(q6*N[1] + q7*N[2] - rr*B[3], N)
vco2 = C.ang_vel(N).cross(CO.rel(N.O))
eq1 = dot(vco1 - vco2, N[1]).expand().subs(N.csqrd_dict).expand()
eq2 = dot(vco1 - vco2, N[2]).expand().subs(N.csqrd_dict).expand()

from sympy import solve
xy_rates = solve([eq1, eq2], qd[6:])

# Create the kinematic differential equations as both a list of Sympy Eq
# objects, and as a dictionary whose keys are the qdots and whose values are
# the right hand sides of the kinematic differential equations.
kindiffs = []
kindiffs_dict = {}
kindiffs_rhs = u_to_qd*Matrix(u_lhs)
for i, rhs in enumerate(kindiffs_rhs):
    kindiffs.append(Eq(qd[i], rhs))
    kindiffs_dict[qd[i]] = rhs

qd67_eqs =[xy_rates[qd[6]].subs({qd[2]:kindiffs_dict[qd[2]]}).expand(),\
           xy_rates[qd[7]].subs({qd[2]:kindiffs_dict[qd[2]]}).expand()]

kindiffs_dict[qd[6]] = qd67_eqs[0]
kindiffs_dict[qd[7]] = qd67_eqs[1]

q6q7rows =  coefficient_matrix(qd67_eqs, u_lhs)
u_to_qd = u_to_qd.col_join(q6q7rows)

# Depends upon yaw, lean, pitch, steer
func_params = (q0, q1, q3, q4)
ds = """\
Linear mapping from generalized speeds to time derivatives of coordinates.
"""
output_string += linear_transform(u_to_qd, func_params,\
        "kindiffs", x=u_lhs, y=qd, docstring=ds)

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
CDO.abs_vel = Vector(u8*D[1] + u9*D[2] + u10*D[3])
EFO.abs_vel = Vector(u11*D[1] + u12*D[2] + u13*D[3])

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
speed_constraint_eqs = [ u8 - dot(vcdon, D[1]),\
                   u9 - dot(vcdon, D[2]),\
                  u10 - dot(vcdon, D[3])]

# Front assembly mass center constraint (incorporates rolling constraints of
# front wheel)
vefon = cross(F.ang_vel(), FO.rel(FN)) + cross(E.ang_vel(), EFO.rel(FO))
speed_constraint_eqs += [u11 - dot(vefon, E[1]),\
                   u12 - dot(vefon, E[2]),\
                   u13 - dot(vefon, E[3])]

# Motion constraints associated with point DE (top of steer axis)
vden1 = cross(C.ang_vel(), CO.rel(N.O).express(D)) + cross(D.ang_vel(),\
        DE.rel(CO))
vden2 = cross(F.ang_vel(), FO.rel(FN)) + cross(E.ang_vel(), DE.rel(FO))
speed_constraint_eqs += [dot(vden1 - vden2, D[1]),\
                   dot(vden1 - vden2, D[2]),\
                   dot(vden1 - vden2, D[3])]

# Motion constraints associated with angular velocities of D and E
speed_constraint_eqs += [u3 - dot(D.ang_vel(), E[1]),\
                   u4 - dot(D.ang_vel(), E[2])]

# Form the constraint matrix for B*u = 0
B_con = coefficient_matrix(speed_constraint_eqs, u)
B_con_s = B_con.subs(symbol_subs_dict)

# Use angular velocity measure numbers of front wheel as the independent speeds
u_indep = [u3, u5, u7]
u_dep = [u0, u1, u2, u4, u6, u8, u9, u10, u11, u12, u13]

# This part is extremely slow.  This is because the adjugate of 13x13 matrix is
# being formed and needs to be expanded, and each expression being expanded is
# very long and messy (before expansion).  Thing simplify greatly after
# expansion, so it is worth the wait.

UD = B_con_s[:,0:3].row_join(B_con_s[:,4]).row_join(B_con_s[:,6]).row_join(B_con_s[:,8:])
UI = B_con_s[:,3].row_join(B_con_s[:,5]).row_join(B_con_s[:,7])
UD_adj = UD.adjugate().expand().subs(trig_subs_dict).expand()
UD_det = UD.det().expand().subs(trig_subs_dict).expand()
stop


T_ud_adj, T_ud_det, dep_ci, indep_ci, T_ud_dict = transform_matrix(B_con_s, u, u_dep,
        subs_dict=True)
stop

#T_ud_adj_sym = T_ud_adj.subs(symbol_subs_dict).expand().subs(trig_subs_dict).expand()
#T_ud_det_sym = factor(T_ud_det.subs(symbol_subs_dict).expand().subs(trig_subs_dict).expand())
#m, n = T_ud_adj.shape
#for i in range(m):
#    for j in range(n):
#        T_ud_adj_sym[i,j] = factor(T_ud_adj_sym[i,j])
stop

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
