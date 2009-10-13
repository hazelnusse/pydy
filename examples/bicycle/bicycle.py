from sympy import collect, Function, solve, asin
from pydy import *

# Declare a NewtonianReferenceFrame
N = NewtonianReferenceFrame('N')

# Declare parameters
params = N.declare_parameters('rr rf lr ls lf l1 l2 l3 l4 mcd mef IC22 ICD11\
        ICD13 ICD33 ICD22 IEF11 IEF13 IEF33 IEF22 IF22 g')
# Declare coordinates and their time derivatives
q, qd = N.declare_coords('q', 8)
# Declare speeds and their time derivatives
u, ud = N.declare_speeds('u', 6)
# Unpack the lists
rr, rf, lr, ls, lf, l1, l2, l3, l4, mcd, mef, IC22, ICD11, ICD13, ICD33,\
        ICD22, IEF11, IEF13, IEF33, IEF22, IF22, g = params
q1, q2, q3, q4, q5, q6, q7, q8 = q
q1d, q2d, q3d, q4d, q5d, q6d, q7d, q8d = qd
u1, u2, u3, u4, u5, u6 = u
u1d, u2d, u3d, u4d, u5d, u6d = ud

# Create variables for to act as place holders in the vector from FO to FN
g31 = Function('g31')(t)
g33 = Function('g31')(t)
g31_s, g33_s = symbols('g31 g33')

# Some simplifying symbols / trig expressions
s1, s2, s3, s4, s5, s6, c1, c2, c3, c4, c5, c6, t2 = symbols('s1 \
        s2 s3 s4 s5 s6 c1 c2 c3 c4 c5 c6 t2')

symbol_subs_dict = {sin(q1): s1,
                    cos(q1): c1,
                    tan(q2): t2,
                    cos(q2): c2,
                    sin(q2): s2,
                    cos(q3): c3,
                    sin(q3): s3,
                    cos(q4): c4,
                    sin(q4): s4,
                    cos(q5): c5,
                    sin(q5): s5,
                    cos(q6): c6,
                    sin(q6): s6,
                    g31    : g31_s,
                    g33    : g33_s,
                    u1     : Symbol('u1'),
                    u2     : Symbol('u2'),
                    u3     : Symbol('u3'),
                    u4     : Symbol('u4'),
                    u5     : Symbol('u5'),
                    u6     : Symbol('u6')}

symbol_subs_dict_back = dict([(v,k) for k,v in symbol_subs_dict.items()])
trig_subs_dict = {c1**2: 1-s1**2,
                  c2**2: 1-s2**2,
                  c3**2: 1-s3**2,
                  c4**2: 1-s4**2,
                  c5**2: 1-s5**2,
                  c6**2: 1-s6**2}

# Reference Frames
# Yaw frame
A = N.rotate('A', 3, q1)
# Lean frame
B = A.rotate('B', 1, q2)
# Bicycle frame with rigidly attached rider
D = N.rotate('D', 'BODY312', (q1, q2, q3), I=(ICD11, ICD22, ICD33, 0, 0, ICD13))
# Rear wheel
C = D.rotate('C', 2, q4, I=(0, IC22, 0, 0, 0, 0), I_frame=D)
# Fork / handle bar assembly
E = D.rotate('E', 3, q5, I=(IEF11, IEF22, IEF33, 0, 0, IEF13))
# Front wheel
F = E.rotate('F', 2, q6, I=(0, IF22, 0, 0, 0, 0), I_frame=E)

g3_num = Vector(N[3] - dot(E[2], N[3])*E[2]).express(E)
g3_den = sqrt(g3_num.dict[E[1]]**2 + g3_num.dict[E[3]]**2)
g3 = Vector({E[1]: g3_num.dict[E[1]] / g3_den, E[3]: g3_num.dict[E[3]] / g3_den})
print g3

# Unit vector in the plane of the front wheel, pointed towards the ground
q9 = asin(dot(E[2], N[3]))
print 'Front assembly lean q9: ', q9
print dot(E[2], N[3])
q10 = asin(-dot(E[1], g3))
print 'Front assembly pitch q10: ', q10
print -dot(E[1], g3)

# Some manual manipulations to express g in the E frame without expanding the
# coefficients
g31_expr = rf*dot(g3, E[1])
g33_expr = rf*dot(g3, E[3])
print g31_expr
print g33_expr
num1 = g3_num.dict[E[1]]
num2 = g3_num.dict[E[3]]
den = g3_den
# Compute time derivatives using quotient rule and prevent unneccessary expansion
g31_expr_dt = (num1.diff(t)*den - num1*den.diff(t))/den**2
g33_expr_dt = (num2.diff(t)*den - num2*den.diff(t))/den**2
fo_fn_subs_dict = {g31: g31_expr, g33: g33_expr, g31.diff(t): g31_expr_dt,
        g33.diff(t): g33_expr_dt}
fo_fn = Vector({E[1]: g31, E[3]: g33})

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
EFO = FO.locate('EFO', l3*E[1] + l4*E[3], E, mass=mef)
# Locate front wheel contact point (fixed in the front wheel)
FN = FO.locate('FN', fo_fn, F)
# Locate another point fixed in N
N1 = CO.locate('N1', rr*B[3] - q7*N[1] - q8*N[2])

###############################################################################
#  Construction of kinematic differential equations ###########################
###############################################################################
# Form mapping from u's to qdots
# Definitions of the generalized speeds in terms of time derivatives of
# coordinates
u_rhs = [dot(D.ang_vel(N), D[1]),   # := u1
         dot(D.ang_vel(N), D[2]),   # := u2
         dot(D.ang_vel(N), D[3]),   # := u3
         dot(C.ang_vel(N), D[2]),   # := u4
         dot(E.ang_vel(N), D[3]),   # := u5
         dot(F.ang_vel(N), E[2])]   # := u6

qd_to_u = coefficient_matrix(u_rhs, qd[:-2])

# Steady turning conditions: q2d = q3d = q5d = 0 rad / sec
u_steady = qd_to_u * Matrix([q1d, 0, 0, q4d, 0, q6d])

u_to_qd = qd_to_u.inverse_ADJ().expand().subs(N.csqrd_dict).expand().\
        subs({sin(q2)/cos(q2):tan(q2)})

print u_to_qd
stop

# Form velocity of rear wheel center in two distinct but equivalent ways.  This
# allows for the rates of the rear wheel coordinates to be determined.
vco1 = dt(CO.rel(N1), N)
vco2 = cross(C.ang_vel(N), CO.rel(N.O))
eq1 = dot(vco1 - vco2, N[1]).expand().subs(N.csqrd_dict).expand()
eq2 = dot(vco1 - vco2, N[2]).expand().subs(N.csqrd_dict).expand()
xy_rates = solve([eq1, eq2], qd[6:])

# Create the kinematic differential equations as both a list of Sympy Eq
# objects, and as a dictionary whose keys are the qdots and whose values are
# the right hand sides of the kinematic differential equations.
kindiff_rhs = u_to_qd*Matrix(u)

kindiff_eqns = [Eq(qdot, qdot_rhs) for qdot, qdot_rhs in zip(qd[:6], kindiff_rhs)] +\
               [Eq(qd[6], xy_rates[q7d]),
                Eq(qd[7], xy_rates[q8d])]

# Depends upon yaw, lean, pitch, steer
func_params = (q1, q2, q4, q5)
ds = """\
Linear mapping from generalized speeds to time derivatives of coordinates.
"""
output_string = "from __future__ import division\n"
output_string += "from math import sin, cos, tan\n\n"

output_string = generate_function("kindiffs", kindiff_eqns, func_params,
        docstring=ds)
###############################################################################


###############################################################################
# Set velocity and angular velocities using only generalized speeds, and
# possibly generalized coordinates, but no time derivatives of coordinates
###############################################################################
# Rigid body angular velocities
# Angular velocity of rear frame with rigidly attached rider
D.abs_ang_vel = Vector(u1*D[1] + u2*D[2] + u3*D[3])
# Angular velocity of fork handlebar assembly
E.abs_ang_vel = Vector(u1*D[1] + u4*D[2] + u3*E[3]).express(E)
# Angular velocity of rear wheel
C.abs_ang_vel = Vector(u1*D[1] + u5*D[2] + u3*D[3])
# Angular velocity of front wheel
F.abs_ang_vel = E.abs_ang_vel - E.abs_ang_vel.dot(E[2])*E[2] + Vector(u6*E[2])

# Mass center velocities
CDO.abs_vel = express(cross(C.ang_vel(), CO.rel(N.O)) + cross(D.ang_vel(),\
                CDO.rel(CO)), D)
EFO.abs_vel = express(cross(F.ang_vel(), FO.rel(FN)) + cross(E.ang_vel(),\
                EFO.rel(FO)), E)
###############################################################################

###############################################################################
# Motion constraints
# We have introduced 6 generalized speeds to describe the angular velocity of
# the rigid bodies and the velocity of the mass centers.  Along with these
# generalized speeds, there are 3 motion constraints which reduces the number
# of degrees of freedom of the system to 3, consistent with what we know about
# the bicycle.
###############################################################################

# Motion constraints associated with point DE (top of steer axis)
# Velocity of point DE must be identical when formulated in the following two
# ways:
# Method 1:
vden1 = cross(C.ang_vel(), CO.rel(N.O).express(D)) +\
        cross(D.ang_vel(), DE.rel(CO))
# Method 2:
vden2 = cross(F.ang_vel(), FO.rel(FN)) +\
        cross(E.ang_vel(), DE.rel(FO))
# Form the constraint equations:
speed_constraint_eqs = [dot(vden1 - vden2, D[1]),\
                        dot(vden1 - vden2, D[2]),\
                        dot(vden1 - vden2, D[3])]

# Form the constraint matrix for B*u = 0
B_con = coefficient_matrix(speed_constraint_eqs, u)
B_con[0, 3] = B_con[0, 3].subs({cos(q5)**2: 1-sin(q5)**2}).expand()
B_con[1, 0] = B_con[1, 0].subs({sin(q5)**2: 1-cos(q5)**2}).expand()

# Use: u1 = dot(W_D_N>, D[1])
#      u4 = dot(W_E_N>, E[3])
#      u6 = dot(W_F_N>, 2[2])
# As the independent speeds
u_indep = [u1, u4, u6]
u_dep = [u2, u3, u5]

# The constraints can be written as:
# B*u = 0
# Bd * ud + Bi * ui = 0
# ud = -inv(Bd)*Bi*ui
# Form inv(Bd), Bi, and a substitution dictionary
Bd_inv, Bi, dep_ci, indep_ci, B_subs_dict = \
        transform_matrix(B_con, u, u_dep, subs_dict=True)
stop
B_subs_dict_time = {}
B_subs_dict_time_rev = {}
B_subs_dict_derivs = {}
B_subs_dict_dt_rhs = {}
for k, v in B_subs_dict.items():
    tv = Function("_"+str(k.name))(t)
    B_subs_dict_time[k] = tv
    B_subs_dict_time_rev[tv] = k
    tvd = Symbol("_"+str(k.name)+"d")
    B_subs_dict_derivs[tv.diff(t)] = tvd
    B_subs_dict_dt_rhs[tvd] = B_subs_dict[k].diff(t).subs({e1c.diff(t):
        Symbol('e1cd'), e3c.diff(t):Symbol('e3cd')})

# Matrix mapping inpependent speeds to dependent speeds
T_ud = -Bd_inv*Bi
T_ud_dt = zeros((3,3))
# Right hand sides of the dependent speeds
u_dep_rhs = T_ud*Matrix(u_indep)
for i, rhs in enumerate(u_dep_rhs):
    rhs = simplify(rhs.subs(N.symbol_dict)).subs(N.symbol_dict_back)
    n, d = rhs.as_numer_denom()
    n = collect(n, u)
    n_u0 = n.coeff(u0).subs(B_subs_dict_time)
    n_u3 = n.coeff(u3).subs(B_subs_dict_time)
    n_u5 = n.coeff(u5).subs(B_subs_dict_time)
    d_t = d.subs(B_subs_dict_time)
    # Quotient rule
    T_ud_dt[i, 0] = ((n_u0.diff(t)*d_t - n_u0*d_t.diff(t))/d_t**2)
    T_ud_dt[i, 1] = ((n_u3.diff(t)*d_t - n_u3*d_t.diff(t))/d_t**2)
    T_ud_dt[i, 2] = ((n_u5.diff(t)*d_t - n_u5*d_t.diff(t))/d_t**2)
    u_dep_rhs[i] = n / d

T_ud_dt = T_ud_dt.subs(B_subs_dict_derivs).subs(B_subs_dict_time_rev)

dep_speed_eqns = [Eq(u_d, u_d_r) for u_d, u_d_r in zip(u_dep, u_dep_rhs)]

# Parameters that the entries of matrix T_ud depend upon
func_params = params[:9] + (q1, q3, q4)

# Dealing with front wheel terms
fw_terms = {e1c_s: e1c_expr, e3c_s: e3c_expr}
for k,v in B_subs_dict.items():
    B_subs_dict[k] = v.subs({e1c: e1c_s, e3c: e3c_s})

# Nested terms
nt = [fw_terms, B_subs_dict]

# Docstring
ds = """\
Linear mapping from independent generalized speeds to dependent generalized
speeds.
"""
output_string += generate_function("speed_transform", dep_speed_eqns, u_indep,
        func_params, nested_terms=nt, docstring=ds)

for k,v in B_subs_dict.items():
    B_subs_dict[k] = v.subs({e1c_s: e1c, e3c_s: e3c})

# Set the speed transform matrix and determine the dependent speeds
T_con_dict, T_con_dt_dict = N.set_motion_constraint_matrix(T_ud, T_ud_dt,
        u_dep, u_indep, dep_ci, indep_ci)

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
fo_fn_subs = {e3c: e3c_s, e1c: e1c_s}
# Form Kane's equations and solve them for the udots
kanes_eqns = N.form_kanes_equations()
kanes_eqns_new = []
for i, ke in enumerate(kanes_eqns):
    lhs = ke.lhs.subs(fo_fn_subs)
    rhs = ke.rhs.subs(fo_fn_subs)
    kanes_eqns_new.append(Eq(lhs, rhs))

N.set_kanes_equations(kanes_eqns_new)

fw_terms[Symbol('e1cd')] = e1c_expr_dt
fw_terms[Symbol('e3cd')] = e3c_expr_dt
dyndiffs, mass_matrix_dict = N.solve_kanes_equations(dummy_mass_matrix=True)

for k, v in B_subs_dict.items():
    B_subs_dict[k] = v.subs(fo_fn_subs)
for k, v in B_subs_dict_dt_rhs.items():
    B_subs_dict_dt_rhs[k] = v.subs(fo_fn_subs)
for k, v in mass_matrix_dict.items():
    mass_matrix_dict[k] = v.subs(fo_fn_subs)

u_dep_dict = eqn_list_to_dict(dep_speed_eqns)
nt = [fw_terms, B_subs_dict, u_dep_dict, B_subs_dict_dt_rhs, T_con_dict, T_con_dt_dict,
        mass_matrix_dict]
ds = """\
Equations of motion for a benchmark bicycle.
"""
output_string += generate_function("eoms", kindiff_eqns+dyndiffs, q+u_indep, params,
        nested_terms=nt, docstring=ds, time=True)

output_string += 'mj_params = {"w": 1.02, "c": 0.08, "lambda": pi/10., "g":\
9.81, "rr": 0.3, "mr": 2.0, "IRxx": 0.0603, "IRyy": 0.12, "xb": 0.3, "zb":\
-0.9, "mb": 85.0, "IBxx": 9.2, "IBxz": 2.4, "IByy": 11.0, "IBzz": 2.8, "xh":\
0.9, "zh": -0.7, "mh": 4.0, "IHxx": 0.05892, "IHxz": -0.00756, "IHyy": 0.06,\
"IHzz": 0.00708, "rf": 0.35, "mf": 3.0, "IFxx": 0.1405, "IFyy": 0.28}\n'

fh = open('bicycle_lib.py', 'w')
fh.write(output_string)
fh.close()

