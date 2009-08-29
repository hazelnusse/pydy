from sympy import collect
from pydy import *

N = NewtonianReferenceFrame('N')

rr, rf, lr, ls, lf, l1, l2, l3, l4, mc, md, me, mf, IC11, IC22, ID11, ID13, ID33, ID22, IE11, IE13, IE33, IE22, IF11, IF22, g = N.declare_parameters('rr rf lr ls lf l1 l2 l3 l4 mc md me mf IC11 IC22 ID11 ID13 ID33 ID22 IE11 IE13 IE33 IE22 IF11 IF22 g')
(q1, q2, q3, q4, q5, q6), q_list, qdot_list = N.declare_coords('q', 6)
(u1, u2, u3, u4, u5, u6), u_list, udot_list = N.declare_speeds('u', 6)

# Reference Frames
# Yaw frame
A = N.rotate('A', 3, q1)
# Lean frame
#B = A.rotate('B', 1, q2)
# Bicycle frame pitch frame
#D = B.rotate('D', 2, q4, I=(ID11, ID22, ID33, 0, 0, ID13))
#print 'B[3].in D=',B[3].express(D)
#stop
D = N.rotate('D', 'BODY312', (q1, q2, q4), I=(ID11, ID22, ID33, 0, 0, ID13))
# Rear wheel spin frame
#C = B.rotate('C', 2, q3, I=(IC11, IC22, IC11, 0, 0, 0), I_frame=D)
C = N.rotate('C', 'BODY312', (q1, q2, q3), I=(IC11, IC22, IC11, 0, 0, 0), I_frame=D)
# Steer frame
E = D.rotate('E', 3, q5, I=(IE11, IE22, IE33, 0, 0, IE13))
# Front wheel spin frame
F = E.rotate('F', 2, q6, I=(IF11, IF22, IF11, 0, 0, 0), I_frame=E)

# In last reference frame, use E[2] instead of F[2] for the angular velocity,
# this prevents the ignorable coordinate q8 from appearing in the nonholonomic
# constraint equations.
#C._wrel = Vector(q3.diff(t)*D[2])
F._wrel = Vector(q6.diff(t)*E[2])

# Unit vector in the plane of the front wheel, pointed towards the ground
fo_fn = Vector(N[3] - (dot(E[2], N[3]))*E[2]).normalized

# Some manual manipulations to express g in the E frame without expanding the
# coefficients
N3inE = N[3].express(E)
e1c = rf*N3inE.dict[E[1]]*fo_fn.dict[N[3]]
e3c = rf*N3inE.dict[E[3]]*fo_fn.dict[N[3]]

# Necessary to prevent the denominatory from getting expanded and making things
# messy
num1, den1 = (rf*N3inE.dict[E[1]]*fo_fn.dict[N[3]]).as_numer_denom()
num2, den2 = (rf*N3inE.dict[E[3]]*fo_fn.dict[N[3]]).as_numer_denom()
assert den1==den2
den = Function('den')(t)
fo_fn = Vector({E[1]:num1/den, E[3]:num2/den})
den_subs_dict = {den: den1}

# Locate rear wheel center relative to point fixed in N, coincident with rear
# wheel contact
#CO = N.O.locate('CO', - rr*B[3], C, mass=mc)
CO = N.O.locate('CO', rr*sin(q4)*D[1] - rr*cos(q4)*D[3], C, mass=mc)
# Locate mass center of ricycle with rigidly attached rider
DO = CO.locate('DO', l1*D[1] + l2*D[3], D, mass=md)
# Locate top of steer axis
DE = CO.locate('DE', lr*D[1], D)
# Locate mass center of fork/handlebar assembly
EO = DE.locate('EO', l3*E[1] + l4*E[3], E, mass=me)
# Locate front wheel center
FO = DE.locate('FO', lf*E[1] + ls*E[3], E, mass=mf)
# Locate front wheel contact point (fixed in the front wheel)
FN = FO.locate('FN', fo_fn, F)

# Definitions of the generalized speeds
u_defs = N.define_speeds(
        [Eq(u_list[i-1], dot(D.ang_vel(), D[i])) for i in (1, 2, 3)] + \
        [Eq(u_list[3], dot(C.ang_vel(N), D[2])),
         Eq(u_list[4], dot(E.ang_vel(N), D[3])),
         Eq(u_list[5], dot(F.ang_vel(N), E[2]))])

# Solve u_defs for the qdots in terms of the u's
T, Tinv, kindiffs = N.form_kindiffs(u_defs, qdot_list)

print 'Kinematic differential equations'
for qd in qdot_list:
    print qd, '=', kindiffs[qd]

#A._wrel = A._wrel.express(D).subs(kindiffs)
#B._wrel = B._wrel.express(D).subs(kindiffs)
C._wrel =\
    C._wrel.express(D).subs(kindiffs).expandv().subs(N.csqrd_dict).expandv()
D._wrel =\
        D._wrel.express(D).subs(kindiffs).expandv().subs(N.csqrd_dict).expandv()
E._wrel = E._wrel.express(E).subs(kindiffs).expandv().subs(N.csqrd_dict).expandv()
F._wrel = F._wrel.express(E).subs(kindiffs).expandv().subs(N.csqrd_dict).expandv()

C.abs_ang_vel = C.ang_vel().express(D).expandv().subs(N.csqrd_dict).expandv()
D.abs_ang_vel = D.ang_vel().express(D).expandv().subs(N.csqrd_dict).expandv()
E.abs_ang_vel = E.ang_vel().express(E).expandv()
F.abs_ang_vel = F.ang_vel().express(E).expandv()

CO._vrel = cross(C.ang_vel(), CO.rel(N.O))
DO._vrel = cross(D.ang_vel(), DO.rel(CO))

DE._vrel = cross(D.ang_vel(), DE.rel(CO)).express(E)
EO._vrel = cross(E.ang_vel(), EO.rel(DE)).express(E)
FO._vrel = cross(E.ang_vel(), FO.rel(DE)).express(E)
FN._vrel = cross(F.ang_vel(), FN.rel(FO)).express(E)

# Ensuring all velocities expressions are in the most 'convenient' frame so
# that when accelerations are formed, things stay somewhat clean
CO.abs_vel = CO._vrel  # In D frame
DO.abs_vel = CO._vrel + DO._vrel  # In D frame
DE.abs_vel = CO._vrel.express(E) + DE._vrel  # In E frame
EO.abs_vel = DE.abs_vel + EO._vrel # In E frame
FO.abs_vel = DE.abs_vel + FO._vrel # In E frame
FN.abs_vel = FO.abs_vel + FN._vrel # In E frame

# Setting the angular accelerations and accelerations
C.abs_ang_acc = dt(C.ang_vel(), D) + cross(D.ang_vel(), C.ang_vel())
D.abs_ang_acc = dt(D.ang_vel(), D)
E.abs_ang_acc = dt(E.ang_vel(), E).subs(kindiffs).expandv()
F.abs_ang_acc = (dt(F.ang_vel(), E).subs(kindiffs) + cross(E.ang_vel(),
    F.ang_vel())).expandv().subs(N.csqrd_dict).expandv()

CO.abs_acc = (dt(CO.vel(), D) + cross(D.ang_vel(),
    CO.vel())).subs(kindiffs).expandv()#.subs(N.csqrd_dict).expandv()
DO.abs_acc = (dt(DO.vel(), D) + cross(D.ang_vel(),
    DO.vel())).subs(kindiffs).expandv()#.subs(N.csqrd_dict).expandv()
DE.abs_acc = (dt(DE.vel(), E).subs(kindiffs) + cross(E.ang_vel(),
    DE.vel())).subs(kindiffs).expandv()#.subs(N.csqrd_dict).expandv()
EO.abs_acc = (dt(EO.vel(), E).subs(kindiffs) + cross(E.ang_vel(),
    EO.vel())).subs(kindiffs).expandv().subs({cos(q5)**2: 1-sin(q5)**2}).expandv()
FO.abs_acc = (dt(FO.vel(), E).subs(kindiffs) + cross(E.ang_vel(),
    FO.vel())).subs(kindiffs).expandv().subs({cos(q5)**2: 1-sin(q5)**2}).expandv()
FN.abs_acc = (dt(FN.vel(), E).subs(kindiffs) + cross(E.ang_vel(),
    FN.vel())).subs(kindiffs).expandv().subs({cos(q5)**2: 1-sin(q5)**2}).expandv()
"""
print 'Absolute angular velocities'
print C.ang_vel()
print D.ang_vel()
print E.ang_vel()
print F.ang_vel()
print 'Absolute velocities'
print CO.vel()
print DO.vel()
print DE.vel()
print EO.vel()
print FO.vel()
print FN.vel()
print 'Absolute angular accelerations'
print C.ang_acc()
print D.ang_acc()
print E.ang_acc()
print F.ang_acc()
print 'Absolute accelerations'
print 'CO', CO.acc()
print 'DO', DO.acc()
print 'DE', DE.acc()
print 'EO', EO.acc()
print 'FO', FO.acc()
print 'FN', FN.acc()
stop
"""

constraint_eqs = [Eq(dot(FN.vel(), A[i]).expand().subs({cos(q1)**2:
    1-sin(q1)**2, cos(q5)**2:1-sin(q5)**2, cos(q4)**2:1-sin(q4)**2}).expand(), 0) for i in (1,2,3)]

# Form the constraint matrix for B*u = 0
B_con = N.form_constraint_matrix(constraint_eqs, u_list)
# Hand simplifications
B_con[0,0] = B_con[0,0].subs({sin(q4)**2:1-cos(q4)**2}).expand()
B_con[1,5] = B_con[1,5].subs({cos(q2)**2:1-sin(q2)**2}).expand()
B_con[2,0] = B_con[2,0].subs({sin(q5)**2:1-cos(q5)**2}).expand()

N.set_constraint_matrix(B_con)

# Form speed transform matrix for ud = (adj/det)*ui = T*ui
adj, det = N.form_speed_transform_matrix(B_con, [u1, u3, u4])

# Hand simplifications on the adjugate matrix
adj[1,0] = adj[1,0].subs({sin(q4)**2:1-cos(q4)**2}).expand()
adj[1,2] = adj[1,2].subs({sin(q4)**2:1-cos(q4)**2}).expand()

# Set the speed transform matrix and determine the dependent speeds
ud = N.set_speed_transform_matrix(adj, det)

# Set the kindiffs and dependent_speeds
N.setkindiffs(kindiffs, ud, acc=False)

# Apply gravity
N.gravity(g*N[3])

# Form Kane's equations and solve them for the udots
kanes_eqns = N.form_kanes_equations()
stop

dyndiffs = N.solve_kanes_equations()
