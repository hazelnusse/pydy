from pydy import *
from sympy import *
from bicycle_lib_hand import mj_params as mj_p

var('rr rrt rf rft xb zb xh zh mc md me mf mr mb mh w c lmbda')
var('xrb zrb xhf zhf')
var('l1 l2 l3 l4 lr lf ls')
var('IC22 IF11 IF22 ICD11 ICD13 ICD22 ICD33 IEF11 IEF13 IEF22 IEF33 IF11')
var('IRxx IRyy IBxx IBxz IByy IBzz IHxx IHxz IHyy IHzz IFxx IFyy')
# Subs dict and a collection list
sin_over_tan = {sin(lmbda)/tan(lmbda): cos(lmbda)}
col = [sin(lmbda), cos(lmbda)]

# Total mass of rear and front assemblies
mc = mr
md = mb
me = mh
mcd = mc + md
mef = me + mf
print 'mcd =',mcd
print 'mef =',mef

output_eqns = []
output_eqns.append(Eq(Symbol('mcd'), mcd))
output_eqns.append(Eq(Symbol('mef'), mef))

# Reference frames
N = NewtonianReferenceFrame('N')
# D Frame in the upright configuration
D = N.rotate('D', 2, lmbda)
# In the upright zero-steer configuration, the D and E frames are aligned
E = D

# Take N.O to be the rear contact point
# Locate all points using my parameters
CO = N.O.locate('CO', -(rr+rrt)*N[3], mass=mc)
CDO = CO.locate('CDO', l1*D[1] + l2*D[3], mass=md)
FO1 = CO.locate('FO1', lr*D[1] + ls*D[3] + lf*E[1], mass=mf)
EFO = FO1.locate('EFO', l3*E[1] + l4*E[3], mass=me)
FN1 = FO1.locate('FN1', (rf+rft)*N[3])

# Locate all points using Meijaard parameters
RO = N.O.locate('RO', -(rr + rrt)*N[3], mass=mr)
BO = N.O.locate('BO', xb*N[1] + zb*N[3], mass=mb)
HO = N.O.locate('HO', xh*N[1] + zh*N[3], mass=mh)
FO2 = N.O.locate('FO2', w*N[1] - (rf+rft)*N[3], mass=mf)
FN2 = N.O.locate('FN2', w*N[1])
CDO2 = N.O.locate('CDO2', xrb*N[1] + zrb*N[3])
EFO2 = N.O.locate('EFO2', xhf*N[1] + zhf*N[3])


# Convert geometric frame parameters
geom_eqs = [dot(FN1.rel(N.O)-FN2.rel(N.O), D[1]),
            dot(FN1.rel(N.O)-FN2.rel(N.O), D[3]),
            rf+rft-lf/sin(lmbda)-c/tan(lmbda)]       # From similar triangles
soln_geom = solve(geom_eqs, [lr, lf, ls])
print 'lr =', collect(soln_geom[lr].expand().subs(sin_over_tan), col)
print 'lf =', collect(soln_geom[lf].expand().subs(sin_over_tan), col)
print 'ls =', collect(soln_geom[ls], col)
for l in soln_geom:
    output_eqns.append(Eq(l, soln_geom[l]))

# Mass center's of Meijaards bodies
RBO = mass_center(N.O, [RO, BO])  # Relative to rear contact point
HFO = mass_center(FN2, [HO, FO2]) # Relative to front contact point

cm_eqs = [dot(CDO.rel(N.O) - RBO, D[1]),
          dot(CDO.rel(N.O) - RBO, D[3]),
          dot(EFO.rel(FN1) - HFO, D[1]),
          dot(EFO.rel(FN1) - HFO, D[3])]

soln_cm = {l1: -cm_eqs[0] + l1,
           l2: -cm_eqs[1] + l2,
           l3: -cm_eqs[2] + l3,
           l4: -cm_eqs[3] + l4}
print 'l1 =', soln_cm[l1]
print 'l2 =', soln_cm[l2]
print 'l3 =', soln_cm[l3]
print 'l4 =', soln_cm[l4]
for l in (l1, l2, l3, l4):
    output_eqns.append(Eq(l, soln_cm[l]))

# Test to ensure that the l1, l2, l3, l4 expressions just solved for result in
# the same mass center location as what is published by Sharp
xrbs = CDO.rel(N.O).dot(N[1])
zrbs = CDO.rel(N.O).dot(N[3])
xhfs = EFO.rel(N.O).dot(N[1])
zhfs = EFO.rel(N.O).dot(N[3])

# Below results match COM locations of Sharp 2008
print 'xrb =', xrbs.subs(soln_cm).subs(mj_p).n()
print 'zrb =', zrbs.subs(soln_cm).subs(mj_p).n()
print 'xhf =', xhfs.subs(soln_cm).subs(soln_geom).subs(mj_p).n()
print 'zhf =', zhfs.subs(soln_cm).subs(soln_geom).subs(mj_p).n()

# Inertia Calculations
# Using my parameters
I_CD_CDO = Inertia(D, [ICD11, ICD22, ICD33, 0, 0, ICD13])
I_EF_CDO = Inertia(E, [IEF11, IEF22, IEF33, 0, 0, IEF13])
# Using Meijaard's parameters
I_R_RO =  Inertia(N, [IRxx, IRyy, IRxx, 0, 0, 0])
I_B_BO =  Inertia(N, [IBxx, IByy, IBzz, 0, 0, IBxz])
I_H_HO =  Inertia(N, [IHxx, IHyy, IHzz, 0, 0, IHxz])
I_F_FO =  Inertia(N, [IFxx, IFyy, IFxx, 0, 0, 0])
# In plane Inertia of rear wheel
I_R_RO_p = Inertia(N, [IRxx, 0, IRxx, 0, 0, 0])
# In plane Inertia of front wheel
I_F_FO_p =  Inertia(N, [IFxx, 0, IFxx, 0, 0, 0])

# Inertia of rear frame with rigidly attached rider relative to rear assembly
# mass center
print BO_rel_CDO
print BO.rel(CDO2)
stop

#I_BO_CDO = inertia_of_point_mass(mb, BO.rel(CDO), N)
I_BO_CDO = inertia_of_point_mass(mb, BO.rel(CDO2), N)
I_B_CDO = I_B_BO + I_BO_CDO
# Planar inertia of rear wheel relative to rear assembly mass center
RO_rel_CDO = Vector(-xrb*N[1] + (-zrb-(rr+rft))*N[3])
#I_RO_CDO = inertia_of_point_mass(mr, RO.rel(CDO), N)
I_RO_CDO = inertia_of_point_mass(mr, RO_rel_CDO, N)
I_R_CDO = I_R_RO_p + I_RO_CDO
# Inertia of front frame relative to front assembly mass center
HO_rel_EFO = Vector((xh - xhf)*N[1] + (zh - zhf)*N[3])
#I_HO_EFO = inertia_of_point_mass(mh, HO.rel(EFO), N)
I_HO_EFO = inertia_of_point_mass(mh, HO_rel_EFO, N)
I_H_EFO = I_H_HO + I_HO_EFO
# Planar inertia of front wheel relative to front assembly mass center
FO_rel_CDO = Vector((w - xhf)*N[1] + (-(rf+rft) - zhf)*N[3])
#I_FO_EFO = inertia_of_point_mass(mf, FO2.rel(EFO), N)
I_FO_EFO = inertia_of_point_mass(mf, FO_rel_EFO, N)
I_F_EFO = I_F_FO_p + I_FO_EFO

# Rear assembly central inertia dyadics
I_RB_CDO = I_R_CDO + I_B_CDO
I_HF_EFO = I_H_EFO + I_F_EFO
# n2>*n2> components don't match Sharps, all others match to all digits.
#print 'I_RB_CDO =', I_RB_CDO.subs(soln_cm).subs(mj_p).n()
#print 'I_HF_EFO =', I_HF_EFO.subs(soln_cm).subs(soln_geom).subs(mj_p).n()
I_RB_CDO_D = I_RB_CDO.express(D)
I_HF_EFO_D = I_HF_EFO.express(D)
#print I_RB_CDO_D
#print I_HF_EFO_D
