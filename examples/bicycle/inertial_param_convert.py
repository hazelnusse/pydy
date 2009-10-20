from pydy import *
from sympy import *
from bicycle_lib_hand import mj_params as mj_p

# Reference frames
N = NewtonianReferenceFrame('N')
(q,), (qd,) = N.declare_coords('q', 1)
lmbda=q


var('rr rrt rf rft xb zb xh zh mc md me mf mr mb mh w c')
var('xrb zrb xhf zhf')
var('l1 l2 l3 l4 lr lf ls')
var('IC22 IF11 IF22 ICD11 ICD13 ICD22 ICD33 IEF11 IEF13 IEF22 IEF33 IF11')
var('IRxx IRyy IBxx IBxz IByy IBzz IHxx IHxz IHyy IHzz IFxx IFyy')
# Subs dict and a collection list
sin_over_tan = {sin(lmbda)/tan(lmbda): cos(lmbda)}
col = [sin(lmbda), cos(lmbda)]

output_eqns = [Eq(rr, rr),
               Eq(rrt, rrt),
               Eq(rf, rf),
               Eq(rft, rft)]
output_eqns = []

# Total mass of rear and front assemblies
mc = mr
md = mb
me = mh

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
RO = N.O.locate('RO', -(rr + rrt)*N[3], mass=mr)  # Same point as CO
BO = N.O.locate('BO', xb*N[1] + zb*N[3], mass=mb)
HO = N.O.locate('HO', xh*N[1] + zh*N[3], mass=mh)
FO2 = N.O.locate('FO2', w*N[1] - (rf+rft)*N[3], mass=mf) # Same point as FO
FN2 = N.O.locate('FN2', w*N[1])
CDO2 = N.O.locate('CDO2', xrb*N[1] + zrb*N[3])
EFO2 = N.O.locate('EFO2', xhf*N[1] + zhf*N[3])


# Convert geometric frame parameters
geom_eqs = [dot(FN1.rel(N.O)-FN2.rel(N.O), D[1]),
            dot(FN1.rel(N.O)-FN2.rel(N.O), D[3]),
            rf+rft-lf/sin(lmbda)-c/tan(lmbda)]       # From similar triangles
soln_geom = solve(geom_eqs, [lr, lf, ls])
for l in (lr, ls, lf):
    output_eqns.append(Eq(l,
        collect(soln_geom[l].expand().subs(sin_over_tan), col)))

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
for l in (l1, l2, l3, l4):
    output_eqns.append(Eq(l, soln_cm[l]))

# Test to ensure that the l1, l2, l3, l4 expressions just solved for result in
# the same mass center location as what is published by Sharp
xrbs = CDO.rel(N.O).dot(N[1])
zrbs = CDO.rel(N.O).dot(N[3])
xhfs = EFO.rel(N.O).dot(N[1])
zhfs = EFO.rel(N.O).dot(N[3])

# Below results match COM locations of Sharp 2008
#print 'xrb =', xrbs.subs(soln_cm).subs(mj_p).n()
#print 'zrb =', zrbs.subs(soln_cm).subs(mj_p).n()
#print 'xhf =', xhfs.subs(soln_cm).subs(soln_geom).subs(mj_p).n()
#print 'zhf =', zhfs.subs(soln_cm).subs(soln_geom).subs(mj_p).n()


# Distances from rear assembly mass center to rear wheel center
l1c = RO.rel(CDO).dot(D[1])
l3c = RO.rel(CDO).dot(D[3])
# Distances from rear assembly mass center to rear frame and rider mass center
l1d = BO.rel(CDO).dot(D[1])
l3d = BO.rel(CDO).dot(D[3])

# Distance from front assembly mass center to front fork/handle bar
l1e = HO.rel(EFO).dot(D[1])
l3e = HO.rel(EFO).dot(D[3])
# Distance from front assembly mass center to front wheel mass center
l1f = FO1.rel(EFO).dot(D[1])
l3f = FO1.rel(EFO).dot(D[3])

nt = {Symbol('l1c'): l1c,
                Symbol('l3c'): l3c,
                Symbol('l1d'): l1d,
                Symbol('l3d'): l3d,
                Symbol('l1e'): l1e,
                Symbol('l3e'): l3e,
                Symbol('l1f'): l1f,
                Symbol('l3f'): l3f}


# Using Meijaard's parameters
I_C_CO =  Inertia(N, [IRxx, IRyy, IRxx, 0, 0, 0])
I_D_DO =  Inertia(N, [IBxx, IByy, IBzz, 0, 0, IBxz])
I_E_EO =  Inertia(N, [IHxx, IHyy, IHzz, 0, 0, IHxz])
I_F_FO =  Inertia(N, [IFxx, IFyy, IFxx, 0, 0, 0])
# In plane Inertia of rear wheel
I_C_CO_p = Inertia(D, [IRxx, 0, IRxx, 0, 0, 0])
# In plane Inertia of front wheel
I_F_FO_p =  Inertia(E, [IFxx, 0, IFxx, 0, 0, 0])


l1c, l3c, l1d, l3d, l1e, l3e, l1f, l3f = symbols('l1c l3c l1d l3d l1e l3e l1f\
        l3f')

# Position from CDO to CO
CO_rel_CDO = Vector(l1c*D[1] + l3c*D[3])
# Position from CDO to DO
DO_rel_CDO = Vector(l1d*D[1] + l3d*D[3])
# Position from EFO to EO
EO_rel_EFO = Vector(l1e*E[1] + l3e*E[3])
# Position from EFO to FO
FO_rel_EFO = Vector(l1f*E[1] + l3f*E[3])

# Inertia of a particle, of mass mc, relative to CDO
I_CO_CDO = inertia_of_point_mass(mc, CO_rel_CDO, D)
# Parallel axis theorem for rear wheel, except out of plane inertia of wheel
I_C_CDO = I_C_CO_p + I_CO_CDO

# Inertia of a particle of mass md relative to the rear assembly mass center
I_DO_CDO = inertia_of_point_mass(md, DO_rel_CDO, D)
# Parallel axis theorem for rider
I_D_CDO = I_D_DO.express(D) + I_DO_CDO

I_CD_CDO = I_C_CDO + I_D_CDO

# Inertia of a particle, of mass me, relative to EFO
I_EO_EFO = inertia_of_point_mass(me, EO_rel_EFO, E)
# Parallel axis theorem for fork handlebar assembly
I_E_EFO = I_E_EO.express(E) + I_EO_EFO

# Inertia of a particle of mass md relative to the rear assembly mass center
I_FO_EFO = inertia_of_point_mass(mf, FO_rel_EFO, E)
# Parallel axis theorem for rider
I_F_EFO = I_F_FO_p + I_FO_EFO
I_EF_EFO = I_E_EFO + I_F_EFO

mcd = mc + md
mef = me + mf
output_eqns.append(Eq(Symbol('mcd'), mcd))
output_eqns.append(Eq(Symbol('mef'), mef))
output_eqns.append(Eq(Symbol('IC22'), IRyy))


ICD11 = dot(D[1], dot(I_CD_CDO, D[1]))
output_eqns.append(Eq(Symbol('ICD11'), ICD11))
ICD22 = dot(D[2], dot(I_CD_CDO, D[2]))
ICD13 = dot(D[1], dot(I_CD_CDO, D[3]))
output_eqns.append(Eq(Symbol('ICD13'), ICD13))
output_eqns.append(Eq(Symbol('ICD22'), ICD22))
ICD33 = dot(D[3], dot(I_CD_CDO, D[3]))
output_eqns.append(Eq(Symbol('ICD33'), ICD33))
ICD13 = dot(D[1], dot(I_CD_CDO, D[3]))
output_eqns.append(Eq(Symbol('ICD13'), ICD13))

IEF11 = dot(D[1], dot(I_EF_EFO, D[1]))
output_eqns.append(Eq(Symbol('IEF11'), IEF11))
IEF22 = dot(D[2], dot(I_EF_EFO, D[2]))
output_eqns.append(Eq(Symbol('IEF22'), IEF22))
IEF33 = dot(D[3], dot(I_EF_EFO, D[3]))
output_eqns.append(Eq(Symbol('IEF33'), IEF33))
IEF13 = dot(D[1], dot(I_EF_EFO, D[3]))
output_eqns.append(Eq(Symbol('IEF13'), IEF13))

output_eqns.append(Eq(Symbol('IF22'), IFyy))
ops = 0
for e in output_eqns:
    print e
    ops += e.rhs.count_ops()
print ops


params = N.declare_parameters('rr rrt rf rft lr ls lf l1 l2 l3 l4 mcd mef IC22\
        ICD11 ICD22 ICD33 ICD13 IEF11 IEF22 IEF33 IEF13 IF22 g')

input = [w, c, lmbda, rr, rrt, rf, rft, xb, zb, xh, zh, mr, mb, mh, mf, IRxx, IRyy, IBxx, IByy, IBzz, IBxz, IHxx, IHyy, IHzz,
        IHxz, IFxx, IFyy]

output_string = "from __future__ import division\n"
output_string += "from math import sin, cos\n\n"

output_string += generate_function("convert_params", output_eqns, input,
        nested_terms=[nt])
print output_string
stop


file = open('convert_parameters.py', 'w')
file.write(output_string)
file.close()
