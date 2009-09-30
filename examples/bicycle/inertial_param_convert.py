from pydy import *
from sympy import *

var('rr rrt rf rft xb zb xh zh mc md me mf w c lmbda')
var('l1 l2 l3 l4')
var('lr lf ls')
var('IC11 IC22 IF11 IF22 ID11 ID13 ID22 ID33 IE11 IE13 IE22 IE33')
var('IBxx IBxz IByy IBzz IHxx IHxz IHyy IHzz')
# Subs dict and a collection list
sin_over_tan = {sin(lmbda)/tan(lmbda): cos(lmbda)}
col = [sin(lmbda), cos(lmbda)]

# Reference frames
N = ReferenceFrame('N')
# D Frame in the upright configuration
D = N.rotate('D', 2, lmbda)
# In the upright zero-steer configuration, the D and E frames are aligned
E = D

# Total mass
mcd = mc + md
mef = me + mf

# Using my parameters
v_no_co = Vector(-(rr+rrt)*N[3])
v_no_do = v_no_co + Vector(l1*D[1] + l2*D[3])
v_fn_fo = Vector((rf + rft)*N[3])
v_fn_eo = Vector(l3*E[1] + l4*E[3])
v_no_fn = Vector(-(rr+rrt)*N[3] + lr*D[1] + ls*D[3] + lf*E[1] + (rf+rft)*N[3])
I_D = Dyad(ID11*D[1]*D[1] + ID22*D[2]*D[2] + ID13*D[1]*D[3] + ID13*D[3]*D[1] +
        ID33*D[3]*D[3])
I_E = Dyad(IE11*E[1]*E[1] + IE22*E[2]*E[2] + IE13*E[1]*E[3] + IE13*E[3]*E[1] +
        IE33*E[3]*E[3])

# Using Meijaard's parameters
v_no_co_mj = Vector(-(rr+rrt)*N[3])
v_no_do_mj = Vector(xb*N[1] + zb*N[3])
v_fn_fo_mj = Vector((rr + rrt)*N[3])
v_fn_eo_mj = Vector((xh-w)*N[1] + zh*N[3])
v_no_fn_mj = Vector(w*N[1])
I_B = Dyad(IBxx*N[1]*N[1] + IByy*N[2]*N[2] + IBxz*N[1]*N[3] + IBxz*N[3]*N[1] +
        IBzz*N[3]*N[3])
I_H = Dyad(IHxx*N[1]*N[1] + IHyy*N[2]*N[2] + IHxz*N[1]*N[3] + IHxz*N[3]*N[1] +
        IHzz*N[3]*N[3])

rear_cm_eqs = [dot(v_no_do-v_no_do_mj, D[1]),dot(v_no_do-v_no_do_mj, D[3])]
soln_rear = solve(rear_cm_eqs, [l1, l2])
print 'l1 =', collect(soln_rear[l1], col)
print 'l2 =', collect(soln_rear[l2], col)

front_cm_eqs = [dot(v_fn_eo-v_fn_eo_mj, D[1]),dot(v_fn_eo-v_fn_eo_mj, D[3])]
soln_front = solve(front_cm_eqs, [l3, l4])
print 'l3 =', collect(soln_front[l3], col)
print 'l4 =', collect(soln_front[l4], col)

geom_eqs = [dot(v_no_fn-v_no_fn_mj, D[1]),
            dot(v_no_fn-v_no_fn_mj, D[3]),
            rf+rft-lf/sin(lmbda)-c/tan(lmbda)]       # From similar triangles
soln_geom = solve(geom_eqs, [lr, lf, ls])
print 'lr =', collect(soln_geom[lr].expand().subs(sin_over_tan), col)
print 'lf =', collect(soln_geom[lf].expand().subs(sin_over_tan), col)
print 'ls =', collect(soln_geom[ls], col)

I_B_mat = zeros((3,3))
for i in (1,2,3):
    for j in (1,2,3):
        I_B_mat[i-1,j-1] = dot(D[j], dot(D[i], I_B)).expand()

soln_inertia_D = {ID11: I_B_mat[0,0], ID22: I_B_mat[1,1], ID33: I_B_mat[2,2],\
                  ID13: I_B_mat[0,2]}
print 'ID11 =', soln_inertia_D[ID11]
print 'ID13 =', soln_inertia_D[ID13]
print 'ID22 =', soln_inertia_D[ID22]
print 'ID33 =', soln_inertia_D[ID33]

I_H_mat = zeros((3,3))
for i in (1,2,3):
    for j in (1,2,3):
        I_H_mat[i-1,j-1] = dot(E[j], dot(E[i], I_H)).expand()
soln_inertia_E = {IE11: I_H_mat[0,0], IE22: I_H_mat[1,1], IE33: I_H_mat[2,2],\
                  IE13: I_H_mat[0,2]}
print 'IE11 =', soln_inertia_E[IE11]
print 'IE13 =', soln_inertia_E[IE13]
print 'IE22 =', soln_inertia_E[IE22]
print 'IE33 =', soln_inertia_E[IE33]

l1_t = md * l1 / mcd
l2_t = md * l2 / mcd
l3_t = me * l3 / mcd
l4_t = me * l4 / mcd

ICD11 = mc*l1_t**2 + IC11 + md*(l1 - l1_t)**2 + ID11
ICD33 = mc*l2_t**2 + IC11 + md*(l2 - l2_t)**2 + ID33
IEF11 = mf*l3_t**2 + IF11 + me*(l3 - l3_t)**2 + IE11
IEF33 = mf*l4_t**2 + IF11 + me*(l4 - l4_t)**2 + IE33

