from pydy import ReferenceFrame, cross, dot
from sympy import symbols, Function, S

m, g, r1, r2, t = symbols("m, g r1 r2 t")
au1, au2, au3 = symbols("au1 au2 au3")
cf1, cf2, cf3 = symbols("cf1 cf2 cf3")
I, J = symbols("I J")

q1 = Function("q1")(t)
q2 = Function("q2")(t)
q3 = Function("q3")(t)
q4 = Function("q4")(t)
q5 = Function("q5")(t)

N = ReferenceFrame('N')
A = N.rotate("A", 3, q3)
B = A.rotate("B", 1, q4)
C = B.rotate("C", 2, q5)

#u3 = q3.diff(t)
u3 = Function("u3")(t)
u4 = Function("u4")(t)
u5 = Function("u5")(t)

P_NC_CO = r2*A[3] - r1*B[3]
print "P_NC_CO> = ", P_NC_CO
P_NO_CO = q1*N[1] + q2*N[2] + P_NC_CO
print "P_NO_CO> = ", P_NO_CO
W_A_N = u3*A[3]
print "W_A_N> = ", W_A_N
W_B_N = W_A_N + u4*A[1]
print "W_B_N> = ", W_B_N

W_C_N = W_B_N + u5*B[2]
print "W_C_N> = ", W_C_N

V_CN_N = au1*A[1] + au2*A[2] + au3*A[3]
print "V_CN_N> = ", V_CN_N

V_CO_N = V_CN_N + cross(W_C_N, P_NC_CO)
print "V_CO_N> = ", V_CO_N

W3C = W_C_N.coeff(u3)
print "W3C> = ", W3C

W4C = W_C_N.coeff(u4)
print "W4C> = ", W4C

W5C = W_C_N.coeff(u5)
print "W5C> = ", W5C


V3CO = V_CO_N.coeff(u3)
print "V3CO> = ", V3CO

V4CO = V_CO_N.coeff(u4)
print "V4CO> = ", V4CO

V5CO = V_CO_N.coeff(u5)
print "V5CO> = ", V5CO


V3CN = V_CN_N.coeff(au1)
print "V3CN> = ", V3CN

V4CN = V_CN_N.coeff(au2)
print "V4CN> = ", V4CN

V5CN = V_CN_N.coeff(au3)
print "V5CN> = ", V5CN

FORCE_CO = g*m*A[3]
print "FORCE_CO> = ", FORCE_CO

FORCE_CN = cf1*A[1] + cf2*A[2] + cf3*A[3]
print "FORCE_CN> = ", FORCE_CN

print "Need to implement a dt() function to form acceleration of CO and \
        angular acceleration of C"
A_CO_N = S(0)
print "A_CO_N> = ", A_CO_N
ALF_C_N = S(0)
print "ALF_C_N> = ", ALF_C_N

RSTAR_C = -m*A_CO_N
print "RSTAR_C> = ", RSTAR_C

I_C_CO = S(0)
TSTAR_C = -dot(ALF_C_N, I_C_CO) - cross(W_C_N, dot(I_C_CO, W_C_N))
print "TSTAR_C> = ", TSTAR_C


