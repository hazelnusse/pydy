from pydy import ReferenceFrame, cross, dot, dt, express, expression2vector, \
    coeff
from sympy import symbols, Function, S, solve, simplify, collect, sin, cos, pi

m, g, r1, r2, t = symbols("m, g r1 r2 t")
au1, au2, au3 = symbols("au1 au2 au3")
cf1, cf2, cf3 = symbols("cf1 cf2 cf3")
I, J = symbols("I J")
u3p, u4p, u5p = symbols("u3p u4p u5p")

q1 = Function("q1")(t)
q2 = Function("q2")(t)
q3 = Function("q3")(t)
q4 = Function("q4")(t)
q5 = Function("q5")(t)

def eval(a):
    subs_dict = {
        u3.diff(t): u3p,
        u4.diff(t): u4p,
        u5.diff(t): u5p,
        r1:1,
        r2:0,
        m:1,
        g:1,
        I:1,
        J:1,
        }
    r = a.subs(subs_dict)
    state = {
        q4: pi/6,
        u3: 1,
        u4: 0,
        u5: 2,
            }
    r = r.subs(state)
    return r

N = ReferenceFrame('N')
A = N.rotate("A", 3, q3)
B = A.rotate("B", 1, q4)
C = B.rotate("C", 2, q5)

#u3 = q3.diff(t)
u3 = Function("u3")(t)
u4 = Function("u4")(t)
u5 = Function("u5")(t)

P_NC_CO = r2*A[3] - r1*B[3]
#print "P_NC_CO> = ", P_NC_CO
P_NO_CO = q1*N[1] + q2*N[2] + P_NC_CO
#print "P_NO_CO> = ", P_NO_CO
A.set_omega(u3*A[3], N)
#print "W_A_N> = ", A.get_omega(N)
B.set_omega(A.get_omega(N) + u4*A[1], N)
#print "W_B_N> = ", B.get_omega(N)

C.set_omega(B.get_omega(N) + u5*B[2], N)
#print "W_C_N> = ", C.get_omega(N)

V_CN_N = au1*A[1] + au2*A[2] + au3*A[3]
#print "V_CN_N> = ", V_CN_N

V_CO_N = V_CN_N + cross(C.get_omega(N), P_NC_CO)
#print "V_CO_N> = ", V_CO_N

W3C = coeff(C.get_omega(N), u3)
#print "W3C> = ", W3C

W4C = coeff(C.get_omega(N), u4)
#print "W4C> = ", W4C

W5C = coeff(C.get_omega(N), u5)
#print "W5C> = ", W5C

WAU1C = coeff(C.get_omega(N), au1)
#print "WAU1C> = ", WAU1C

WAU2C = coeff(C.get_omega(N), au2)
#print "WAU2C> = ", WAU2C

WAU3C = coeff(C.get_omega(N), au3)
#print "WAU3C> = ", WAU3C

V3CO = coeff(V_CO_N, u3)
#print "V3CO> = ", V3CO

V4CO = coeff(V_CO_N, u4)
#print "V4CO> = ", V4CO

V5CO = coeff(V_CO_N, u5)
#print "V5CO> = ", V5CO

VAU1CO = coeff(V_CO_N, au1)
#print "VAU1CO> = ", VAU1CO

VAU2CO = coeff(V_CO_N, au2)
#print "VAU2CO> = ", VAU2CO

VAU3CO = coeff(V_CO_N, au3)
#print "VAU3CO> = ", VAU3CO

V3CN = coeff(V_CN_N, u3)
#print "V3CN> = ", V3CN

V4CN = coeff(V_CN_N, u4)
#print "V4CN> = ", V4CN

V5CN = coeff(V_CN_N, u5)
#print "V5CN> = ", V5CN

VAU1CN = coeff(V_CN_N, au1)
#print "VAU1CN> = ", VAU1CN

VAU2CN = coeff(V_CN_N, au2)
#print "VAU2CN> = ", VAU2CN

VAU3CN = coeff(V_CN_N, au3)
#print "VAU3CN> = ", VAU3CN

FORCE_CO = g*m*A[3]
#print "FORCE_CO> = ", FORCE_CO

FORCE_CN = cf1*A[1] + cf2*A[2] + cf3*A[3]
#print "FORCE_CN> = ", FORCE_CN

#print "Need to implement a dt() function to form acceleration of CO and \
        #        angular acceleration of C"

V_CO_N = V_CO_N.subs({au1: 0, au2: 0, au3: 0})
C.set_omega(C.get_omega(N).subs({au1: 0, au2: 0, au3: 0}), N, force=True)

A_CO_N = dt(V_CO_N, N, t)
A_CO_N = A_CO_N.subs({q3.diff(t): u3, q4.diff(t): u4, q5.diff(t): u5, \
        u3.diff(t): u3p, u4.diff(t): u4p, u5.diff(t): u5p})
#A_CO_N = express(A_CO_N, B)
#e = expression2vector(A_CO_N, B)
#print "b1> ", e[0]
#print "b2> ", e[1]
#print "b3> ", e[2]
#stop
#print "A_CO_N> = ", A_CO_N
#stop
ALF_C_N = dt(C.get_omega(N), N, t)
ALF_C_N = ALF_C_N.subs({q3.diff(t): u3, q4.diff(t): u4, q5.diff(t): u5, \
        u3.diff(t): u3p, u4.diff(t): u4p, u5.diff(t): u5p})
#ALF_C_N = express(ALF_C_N, B)
A_CO_N = A_CO_N.subs({q3.diff(t): u3, q4.diff(t): u4, q5.diff(t): u5})
#print "ALF_C_N> = ", ALF_C_N
#stop

#print A.get_omega(N)
#print dt(A.get_omega(N), N, t)
#print dt(B.get_omega(N), N, t)
#print dt(C.get_omega(N), N, t)
#print express(dt(C.get_omega(N), N, t), B)
#stop

RSTAR_C = -m*A_CO_N
#print "RSTAR_C> = ", RSTAR_C

I_C_CO = S(0)
I_alfa = dot(ALF_C_N, J*B[1])*B[1] + dot(ALF_C_N, I*B[2])*B[2] + \
        dot(ALF_C_N, J*B[3])*B[3]
#print "I_alfa", I_alfa
#stop

W_C_N = C.get_omega(N)
I_omega = dot(J*B[1], W_C_N)*B[1] + dot(I*B[2], W_C_N)*B[2] + \
        dot(J*B[3], W_C_N)*B[3]
#TSTAR_C = -dot(ALF_C_N, I_C_CO) - cross(W_C_N, dot(I_C_CO, W_C_N))
#print "cross(W_C_N, I_omega)*B[1]", dot(cross(W_C_N, I_omega),B[1])
#print "cross(W_C_N, I_omega)*B[2]", dot(cross(W_C_N, I_omega),B[2])
#print "cross(W_C_N, I_omega)*B[3]", dot(cross(W_C_N, I_omega),B[3])

#stop
TSTAR_C = -I_alfa - cross(W_C_N, I_omega)
#print "TSTAR_C> = ", TSTAR_C
#stop

F3 = dot(V3CO, FORCE_CO) + dot(V3CN, FORCE_CN)
F4 = dot(V4CO, FORCE_CO) + dot(V4CN, FORCE_CN)
F5 = dot(V5CO, FORCE_CO) + dot(V5CN, FORCE_CN)
F6 = dot(VAU1CO, FORCE_CO) + dot(VAU1CN, FORCE_CN)
F7 = dot(VAU2CO, FORCE_CO) + dot(VAU2CN, FORCE_CN)
F8 = dot(VAU3CO, FORCE_CO) + dot(VAU3CN, FORCE_CN)
#print "VAU3CO", VAU3CO
#print FORCE_CO
#print dot(VAU3CO, FORCE_CO)
#stop

#print "F3 = ", F3
#print "F4 = ", F4
#print "F5 = ", F5
#print "F6 = ", F6
#print "F7 = ", F7
#print "F8 = ", F8

#print "1", V3CO
#print "2", RSTAR_C
#print dot(V3CO, RSTAR_C)
#stop
RSTAR_C = RSTAR_C.expand()
F3S = dot(W3C, TSTAR_C) + dot(V3CO, RSTAR_C)
F4S = dot(W4C, TSTAR_C) + dot(V4CO, RSTAR_C)
F5S = dot(W5C, TSTAR_C) + dot(V5CO, RSTAR_C)
F6S = dot(WAU1C, TSTAR_C) + dot(VAU1CO, RSTAR_C)
F7S = dot(WAU2C, TSTAR_C) + dot(VAU2CO, RSTAR_C)
F8S = dot(WAU3C, TSTAR_C) + dot(VAU3CO, RSTAR_C)

#print "F3S = ", F3S
#print "F4S = ", F4S
#print "F5S = ", F5S
#print "F6S = ", F6S
#print "F7S = ", F7S
#print "F8S = ", F8S
#print eval(dot(W3C, TSTAR_C))
#print "-"*40
#print V3CO
#print RSTAR_C
#print eval(dot(V3CO, RSTAR_C))
#print eval(F3S)
#print eval(F4S)
#print eval(F5S)
#print eval(F6S)
#print eval(F7S)
#print eval(F8S)
#stop

EOM1 = F3 + F3S
EOM2 = F4 + F4S
EOM3 = F5 + F5S
EOM4 = F6 + F6S
EOM5 = F7 + F7S
EOM6 = F8 + F8S


#print "-"*80
#print "eom1 = ", EOM1
#print "eom2 = ", EOM2
#print "eom3 = ", EOM3
#print "eom4 = ", EOM4
#print "eom5 = ", EOM5
#print "eom6 = ", EOM6

#print "-"*80
#subs_dict = {
#    u3.diff(t): u3p,
#    u4.diff(t): u4p,
#    u5.diff(t): u5p,
#    r1:1,
#    r2:0,
#    m:1,
#    g:1,
#    I:1,
#    J:1,
#    }
#e1 = EOM1.subs(subs_dict)
#e2 = EOM2.subs(subs_dict)
#e3 = EOM3.subs(subs_dict)
#print "-"*80
#print eval(e1)
#print eval(e2)
#print eval(e3)

print EOM1.coeff(u3p)

stop
r = solve([EOM1, EOM2, EOM3], [u3p, u4p, u5p])
#print r
#stop
print "u3p = ", r[u3p]
print "u4p = ", r[u4p]
print "u5p = ", r[u5p]
