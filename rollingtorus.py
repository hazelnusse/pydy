from math import sin, cos, pi

from numpy import array, arange
from sympy import symbols, Function, S, solve, simplify, \
        collect, Matrix, lambdify, trigsimp, expand, Eq, pretty_print

from pydy import ReferenceFrame, cross, dot, dt, express, \
        Vector, UnitVector, Dyad, Inertia, InertiaForce, \
        InertiaTorque

# Constants
m, g, r1, r2, t, I, J= symbols("m g r1 r2 t I J")

#I = (r1**2/2 + 5*r2**2/8)*m  # Central moment of inertia about any diameter
#J = (r1**2 + 3*r2**2/4)*m    # Central moment of inertia about normal axis

# Auxilliary generalized speeds and contact forces
au1, au2, au3, cf1, cf2, cf3 = symbols("au1 au2 au3 cf1 cf2 cf3")
u3p, u4p, u5p = symbols("u3p u4p u5p")


q1, q2, q3, q4, q5 = [Function("q%i"%j)(t) for j in range(1, 6)]
u1, u2, u3 = [Function("u%i"%j)(t) for j in (1,2,3)]

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

# Locate the mass center of torus
CO = N.O.locate('CO', q1*N[1] + q2*N[2] - r2*N[3] - r1*B[3])

# Contact point, taken to be fixed on C
CN = CO.locate('CN', r1*B[3] + r2*N[3], C)

# Nonholonomic rolling constraints, imposing the requirement that the velocity
# of the point fixed to C, instantaneously in contact with the ground, must be
# zero.
nh1, nh2 = dot(CN.vel, N[1]), dot(CN.vel, N[2])
dependent_rates = solve([nh1,nh2], (q1.diff(t), q2.diff(t)))
print 'Dependent rates', dependent_rates

# Imposing the constraints upon the velocity of the two points CO, CN.  Note
# that CN.vel == zero after inposing the constraints.
CN.vel = express(CN.vel.subs(dependent_rates), N)
CO.vel = CO.vel.subs(dependent_rates)

# Define lists for the time derivatives of the coordinates, the generalized
# speeds
qdot_list = [q3.diff(t), q4.diff(t), q5.diff(t)]
u_lhs = [Function('u%i'%i)(t) for i in (1,2,3)]

# Take the generalized speeds to be the B frame measure numbers of the angular
# velocity of the torus C in N
u_rhs = [dot(C.get_omega(N), B[i]) for i in (1,2,3)]

# Simplest definition of generalized speeds
#u_rhs = qdot_list

# Create the equations that define the generalized speeds, then solve them for
# the time derivatives of the generalized coordinates
u_definitions = [Eq(u_l, u_r) for u_l, u_r in zip(u_lhs, u_rhs)]
kindiffs = solve(u_definitions, qdot_list)

print 'Kinematic differential equations', kindiffs
# Replace occurances of q_dots with the corresponding definition from the
# kinematic differential equations.
CO.vel = CO.vel.subs(kindiffs)
C.W[N] = C.W[N].subs(kindiffs)
print 'v_co_n = ', CO.vel
print 'w_c_n = ', C.W[N]

# Compute the partial velocity of CO
partial_v = CO.vel.partials(u_lhs)
# Compute the partial angular velocity of C
partial_w = C.W[N].partials(u_lhs)

print 'Partial velocities of CO'
pretty_print([pydy_str(p) for p in partial_v])
print 'Partial angular velocities of C'
pretty_print([pydy_str(p) for p in partial_w])

CO.acc = dt(CO.vel, N).subs(kindiffs)
C.alpha = express(dt(C.W[N], N).subs(kindiffs), B)
print 'a_co_n =', CO.acc
print 'alf_c_n =', C.alpha
stop

# Force acting at CO, mass of CO
CO.force = Vector(m*g*N[3])
CO.mass = m

# Central inertia dyad of C
C.inertia = Inertia(B, (I, I, J, 0, 0, 0))

# Inertia Force and Inertia Torque
R_star = InertiaForce(m, CO.acc)
#print R_star
C.W[N] = express(C.W[N], B)
T_star = InertiaTorque(C.inertia, C.W[N], C.alpha)

#print T_star

# Generalized active forces
GAF = [dot(CO.force, p_v) for p_v in partial_v]

# Generalized inertia forces
GIF = [dot(p_w, T_star) + dot(p_v, R_star) for p_w, p_v in zip(partial_w,
    partial_v)]

EOMS = [AF + IF for AF, IF in zip(GAF, GIF)]

#print EOMS
#print trigsimp(EOMS[0].coeff(u1.diff(t)))
#print EOMS[1].coeff(u1.diff(t))
#print EOMS[2].coeff(u1.diff(t))
#print EOMS[0].coeff(u2.diff(t))
#print EOMS[1].coeff(u2.diff(t))
#print trigsimp(EOMS[2].coeff(u2.diff(t)))
#print EOMS[0].coeff(u3.diff(t))
#print trigsimp(EOMS[1].coeff(u3.diff(t)))
#print trigsimp(EOMS[2].coeff(u3.diff(t)))

dyn_eqns = solve(EOMS, u1.diff(t), u2.diff(t), u3.diff(t))
print dyn_eqns
stop






P_NC_CO = r2*A[3] - r1*B[3]
#print "P_NC_CO> = ", P_NC_CO
P_NO_CO = q1*N[1] + q2*N[2] + P_NC_CO

V_CN_N = au1*A[1] + au2*A[2] + au3*A[3]
#print "V_CN_N> = ", V_CN_N

V_CO_N = V_CN_N + cross(C.get_omega(N), P_NC_CO)
#print "V_CO_N> = ", V_CO_N

qdots = [q3.diff(t), q4.diff(t), q5.diff(t), au1, au2, au3, u3.diff(t), \
        u4.diff(t), u5.diff(t)]
us = [u3, u4, u5, au1, au2, au3, u3p, u4p, u5p]
print qdots
print us
gen_speeds = dict(zip(qdots, us))
print gen_speeds
C.set_omega(C.get_omega(N).subs(gen_speeds), N, force = True)

V_CO_N = V_CO_N.subs(gen_speeds)
V_CN_N = V_CN_N.subs(gen_speeds)


WC = coeff(C.get_omega(N), us[:-3])
VC = coeff(V_CO_N, us[:-3])
VCN = coeff(V_CN_N, us[:-3])

#print WC
#print VC
#print VCN
#stop

FORCE_CO = g*m*A[3]
#print "FORCE_CO> = ", FORCE_CO

FORCE_CN = cf1*A[1] + cf2*A[2] + cf3*A[3]

# Set the auxilliary generalized speeds to zero before forming the 
# acceleration of CO in N  and angular acceleration of C in N.
V_CO_N = V_CO_N.subs({au1: 0, au2: 0, au3: 0})
C.set_omega(C.get_omega(N).subs({au1: 0, au2: 0, au3: 0}), N, force=True)

A_CO_N = dt(V_CO_N, N, t)
A_CO_N = A_CO_N.subs(gen_speeds)
A_CO_N = express(A_CO_N,B)

ALF_C_N = dt(C.get_omega(N), N, t)
ALF_C_N = ALF_C_N.subs(gen_speeds)
ALF_C_N = express(ALF_C_N,B)

#print "A_CO_N> * B[1]= ", dot(A_CO_N,B[1])
#print "A_CO_N> * B[2]= ", dot(A_CO_N,B[3])
#print "A_CO_N> * B[3]= ", dot(A_CO_N,B[3])

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
#print "RSTAR_C> * B[1]= ", dot(RSTAR_C,B[1])
#print "RSTAR_C> * B[2]= ", dot(RSTAR_C,B[2])
#print "RSTAR_C> * B[3]= ", dot(RSTAR_C,B[3])

#stop
I_C_CO = S(0)
I_alfa = dot(ALF_C_N, J*B[1])*B[1] + dot(ALF_C_N, I*B[2])*B[2] + \
        dot(ALF_C_N, J*B[3])*B[3]
#print "I_alfa", I_alfa
#print "I_alfa * B[1]", dot(I_alfa, B[1])
#print "I_alfa * B[2]", dot(I_alfa, B[2])
#print "I_alfa * B[3]", dot(I_alfa, B[3])
# above matches autolev results
#stop

W_C_N = C.get_omega(N)
#print "W_C_N = ", W_C_N
#stop
I_omega = dot(J*B[1], W_C_N)*B[1] + dot(I*B[2], W_C_N)*B[2] + \
        dot(J*B[3], W_C_N)*B[3]

#print "I_omega * B[1]: ", dot(I_omega,B[1])
#print "I_omega * B[2]: ", dot(I_omega,B[2])
#print "I_omega * B[3]: ", dot(I_omega,B[3])
# Above matches autlev code
#stop



#TSTAR_C = -dot(ALF_C_N, I_C_CO) - cross(W_C_N, dot(I_C_CO, W_C_N))
#print "cross(W_C_N, I_omega)*B[1]", dot(cross(W_C_N, I_omega),B[1])
#print "cross(W_C_N, I_omega)*B[2]", dot(cross(W_C_N, I_omega),B[2])
#print "cross(W_C_N, I_omega)*B[3]", dot(cross(W_C_N, I_omega),B[3])

#stop
#print "w x (I * w) * B[1]: ", dot(cross(W_C_N, I_omega), B[1])
#print "w x (I * w) * B[2]: ", dot(cross(W_C_N, I_omega), B[2])
#print "w x (I * w) * B[3]: ", dot(cross(W_C_N, I_omega), B[3])
# above code matches with Autolev
#stop

TSTAR_C = -I_alfa - cross(W_C_N, I_omega)
#print "TSTAR_C> = ", TSTAR_C
#stop

#print "VC = ", VC
#print "VCN = ", VCN
#stop

FR = Matrix([dot(ViCO, FORCE_CO) for ViCO in VC]) + \
        Matrix([dot(ViCN, FORCE_CN) for ViCN in VCN])

#print "FR = \n", FR
# above code matches autolev
#stop


#RSTAR_C = RSTAR_C.expand()

#print "dot(WC3>, TSTAR_C) = ",dot(WC[0],TSTAR_C)
#stop
#print dot(VC[0],TSTAR_C)
#stop

FRSTAR = Matrix([dot(WiC, TSTAR_C) for WiC in WC]) + \
        Matrix([dot(ViCO, RSTAR_C) for ViCO in VC])

#print "FRSTAR = \n", FRSTAR
#stop
#print "F3S = ", FRSTAR[0]
#stop
#print "F4S = ", FRSTAR[1]
#print "F5S = ", FRSTAR[2]
#print "F6S = ", FRSTAR[3]
#print "F7S = ", FRSTAR[4]
#print "F8S = ", FRSTAR[5]
#stop
#F3S = dot(W3C, TSTAR_C) + dot(V3CO, RSTAR_C)
#F4S = dot(W4C, TSTAR_C) + dot(V4CO, RSTAR_C)
#F5S = dot(W5C, TSTAR_C) + dot(V5CO, RSTAR_C)
#F6S = dot(WAU1C, TSTAR_C) + dot(VAU1CO, RSTAR_C)
#F7S = dot(WAU2C, TSTAR_C) + dot(VAU2CO, RSTAR_C)
#F8S = dot(WAU3C, TSTAR_C) + dot(VAU3CO, RSTAR_C)

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


ZERO = FR + FRSTAR
#EOM1 = F3 + F3S
#EOM2 = F4 + F4S
#EOM3 = F5 + F5S
#EOM4 = F6 + F6S
#EOM5 = F7 + F7S
#EOM6 = F8 + F8S


#print "-"*80
#print "eom1 = ", EOM1
#print "eom2 = ", EOM2
#print "eom3 = ", EOM3
#print "eom4 = ", EOM4
#print "eom5 = ", EOM5
#print "eom6 = ", EOM6

print "-"*80
subs_dict = {
    r1:1,
    r2:0,
    m:1,
    g:1,
    I:1,
    J:1,
    q4:pi/6,
    g:1,
    u3:1,
    u4:0,
    u5:2,
    }

#E = ZERO.subs(subs_dict)
#print E
#stop

#e1 = EOM1.subs(subs_dict)
#e2 = EOM2.subs(subs_dict)
#e3 = EOM3.subs(subs_dict)
#print "-"*80
#print eval(e1)
#print eval(e2)
#print eval(e3)

#print EOM1.coeff(u3p)

#stop
#print ZERO[0]
#print ZERO[1]
#print ZERO[2]
subs_dict = {
    r1:1,
    r2:0,
    m:1,
    g:1,
    I:1,
    J:1,
    }
eqs = [x.subs(subs_dict) for x in ZERO[:3]]
r = solve(eqs, [u3p, u4p, u5p])
print "Solution:"
print r
print "_"*80
print "Integrating the equations of motion"
def eval_sol(sol, y):
    r = sol.subs({
        q3: y[2],
        q4: y[3],
        q5: y[4],
        u3: y[5],
        u4: y[6],
        u5: y[7],
        })
    return lambdify(y, r)

y = symbols("q1 q2 q3 q4 q5 u3 u4 u5")
u3p = eval_sol(r[u3p], y)
u4p = eval_sol(r[u4p], y)
u5p = eval_sol(r[u5p], y)

from solver import rk4int
def derivs(y, x):
    q1, q2, q3, q4, q5, u3, u4, u5 = y
    r1 = 1.0
    r2 = 0.
    return array([
        -cos(q3)*(r1-r2*cos(q4))*u5, #q1'
        -sin(q3)*(r1-r2*cos(q4))*u5, #q2'
        u3, #q3'
        u4, #q4'
        u5, #q5'
        u3p(*y), #u3'
        u4p(*y), #u4'
        u5p(*y), #u5'
        ])

#y0 = (0, 0, 0, pi/10, 0, 1.0, -0.1, 3.)
y0 = (0, 0, 0, pi/10, 0, 1.0, -0.1, 3.)
#y0 = (0., 0., 0., 0., 0., 0., 0., 3.)
#y0 = (0., 0., 0., 0., 0., 0., 0.2, 3.)
h = 0.01
t = arange(0., 40, h)
print "start"
sol = rk4int(derivs, y0, t, h)
print "finished"
from pylab import plot, show, legend
plot(sol[:, 0], sol[:, 1], "-", lw=2, label="x-y")
#plot(t, sol[:, 0], "-", lw=2, label="q1 (x)")
#plot(t, sol[:, 1], "-", lw=2, label="q2 (y)")
#plot(t, sol[:, 2], "-", label="q3")
#plot(t, sol[:, 3], "-", label="q4")
#plot(t, sol[:, 4], "-", label="q5")
#plot(t, sol[:, 5], "--", label="u3")
#plot(t, sol[:, 6], "--", label="u4")
#plot(t, sol[:, 7], "--", label="u5")
legend()
show()
