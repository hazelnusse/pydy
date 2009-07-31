# This program calculates the full nonlinear equations of motion (EOM) of
# Garcia's simplest point foot walker model.  The complete paper can be found:
# http://ruina.tam.cornell.edu/research/topics/locomotion_and_robotics/papers/some_results_passive_dynamic.pdf


from sympy import *
from pydy import *

# Define Model Parameters
t,Mh,mf,g,L,gamma = symbols('t Mh mf g L gamma')
# Mh = hip mass, (kg)
# mf = foot mass, (kg)
# g = gravitational acceleration constant, (m/s/s)
# L = leg length
# gamma = ground slope, (rad)

# Define Variables (generalized coordinates, (GC))
th = Function('theta')(t)
ph = Function('phi')(t)

# Define generalized speeds (GS)
u1 = Function('u1')(t)
u2 = Function('u2')(t)

# Reference Frames
N = NewtonianReferenceFrame('N')
F = N.rotate('F', 3, -gamma)
A = F.rotate('A', 3, th)
B = A.rotate('B', 3, -ph)

P = N.O.locate('P', L*A[2], frame=A, mass=Mh)
Q = P.locate('Q', -L*B[2], frame=B, mass=mf)

# Form kinematical differential equations, i.e., a linear system relating the
# time derivatives of the GC to the GS
u_list = [u1, u2]
kindiffs = {th.diff(t): u1, ph.diff(t): u2}
N.setgenspeeds(u_list)
N.setkindiffs(kindiffs)


N.apply_gravitational_force(g*N[2])

print N.O._vrel
print N.O._partialvrel
print P._vrel
print P._partialvrel
print Q._vrel
print Q._partialvrel

stop


"""
print N.W
print F.W[N]        # omega of F in N
print A.W[F]        # omega of A in F
print B.W[A]        # omega of B in A
"""

"""
print '\n'
print P.vel
print express(P.vel,F)
"""


# Calculate the partial derivatives wrt the GS
partials_p = P.vel.partials(u_lhs)
partials_q = Q.vel.partials(u_lhs)
print 'v_p_1 = ', partials_p[0], '\nv_p_2 = ', partials_p[1]
print 'v_q_1 = ', partials_q[0], '\nv_q_2 = ', partials_q[1]

# Once trigsimp() is fixed, this line will be unnecessary
partials_q[0] = Vector(-L*A[1] + L*B[1])

# Defining mass and gravitational forces acting on the particles
P.mass = Mh
P.force = Vector(-P.mass*g*N[2])
Q.mass = mf
Q.force = Vector(-Q.mass*g*N[2])

# Generalized Inertial Forces, (GIF)
R_star_p = InertiaForce(P.mass, P.acc)
R_star_q = InertiaForce(Q.mass, Q.acc)

# Generalized Active Forces, (GAF)
GAF = [dot(P.force, partials_p[i]) + dot(Q.force, partials_q[i]) for i in
        (0,1)]

print 'F1\n', GAF[0]
print 'F2\n', GAF[1]

GIF = [dot(R_star_p, partials_p[i]) + dot(R_star_q, partials_q[i]) for i in
        (0,1)]

# Equations of Motion = GAF + GIF = Fr + Frstar
EOMS = [AF + IF for AF, IF in zip(GAF, GIF)]
print 'EOM1: \n', EOMS[0]
print 'EOM2: \n', EOMS[1]
