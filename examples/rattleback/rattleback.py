from pydy import *
from sympy import symbols, sqrt

N = NewtonianReferenceFrame('N')


params = N.declare_parameters('a b c A B C D m g h')
q, qd = N.declare_coords('q', 5)
u, ud = N.declare_speeds('u', 3)

# Unpack the lists
q1, q2, q3, q4, q5 = q
q1d, q2d, q3d, q4d, q5d = qd
u1, u2, u3 = u
u1d, u2d, u3d = ud
a, b, c, A, B, C, D, m, g, h = params

R = N.rotate('R', 'BODY312', (q1, q2, q3), I=(A, B, C, D, 0, 0))

# Right hand sides of u's using time derivatives of coordinates
u_rhs = [dot(R.ang_vel(N), R[i]) for i in (1, 2, 3)]
# Matrix mapping qd to u
B = coefficient_matrix(u_rhs, qd[:3])

mu1 = dot(R[1], N[3])
mu2 = dot(R[2], N[3])
mu3 = dot(R[3], N[3])
eps = sqrt((a*mu1)**2 + (b*mu2)**2 + (c*mu3)**2)

x1 = a**2 * mu1 / eps
x2 = b**2 * mu2 / eps
x3 = c**2 * mu3 / eps

print 'x1 =', x1
print 'x2 =', x2
print 'x3 =', x3
print x1.diff(t)
stop

x1,x2,x3 = symbols('x1 x2 x3')
RO = N.O.locate('RO', -x1*R[1] - x2*R[2] + (h - x3)*R[3], frame=R, mass=m)

N1 = RO.locate('N1', x1*R[1] + x3*R[2] + (x3 - h)*R[3] - q4*N[1] - q5*N[2])

print RO.rel(N1)
stop

vron1 = dt(RO.rel(N1), N)
vron2 = cross(R.ang_vel(N), RO.rel(N.O))
eq1 = dot(vron1 - vron2, N[1])
eq2 = dot(vron1 - vron2, N[2])
print eq1
