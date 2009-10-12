from sympy import Matrix, symbols, collect, integrate, pi
from pydy import *
Int = integrate
a, b, c, r1, r2, rho, m = symbols('a b c r1 r2 rho m')

# Volume of a torus
V = 2*pi**2*r1*r2**2

# A is the frame fixed in the torus, think of these as the global coordinates
A = NewtonianReferenceFrame('A')
# B is a frame that rotates about the axis out of the plane of symmetry of the
# torus
B = A.rotate('B', 3, a)
# C rotates about an axis along the tangent to the center of the circular cross
# section
C = B.rotate('C', 1, b)

# Locate an arbitrary point on the torus
p = Vector(r1*B[2] + c*C[2])

# Express it in the global coordinate system
x = dot(p, A[1])
y = dot(p, A[2])
z = dot(p, A[3])

# Determine the determinant of the Jacobian of the mapping from a, b, c to x,
# y, z
X = [a, b, c]
F = Matrix([x, y, z])
J = F.jacobian(X)
dv = (J.det(method='berkowitz')).expand().subs({cos(a)**2:1-sin(a)**2,
    cos(b)**2:1-sin(b)**2}).expand()

print 'dx*dy*dz =(', dv, ')*da*db*dc'

i11 = rho * dot(cross(p, A[1]), cross(p, A[1])) * dv
i22 = rho * dot(cross(p, A[2]), cross(p, A[2])) * dv
i33 = rho * dot(cross(p, A[3]), cross(p, A[3])) * dv

I11 = Int(Int(Int(i11, (a, 0, 2*pi)), (b, 0, 2*pi)), (c, 0, r2))
I22 = Int(Int(Int(i22, (a, 0, 2*pi)), (b, 0, 2*pi)), (c, 0, r2))
I33 = Int(Int(Int(i33, (a, 0, 2*pi)), (b, 0, 2*pi)), (c, 0, r2))

print 'I11 =', I11.subs({rho: m/V})
print 'I22 =', I22.subs({rho: m/V})
print 'I33 =', I33.subs({rho: m/V})
