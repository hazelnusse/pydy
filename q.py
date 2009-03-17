from pydy import *
from sympy import *
q1, q2, q3 = symbols('q1 q2 q3')
N = ReferenceFrame('N')
A = N.rotate('A', 3, q1)
B = A.rotate('B', 1, q2)
C = B.rotate('C', 2, q3)

#print express(A[1], B)
#print express(A[2], B)
#print express(A[1]+A[2], B)
for i in range(1, 4):
    for j in range(1, 4):
        a = dot(N[i], B[j])
        b = dot(B[i], N[j])
        #b = express(cross(B[j], N[i]), B)
        #b = b.expand()
        print i, j, a, "|", b, a == b

#print A[1]*A[2]*B[1]
#print B[1]*A[2]*A[1]
