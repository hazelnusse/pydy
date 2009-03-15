from PyDy import ReferenceFrame
from sympy import symbols, Function

N = ReferenceFrame('N')
A = ReferenceFrame('A')
B = ReferenceFrame('B')
C = ReferenceFrame('C')

g, r1, r2, t = symbols("g r1 r2 t")
q1 = Function("q1")(t)
q2 = Function("q2")(t)
q3 = Function("q3")(t)
q4 = Function("q4")(t)
q5 = Function("q5")(t)

print A[3]
P_NC_CO = r2*A[3] - r1*B[3]
print P_NC_CO
P_NO_CO = q1*N[1] + q2*N[2] + P_NC_CO

print P_NO_CO
