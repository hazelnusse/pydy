from PyDy import ReferenceFrame, cross
from sympy import symbols, Function

g, r1, r2, t = symbols("g r1 r2 t")
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
print P_NC_CO
P_NO_CO = q1*N[1] + q2*N[2] + P_NC_CO
print P_NO_CO
W_A_N = u3*A[3]
print W_A_N
W_B_N = W_A_N + u4*A[1]
print W_B_N

W_C_N = W_B_N + u5*B[2]
print W_C_N

V_CO_N = cross(W_C_N, P_NC_CO)
print V_CO_N
