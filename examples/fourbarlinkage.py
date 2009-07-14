from sympy import *
from pydy import *

la, lb, lc, ln = symbols('la lb lc ln')
(q1, q2, q3), q_list, qdot_list = gcs('q', 3, list=True)

N = NewtonianReferenceFrame('N')
N.q_list = q_list
N.qdot_list = qdot_list
A = N.rotate('A', 3, q1)
B = N.rotate('B', 3, q2)
C = N.rotate('C', 3, q3)

AB = N.O.locate('AB', la*A[1], A)
BC = AB.locate('BC', lb*B[1], B)
CD = BC.locate('CD', lc*C[1])

kinematic_chain(CD, N.O, -ln*N[2])

for kc, dkc in (N.kc_eqs, N.dkc_eqs):
    print kc
    print dkc

#kindiffs = solve(N.dkc_eqs, q2.diff(t), q3.diff(t))
#print kindiffs
