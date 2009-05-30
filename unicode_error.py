from sympy import *
from pydy import *
from sympy.printing.pretty.pretty import PrettyPrinter, xsym

t = Symbol('t')
q1 = Symbol('q1')
q2 = Symbol('q2')
u1 = Function('u1')(t)
A = ReferenceFrame('A')
v2 = Vector(q1*u1*A[1] + q2*t*sin(t)*A[2])
print v2
