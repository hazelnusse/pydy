from pydy import *
from sympy import *
import sympy
print sympy.__file__
t, x, y = symbols('t x y')
q1 = Function("q1")(t)
q2 = Function("q2")(t)
q3 = Function("q3")(t)

N = ReferenceFrame('N')
print N.ref_frame_list
#stop

#print N.W
#stop
A = N.rotate('A', 3, q1)
print A.ref_frame_list

#print "W_A_N> = ", A.W
B = A.rotate('B', 1, q2)
print B.ref_frame_list
#print "W_B_A> = ", B.W
C = B.rotate('C', 2, q3)
print C.ref_frame_list
stop
#print "W_C_B> = ", C.W
#stop


#print coeff(C.get_omega(N),[q1.diff(t),q2.diff(t),q3.diff(t)])
#print N.get_omega(C)

##
e = y + x*A[1] + x + A[2]
print e
print coeff(e,[x,y])
print [A[1]+1, 1]

stop
#print A.W
#stop
#print B.W
#print C.W

#print B.get_rot_matrices(A)
#print A.get_rot_matrices(B)
#print A.get_rot_matrices(C)

#e = B[1]*q1.diff(t)*B[1]
#e = B[1]*q1*B[1]
#print B[1].is_commutative
#print q2.is_commutative
#print q1.is_commutative
#print q1.diff(t).is_commutative
#print sin(q2).is_commutative
#print e

stop
e = cross(N[1], A[1])
print e
e = cross(N[1], A[2])
print e
e = cross(N[1], A[3])
print e
e = cross(N[2], A[1])
print e
e = cross(N[2], A[2])
print e
e = cross(N[2], A[3])
print e
e = cross(N[3], A[1])
print e
e = cross(N[3], A[2])
print e
e = cross(N[3], A[3])
print e


e = cross(A[1], A[1])
print e
