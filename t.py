from pydy import *
from sympy import *
import sympy
print sympy.__file__
t, x, y = symbols('t x y')
q1 = Function("q1")(t)
q2 = Function("q2")(t)
q3 = Function("q3")(t)
q4 = Function("q4")(t)
q5 = Function("q5")(t)
q6 = Function("q6")(t)

N = ReferenceFrame('N')
A = N.rotate('A', 3, q1)
B = A.rotate('B', 1, q2)
C = B.rotate('C', 2, q3)
D = A.rotate('D', 2, q4)
E = D.rotate('E', 1, q5)
F = E.rotate('F', 3, q6)

'''
print "N.ref_frame_list", N.ref_frame_list
print "A.ref_frame_list", A.ref_frame_list
print "N to N", N.get_frames_list(N)
print "N to A", N.get_frames_list(A)
print "A to A", A.get_frames_list(A)
print "A to N", A.get_frames_list(N)
print "B to B", B.get_frames_list(B)
print "B to A", B.get_frames_list(A)
print "B to N", B.get_frames_list(N)
print "A to B", A.get_frames_list(B)
print "N to B", N.get_frames_list(B)
print "C to D: ", C.get_frames_list(D)
print "C to E: ", C.get_frames_list(E)
print "C to F: ", C.get_frames_list(F)
print "D to C: ", D.get_frames_list(C)
print "E to C: ", E.get_frames_list(C)
print "F to C: ", F.get_frames_list(C)
'''

'''
Testing the parse_terms(function)
print "A[1]: ", type(A[1])
print parse_terms(A[1])
print "S(0): ", type(S(0))
print parse_terms(S(0))
print "sin(q1)*q2*A[3]: ", type(sin(q1)*q2*A[3])
print parse_terms(sin(q1)*q2*A[3])
print "sin(q1)*sin(q1)*q2*A[3]: ", type(sin(q1)*sin(q1)*q2*A[3])
print parse_terms(sin(q1)*sin(q1)*q2*A[3])
test = sin(q1)*sin(q1)*q2*A[3] + q1*A[2]
print test, type(test)
print parse_terms(test)
test = sin(q1)*sin(q1)*q2*A[3] + q1*A[2] + S(0) + cos(q3)*A[2]
print test, type(test)
print parse_terms(test)
stop
'''

'''
x, y, z = symbols('x y z')
q1 = Function('q1')(t)

v = Vector(x*A[1])
print "v = ", v, "type(v) = ", type(v)
print "v.dict = ", v.dict

v = Vector(x*A[2]+y*A[3]+x*z*B[2])
print "v = ", v
print "v.dict = ", v.dict

v2 = Vector(y*A[1] + x*B[3] + z*A[2])
print "v2 = ", v2

v3 = v + v2
print "v3 = v + v2 = ", v3
print type(v3)
print v3.dict

v1 = Vector({A[1]: x, B[2]: y})
v2 = Vector({B[2]: y, A[1]: x})
print "v1 == v2 : ", v2 == Vector({A[1]: x, B[2]: y}) 
'''

'''
print "N[1].express(N) = ", N[1].express(N)
print "N[1].express(A) = ", N[1].express(A)
print "N[2].express(N) = ", N[2].express(N)
print "N[2].express(A) = ", N[2].express(A)
print "N[3].express(N) = ", N[3].express(N)
print "N[3].express(A) = ", N[3].express(A)
print "N[1].dot(A[1]) = ", N[1].dot(A[1])
print "N[1].dot(A[2]) = ", N[1].dot(A[2])
print "N[1].dot(A[3]) = ", N[1].dot(A[3])
print "N[2].dot(A[1]) = ", N[2].dot(A[1])
print "N[2].dot(A[2]) = ", N[2].dot(A[2])
print "N[2].dot(A[3]) = ", N[2].dot(A[3])
print "N[3].dot(A[1]) = ", N[3].dot(A[1])
print "N[3].dot(A[2]) = ", N[3].dot(A[2])
print "N[3].dot(A[3]) = ", N[3].dot(A[3])


print "A[1].express(N) = ", A[1].express(N)
print "A[1].express(A) = ", A[1].express(N)
print "A[1].dot(N[1]) = ", A[1].dot(N[1])
print "A[1].dot(N[2]) = ", A[1].dot(N[2])
print "A[1].dot(N[3]) = ", A[1].dot(N[3])
print "A[2].dot(N[1]) = ", A[2].dot(N[1])
print "A[2].dot(N[2]) = ", A[2].dot(N[2])
print "A[2].dot(N[3]) = ", A[2].dot(N[3])
print "A[3].dot(N[1]) = ", A[3].dot(N[1])
print "A[3].dot(N[2]) = ", A[3].dot(N[2])
print "A[3].dot(N[3]) = ", A[3].dot(N[3])
'''

n_one = UnitVector(N,1)

#print "2*N[1] =", 2*n_one
#print 'type(N[1]) =', type(N[1])
#print 'q1*N[1] =', q1*N[1]
#print 'type(q1*N[1]) =', type(q1*N[1])
v = Vector(q1*N[1] + q2*N[2])
#print "type(v) ", type(v)
#print "v =", v

v2 = q1*v
#print 'v2 = q1*v =', v2
#print 'type(v2) ', type(v2)
#print 'v2.args = ', v2.args, 'type(v2.args) =', type(v2.args)
#print v2.args[0]

v3 = Vector(q1*v)
#print "type(v3) =", type(v3)
#print "v3 =", v3

#v4 = Vector(q1*v2) + N[3] + q1*v2
#print "v4 = Vector(q1*v2) + N[3] + q1*v2", v4
#print "type(v4)", type(v4)

#stop


v = Vector(q1*sin(q2)*A[1]+ q2*cos(q3)*F[3])

#print "v", v
#print "A[1].dot(N[1]) = ", A[1].dot(N[1])
z = Vector({})
#print "z", z, "z.dict", z.dict


# Expressing Vector instances in different frames
A = N.rotate('A',3,q1)
B = A.rotate('B',1,q2)
v1 = Vector(sin(q1)*A[2])
print "A[2].express(B)", A[2].express(B)
print "v1 = ", v1
print "v1.dict", v1.dict
print "v1.express(A) = ", v1.express(A), "v1.express(B) = ", v1.express(B)
v2 = v1.express(B)
print "v2 =", v2
print "v2.dict", v2.dict
stop



#print "v = ", v
#print "A[1].dot(v) = ", A[1].dot(v)

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
