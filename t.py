"""
from sympy import Basic
from sympy.printing.str import StrPrinter
from sympy.printing.pretty.pretty import xsym

class GreenEggsAndHam(Basic):
    def __init__(self, string):
        self.s = string

    def __str__(self):
        return print_GreenEggsAndHam(self)

class HamPrinter(StrPrinter):
    def _print_GreenEggsAndHam(self, e):
        # The xsym('*') causes a UnicodeEncodeError
        return e.s.lower() + '*' + "\xC2\xB7"

def print_GreenEggsAndHam(e):
    pp = HamPrinter()
    return pp.doprint(e)

MyBreakfast = GreenEggsAndHam('I LOVE SYMPY')
print MyBreakfast
stop
"""

from sympy import *
from pydy import *
from sympy.printing.pretty.pretty import PrettyPrinter, xsym

q1 = Function('q1')(t)
q2 = Function('q2')(t)
q3 = Function('q3')(t)
A = ReferenceFrame('A')
B = A.rotate('B', 'BODY312', (q1, q2, q3))
print B.get_omega(A)
stop


class GeneralizedCoordinate(Symbol):
    def __init__(self, name, depends_on=Symbol('t'), *args):
        self.dv = depends_on
    """
    def diff(self, var=Symbol('t')):
        if var == self:
            print 'var == self'
            return 1
        elif var == self.dv:
            print 'var != self'
            return GeneralizedCoordinate(self.name+'\'')
        else:
            print 'none'

        if n == 1:
            print 'n = ', n
            return GeneralizedCoordinate(self.name+'\'')
        else:
            print 'n = ', n
            d = GeneralizedCoordinate(self.name+'\'')
            return d.dt(n=n-1)
    """
    def _eval_derivative(self, s):
        if s == self:
            return S.One
        elif s == self.dv:
            return GeneralizedCoordinate(self.name+'\'')
        else:
            return S.Zero

    def _latex_(self):
        return '\dot{'+str(self.name)+'}'


global _compact_trig
_compact_trig = True
t = Symbol('t')
q1 = GeneralizedCoordinate('q1')
q2 = GeneralizedCoordinate('q2')
u1 = Function('u1')(t)
#print 'q1 = ',q1
#print 'diff(q1, t) = ', q1.diff(t)
#print 'dt(sin(q1))', sin(q1).diff(t)
e = sin(q1).diff(t)
#print 'type(sin(q1).diff(t))', type(e)
#print latex(q1)
#print 'print_pydy(sin(q1)', print_pydy(sin(q1)), 'print sin(q1)', sin(q1)
#print print_pydy(cos(q1))
A = ReferenceFrame('A')
#print 'print_pydy(A[1])', print_pydy(A[1])
#print 'print A[1]', A[1]
#print 'print A[1]', Vector((1+sin(q1))*A[1] + cos(q1)*A[2])
#print 'Vector({})', Vector({})
v2 = Vector(q1*u1*A[1] + q2*t*sin(t)*A[2])
print v2
stop

print print_pydy(sin(q1)+cos(q1)+tan(q1))
print 'type(xsym("*"))', type(xsym("*")+'blah')
print 'xsym("*")', xsym("*"), type(print_pydy(sin(q1)))
stop




t = Symbol('t')
q1 = Function('q1')(t)
q2 = Function('q2')(t)
q3 = Function('q3')(t)

A = ReferenceFrame('A')
B = A.rotate('B', 'BODY123', (q1, q2, q3))


stop





print A.get_omega(N)
print N.get_omega(A)
stop



a = Function("a")(t)
b = Function("b")(t)
c = Function("c")(t)

A = ReferenceFrame('A')
B = A.rotate('B', 'BODY123', (a, b, c))
A_B = B.get_rot_matrices(A)[0]
A_B = A_B.subs({sin(a): Symbol('sa'), cos(a): Symbol('ca'),
                sin(b): Symbol('sb'), cos(b): Symbol('cb'),
                sin(c): Symbol('sc'), cos(c): Symbol('cc')})
print A_B
stop











t, x, y = symbols('t x y')
e1 = Matrix([1, 0, 0])
e2 = Matrix([0, 1, 0])
e3 = Matrix([0, 0, 1])
e1n = Matrix([-1, 0, 0])
e2n = Matrix([0, -1, 0])
e3n = Matrix([0, 0, -1])
zero = Matrix([0, 0, 0])

q1 = Function("q1")(t)
q2 = Function("q2")(t)
q3 = Function("q3")(t)
q4 = Function("q4")(t)
q5 = Function("q5")(t)
q6 = Function("q6")(t)

N = ReferenceFrame('N')
A = N.rotate('A', 3, q3)
B = A.rotate('B', 1, q4)
C = B.rotate('C', 2, q5)
#D = A.rotate('D', 2, q4)
#E = D.rotate('E', 1, q5)
#F = E.rotate('F', 3, q6)




r1, r2 = symbols('r1 r2')
p_no_co = Vector(q1*N[1] + q2*N[2] - r2*N[3] - r1*B[3])
p_co_cn = Vector(r2*N[3] + r1*B[3])
print 'p_no_co> =', p_no_co
print 'w_a_n> =', A.get_omega(N)
print 'w_b_a> =', B.get_omega(A)
print 'w_c_b> =', C.get_omega(B)
test = C.get_omega(N).cross(p_co_cn)
test_a = express(test, A)
test_b = express(test, B)
test_c = express(test, C)
test_da1 = dot(test, A[1])
test_da2 = dot(test, A[2])
test_da3 = dot(test, A[3])
print 'cross(W_C_N>, P_CO_CN>) = ', test
print 'dot(test,a1>): ', test_da1
print 'dot(test,a2>): ', test_da2
print 'dot(test,a3>): ', test_da3
print 'a1>: ', test_a.dict[A[1]]
print 'a2>: ', test_a.dict[A[2]]
print 'a3>: ', test_a.dict[A[3]]
print 'cross(W_C_N>, P_CO_CN>) (in A) = ', test_a
print 'cross(W_C_N>, P_CO_CN>) (in B) = ', test_b
print 'cross(W_C_N>, P_CO_CN>) (in C) = ', test_c

stop
CO = N._origin.locate('CO', Vector(q1*N[1] + q2*N[2] - r2*N[3] - r1*B[3]))
print 'position', CO.pos
print 'velocity', CO.vel
CN = CO.locate('CN', Vector(r1*B[3] + r2*N[3]), C)
print 'position', CN.pos
print 'velocity', CN.vel
stop



print 'A[1].dt = ', A[1].dt(N)
print 'A[2].dt = ', A[2].dt(N)
print 'A[3].dt = ', A[3].dt(N)

v1 = Vector(5*q1*A[1] + A[2])
print v1.dt(N)

v2 = Vector(5*q1*q2*A[1] + 3*q2*q3*sin(q1)*A[2])
print v2.dt(N)
print dot(v2.dt(N), A[1])
print dot(v2.dt(N), A[2])

test = [A[3], A[1], A[2]]
print test
test.sort()
print test
stop

'''
print "N[1].cross(N[1])", N[1].cross(N[1])
print "N[1].cross(N[2])", N[1].cross(N[2])
print "N[1].cross(N[3])", N[1].cross(N[3])
print N[1].cross(N[1]) == Vector(0)
print N[1].cross(N[2]) == N[3]
print N[1].cross(N[3]) == Vector({N[2]: -1})

print "N[2].cross(N[1])", N[2].cross(N[1])
print "N[2].cross(N[2])", N[2].cross(N[2])
print "N[2].cross(N[3])", N[2].cross(N[3])
print N[2].cross(N[1]) == Vector({N[3]: -1}) 
print N[2].cross(N[2]) == Vector(0)
print N[2].cross(N[3]) == N[1]

print "N[3].cross(N[1])", N[3].cross(N[1])
print "N[3].cross(N[2])", N[3].cross(N[2])
print "N[3].cross(N[3])", N[3].cross(N[3])
print N[3].cross(N[1]) == N[2]
print N[3].cross(N[2]) == Vector({N[1]: -1})
print N[3].cross(N[3]) == Vector(0)
'''

N = ReferenceFrame('N')
A = N.rotate('A', 3, q1)
B = A.rotate('B', 1, q2)

b = Vector({N[3]: -sin(q1)})
print '-b=',-b
print 'type(-b)=',type(-b)
negb = Vector(-b)

stop


print 'type(-B[1]) =', type(-B[1])
negb1 =Vector(-B[1])
print 'type(negb1) = ',type(negb1), 'negb1.dict =', negb1.dict,\
        'negb1._sympystr_() =',negb1._sympystr_()

negnegb1 = Vector(-negb1)
print 'type(negnegb1) =', type(negnegb1), 'negnegb1.dict =', negnegb1.dict, \
        'negnegb1._sympystr_() =', negnegb1._sympystr_()

test2 = negnegb1 + negb1
print 'type(test2) =', type(test2)
print 'test2.dict =',test2.dict
print 'test2 =',test2
print B[1] - B[1]
stop


print hash(A[1])
print hash(-A[1]+B[2])
print (-A[1])._mhash, type((-A[1])._mhash)
stop

for i in range(1, 4):
    for j in range(1, 4):
        a = cross(N[i], B[j])
        b = express(cross(B[j], N[i]), B)
        print "a=",a,"b=",b,type(b), type(-b)
        print a == -b

stop

A = N.rotate('A', 1, q1)
B = N.rotate('B', 2, q2)
C = N.rotate('C', 3, q3)

x = cross(A[1], A[3])
print "x =", x, "type(x) =", type(x), "x.dict =", x.dict, x.dict[A[2]]
y = Vector(-A[2])
print "y =", y, "type(y) =", type(y), "y.dict =", y.dict, y.dict[A[2]]
print "x == y: ", x==y

stop

print "N[1].cross(A[1])", N[1].cross(A[1])
print "N[1].cross(A[2])", N[1].cross(A[2])
print "N[1].cross(A[3])", N[1].cross(A[3])
print N[1].cross(A[1]) == Vector(0)
print N[1].cross(A[2]) == A[3]
print N[1].cross(A[3]) == Vector(-A[2])

print "N[1].cross(B[1])", N[1].cross(B[1])
print "N[1].cross(B[2])", N[1].cross(B[2])
print "N[1].cross(B[3])", N[1].cross(B[3])
print N[1].cross(B[1]) == Vector(sin(q2)*N[2])
print N[1].cross(B[2]) == N[3]
print N[1].cross(B[3]) == Vector(-cos(q2)*N[2])


#== Vector({A[2]: -1})

stop

print "N[2].cross(N[1])", N[2].cross(N[1])
print "N[2].cross(N[2])", N[2].cross(N[2])
print "N[2].cross(N[3])", N[2].cross(N[3])
print N[2].cross(N[1]) == Vector({N[3]: -1}) 
print N[2].cross(N[2]) == Vector(0)
print N[2].cross(N[3]) == N[1]

print "N[3].cross(N[1])", N[3].cross(N[1])
print "N[3].cross(N[2])", N[3].cross(N[2])
print "N[3].cross(N[3])", N[3].cross(N[3])
print N[3].cross(N[1]) == N[2]
print N[3].cross(N[2]) == Vector({N[1]: -1})
print N[3].cross(N[3]) == Vector(0)

stop

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
#v1 = Vector(sin(q1)*A[2])
#print "A[2].express(B)", A[2].express(B)
#print "v1 = ", v1
#print "v1.dict", v1.dict
#print "v1.express(A) = ", v1.express(A), "v1.express(B) = ", v1.express(B)
#v2 = v1.express(B)
#print "v2 =", v2
#print "v2.dict", v2.dict
#print type(A[3].express(N))
a1 = Vector(cos(q1)*N[1] + sin(q1)*N[2])
print "a1 =", a1, "a1.express(A) = ", a1.express(A)
print "type(a1.express(A))", type(a1.express(A))
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
