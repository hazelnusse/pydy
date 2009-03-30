from pydy import ReferenceFrame, dot, cross, UnitVector, identify, express, \
        coeff, Vector
from sympy import symbols, S, Symbol, Function, sin, cos, Matrix, eye, pprint
import numpy as N
A = ReferenceFrame('A')
e1 = UnitVector(A, 1)
e2 = UnitVector(A, 2)
e3 = UnitVector(A, 3)
#zero = UnitVector(A, 0)
zero = 0

def test_Mul_order():
    assert e1*e2 == e1*e2
    assert e1*e2 != e2*e1

    assert e2*e1*e3 == e2*e1*e3
    assert e2*e1*e3 != e3*e2*e1

def eq(a, b, eps=1e-8):
    return abs(a-b) < eps

def vec_eq(a, b, eps=1e-8):
    return N.dot((a.v['num'] - b),(a.v['num'] - b)) < eps

def test_UnitVector():
    a = UnitVector(A, 3)

def test_dot_cross():
    assert eq(dot(A[1], A[1]), 1.0)
    assert eq(dot(A[1], A[2]), 0.0)
    assert eq(dot(A[1], A[3]), 0.0)
    assert eq(dot(A[2], A[1]), 0.0)
    assert eq(dot(A[2], A[2]), 1.0)
    assert eq(dot(A[2], A[3]), 0.0)
    assert eq(dot(A[3], A[1]), 0.0)
    assert eq(dot(A[3], A[2]), 0.0)
    assert eq(dot(A[3], A[3]), 1.0)
    assert cross(A[1], A[1]) == zero
    assert cross(A[1], A[2]) == e3
    assert cross(A[1], A[3]) == -e2
    assert cross(A[2], A[1]) == -e3
    assert cross(A[2], A[2]) == zero
    assert cross(A[2], A[3]) == e1
    assert cross(A[3], A[1]) == e2
    assert cross(A[3], A[2]) == -e1
    assert cross(A[3], A[3]) == zero

def test_expressions():
	x, y = symbols("x y")
	e = x+x*A[1]+y+A[2]
	assert e == x+x*A[1]+y+A[2]
	assert e != x+x*A[1]+x+A[2]

def test_coeff():
    x, y = symbols('x y')
    e = y+x*A[1]+x+A[2]
    assert e.coeff(x) == 1+A[1]
    assert e.coeff(y) == 1
    assert e.coeff(A[1]) == x
    assert e.coeff(A[2]) == 1
    assert coeff(e,[x,y]) == [A[1]+1, 1]

def test_identify():
    x = symbols("x")
    assert set(identify(A[1]*A[2])) == set((S(1), A[1], A[2]))
    assert set(identify(x*A[1]*A[2])) == set((x, A[1], A[2]))
    assert set(identify(x*A[1]*A[1])) == set((x, A[1], A[1]))

def test_ReferenceFrame():
    phi = Symbol("phi")
    B = A.rotate("B", 1, phi)
    assert B.transforms[A] is not None
    B = A.rotate("B", 2, phi)
    assert B.transforms[A] is not None
    B = A.rotate("B", 3, phi)
    assert B.transforms[A] is not None

def test_cross_different_frames1():
    q1, q2, q3, t = symbols('q1 q2 q3 t')
    u1 = Function("u1")(t)
    N = ReferenceFrame('N')
    A = N.rotate('A', 3, q1)
    assert cross(N[1], A[1]) == sin(q1)*A[3]
    assert cross(N[1], A[2]) == cos(q1)*A[3]
    assert cross(N[1], A[3]) == -sin(q1)*A[1]-cos(q1)*A[2]
    assert cross(N[2], A[1]) == -cos(q1)*A[3]
    assert cross(N[2], A[2]) == sin(q1)*A[3]
    assert cross(N[2], A[3]) == cos(q1)*A[1]-sin(q1)*A[2]
    assert cross(N[3], A[1]) == A[2]
    assert cross(N[3], A[2]) == -A[1]
    assert cross(N[3], A[3]) == 0
    assert cross(u1*N[1], A[1]) == sin(q1)*u1*A[3]

def test_cross_different_frames3():
    q1, q2, q3 = symbols('q1 q2 q3')
    N = ReferenceFrame('N')
    A = N.rotate('A', 3, q1)
    B = A.rotate('B', 1, q2)
    C = B.rotate('C', 2, q3)
    assert cross(A[1], C[1]) == sin(q3)*C[2]
    assert cross(A[1], C[2]) == -sin(q3)*C[1]+cos(q3)*C[3]
    assert cross(A[1], C[3]) == -cos(q3)*C[2]

def test_express1():
    q1, q2, q3 = symbols('q1 q2 q3')
    N = ReferenceFrame('N')
    A = N.rotate('A', 3, q1)
    B = A.rotate('B', 1, q2)
    C = B.rotate('C', 2, q3)
    assert express(A[1], C) == cos(q3)*C[1] + sin(q3)*C[3]
    assert express(A[2], C) == sin(q2)*sin(q3)*C[1] + cos(q2)*C[2] - \
            sin(q2)*cos(q3)*C[3]
    assert express(A[3], C) == -sin(q3)*cos(q2)*C[1] + sin(q2)*C[2] + \
            cos(q2)*cos(q3)*C[3]

def test_express2():
    q1, q2, q3 = symbols('q1 q2 q3')
    N = ReferenceFrame('N')
    A = N.rotate('A', 3, q1)
    B = A.rotate('B', 1, q2)
    C = B.rotate('C', 2, q3)
    print A[1].express(C)
    print cos(q3)*C[1] + sin(q3)*C[3]
    assert A[1].express(C) == Vector(cos(q3)*C[1] + sin(q3)*C[3])
    assert A[2].express(C) == Vector(sin(q2)*sin(q3)*C[1] + cos(q2)*C[2] - \
            sin(q2)*cos(q3)*C[3])
    assert A[3].express(C) == Vector(-sin(q3)*cos(q2)*C[1] + sin(q2)*C[2] + \
            cos(q2)*cos(q3)*C[3])

def test_cross_different_frames2():
    q1, q2, q3 = symbols('q1 q2 q3')
    N = ReferenceFrame('N')
    A = N.rotate('A', 3, q1)
    assert cross(N[1], A[1]) == sin(q1)*A[3]
    assert cross(N[1], A[2]) == cos(q1)*A[3]
    assert cross(N[1], A[1] + A[2]) == sin(q1)*A[3] + cos(q1)*A[3]
    #assert cross(A[1] + A[2], N[1]) == -sin(q1)*A[3] - cos(q1)*A[3]

def test_get_frames_list1():
    q1, q2, q3 = symbols('q1 q2 q3')
    N = ReferenceFrame('N')
    A = N.rotate('A', 3, q1)
    B = A.rotate('B', 1, q2)
    C = B.rotate('C', 2, q3)

    assert B.get_frames_list(A) == [B, A]
    assert A.get_frames_list(B) == [A, B]
    assert A.get_frames_list(C) == [A, B, C]
    assert A.get_frames_list(C) == [A, B, C]
    assert C.get_frames_list(A) == [C, B, A]

def test_get_frames_list2():
    q1, q2, q3 = symbols('q1 q2 q3')
    N = ReferenceFrame('N')
    A = N.rotate('A', 3, q1)
    B = A.rotate('B', 1, q2)
    C = A.rotate('C', 2, q3)

    assert B.get_frames_list(C) == [B, A, C]
    assert C.get_frames_list(B) == [C, A, B]

def test_get_frames_list3():
    q1, q2, q3 = symbols('q1 q2 q3')
    N = ReferenceFrame('N')
    A = N.rotate('A', 3, q1)
    B = A.rotate('B', 1, q2)
    C = B.rotate('C', 2, q3)
    D = A.rotate('D', 2, q3)
    E = D.rotate('E', 2, q3)
    F = E.rotate('F', 2, q3)
    assert N.get_frames_list(N) == [N]
    assert N.get_frames_list(A) == [N, A]
    assert A.get_frames_list(A) == [A]
    assert A.get_frames_list(N) == [A, N]
    assert B.get_frames_list(B) == [B]
    assert B.get_frames_list(A) == [B, A]
    assert B.get_frames_list(N) == [B, A, N]
    assert A.get_frames_list(B) == [A, B]
    assert N.get_frames_list(B) == [N, A, B]
    assert C.get_frames_list(D) == [C, B, A, D]
    assert C.get_frames_list(E) == [C, B, A, D, E]
    assert C.get_frames_list(F) == [C, B, A, D, E, F]
    assert D.get_frames_list(C) == [D, A, B, C]
    assert E.get_frames_list(C) == [E, D, A, B, C]
    assert F.get_frames_list(C) == [F, E, D, A, B, C]

def test_get_frames_list4():
    q1, q2, q3 = symbols('q1 q2 q3')
    N = ReferenceFrame('N')
    A = N.rotate('A', 3, q1)
    B = A.rotate('B', 1, q2)
    C = B.rotate('C', 2, q3)
    D = A.rotate('D', 2, q3)
    E = D.rotate('E', 2, q3)
    F = E.rotate('F', 2, q3)

    assert B.get_frames_list(N) == [B, A, N]
    assert N.get_frames_list(B) == [N, A, B]
    assert C.get_frames_list(N) == [C, B, A, N]
    assert N.get_frames_list(C) == [N, A, B, C]

def test_get_rot_matrices1():
    q1, q2, q3 = symbols('q1 q2 q3')
    N = ReferenceFrame('N')
    A = N.rotate('A', 3, q1)
    B = A.rotate('B', 1, q2)
    C = B.rotate('C', 2, q3)

    B_A = Matrix([
        [1, 0, 0],
        [0, cos(q2), sin(q2)],
        [0, -sin(q2), cos(q2)]
        ])
    A_B = Matrix([
        [1, 0, 0],
        [0, cos(q2), -sin(q2)],
        [0, sin(q2), cos(q2)]
        ])
    B_C = Matrix([
        [cos(q3), 0, sin(q3)],
        [0, 1, 0],
        [-sin(q3), 0, cos(q3)]
        ])
    C_B = Matrix([
        [cos(q3), 0, -sin(q3)],
        [0, 1, 0],
        [sin(q3), 0, cos(q3)]
        ])

    assert B.get_rot_matrices(B) == [eye(3)]
    assert B.get_rot_matrices(A) == [A_B]
    assert A.get_rot_matrices(B) == [B_A]
    assert A.get_rot_matrices(C) == [C_B, B_A]
    assert C.get_rot_matrices(A) == [A_B, B_C]

def test_get_rot_matrices2():
    q1, q2, q3, q4, q5, q6 = symbols('q1 q2 q3 q4 q5 q6')
    N = ReferenceFrame('N')
    A = N.rotate('A', 3, q1)
    B = A.rotate('B', 1, q2)
    C = B.rotate('C', 2, q3)
    D = A.rotate('D', 2, q4)
    E = D.rotate('E', 1, q5)
    F = E.rotate('F', 3, q6)

    A_B = Matrix([
        [1, 0, 0],
        [0, cos(q2), -sin(q2)],
        [0, sin(q2), cos(q2)],
        ])
    B_A = A_B.T

    B_C = Matrix([
        [cos(q3), 0, sin(q3)],
        [0, 1, 0],
        [-sin(q3), 0, cos(q3)],
        ])
    C_B = B_C.T

    A_D = Matrix([
        [cos(q4), 0, sin(q4)],
        [0, 1, 0],
        [-sin(q4), 0, cos(q4)],
        ])
    D_A = A_D.T

    D_E = Matrix([
        [1, 0, 0],
        [0, cos(q5), -sin(q5)],
        [0, sin(q5), cos(q5)],
        ])
    E_D = D_E.T

    E_F = Matrix([
        [cos(q6), -sin(q6), 0],
        [sin(q6), cos(q6), 0],
        [0, 0, 1],
        ])
    F_E = E_F.T

    assert C.get_rot_matrices(F) == [F_E, E_D, D_A, A_B, B_C]
    assert F.get_rot_matrices(C) == [C_B, B_A, A_D, D_E, E_F]

def test_cross2():
    q1, q2, q3 = symbols('q1 q2 q3')
    N = ReferenceFrame('N')
    A = N.rotate('A', 3, q1)
    B = A.rotate('A', 1, q2)
    for i in range(1, 4):
        for j in range(1, 4):
            a = cross(N[i], A[j])
            b = express(cross(A[j], N[i]), A)
            assert a == -b

    for i in range(1, 4):
        for j in range(1, 4):
            a = cross(N[i], B[j])
            b = express(cross(B[j], N[i]), B)
            assert a == -b

def test_dot2():
    q1, q2, q3 = symbols('q1 q2 q3')
    print "q1", q1
    N = ReferenceFrame('N')
    A = N.rotate('A', 3, q1)
    B = A.rotate('A', 1, q2)
    for i in range(1, 4):
        for j in range(1, 4):
            a = dot(N[i], A[j])
            b = dot(A[j], N[i])
            assert a == b

    for i in range(1, 4):
        for j in range(1, 4):
            a = dot(N[i], B[j])
            b = dot(B[j], N[i])
            assert a == b

def test_Vector_class():
    q1, q2, q3, t = symbols('q1 q2 q3 t')
    u1 = Function('u1')(t)
    N = ReferenceFrame('N')
    A = N.rotate('A', 3, q1)
    B = A.rotate('B', 1, q2)
    v1 = Vector(0)
    v2 = Vector(q1*u1*A[1] + q2*t*sin(t)*A[2])
    v3 = Vector({B[1]: q1*sin(q2), B[2]: t*u1*q1*sin(q2)})
    # Basic functionality tests
    assert v1.parse_terms(A[1]) == {A[1]: 1}
    assert v1.parse_terms(0) == {}
    assert v1.parse_terms(S(0)) == {}
    assert v1.parse_terms(q1*sin(t)*A[1] + A[2]*cos(q2)*u1) == {A[1]: \
            q1*sin(t), A[2]: cos(q2)*u1}
    test = sin(q1)*sin(q1)*q2*A[3]  + q1*A[2] + S(0) + cos(q3)*A[2]
    assert v1.parse_terms(test) == {A[3]: sin(q1)*sin(q1)*q2, \
            A[2]:cos(q3) + q1}
    # Equality tests
    v4 = v2 + v3
    assert v4 == v2 + v3
    v3 = Vector({B[1]: q1*sin(q2), B[2]: t*u1*q1*sin(q2)})
    v5 = Vector({B[1]: q1*sin(q2), B[2]: t*u1*q1*sin(q2)})
    assert v3 == v5
    # Another way to generate the same vector
    v5 = Vector(q1*sin(q2)*B[1] + t*u1*q1*sin(q2)*B[2])
    assert v3 == v5
    assert v5.dict == {B[1]: q1*sin(q2), B[2]: t*u1*q1*sin(q2)}
