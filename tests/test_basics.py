from PyDy import ReferenceFrame, dot, cross, UnitVector, identify
from sympy import symbols, S, Symbol, sin, cos, Matrix, eye, pprint
import numpy as N
A = ReferenceFrame('A')
e1 = UnitVector(A, 1)
e2 = UnitVector(A, 2)
e3 = UnitVector(A, 3)
#zero = UnitVector(A, 0)
zero = 0

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
    q1, q2, q3 = symbols('q1 q2 q3')
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

def test_cross_different_frames2():
    q1, q2, q3 = symbols('q1 q2 q3')
    N = ReferenceFrame('N')
    A = N.rotate('A', 3, q1)
    assert cross(N[1], A[1]) == sin(q1)*A[3]
    assert cross(N[1], A[2]) == cos(q1)*A[3]
    assert cross(N[1], A[1] + A[2]) == sin(q1)*A[3] + cos(q1)*A[3]

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

    assert C.get_frames_list(F) == [C, B, A, D, E, F]
    assert F.get_frames_list(C) == [F, E, D, A, B, C]

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
