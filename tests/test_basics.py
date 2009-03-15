from PyDy import ReferenceFrame, dot, cross, UnitVector, identify
from sympy import symbols, S, Symbol, sin, cos, Matrix
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

def test_cross_different_frames():
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

def test_get_rot_matrices():
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
    #B_C = Matrix([
    #    [cos(q3), 0, -sin(q3)],
    #    [0, 1, 0],
    #    [sin(q3), 0, cos(q3)]
    #    ])

    assert B.get_rot_matrices(A) == [B_A]
    assert A.get_rot_matrices(B) == [A_B]
    assert A.get_rot_matrices(C) == [A_B, B_C]
