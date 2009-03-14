from PyDy import ReferenceFrame, dot, cross, UnitVector, identify
from sympy import symbols, S
import numpy as N
e1 = UnitVector('A', 1)
e2 = UnitVector('A', 2)
e3 = UnitVector('A', 3)
zero = UnitVector("A", 0)

def eq(a, b, eps=1e-8):
    return abs(a-b) < eps

def vec_eq(a, b, eps=1e-8):
    return N.dot((a.v['num'] - b),(a.v['num'] - b)) < eps

def test_UnitVector():
    a = UnitVector("A", 3)

def test_dot_cross():
    A = ReferenceFrame('A')
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
	A = ReferenceFrame('A')
	e = x+x*A[1]+y+A[2]
	assert e == x+x*A[1]+y+A[2]
	assert e != x+x*A[1]+x+A[2]

def test_coeff():
    x, y = symbols('x y')
    A=ReferenceFrame('A')
    e = y+x*A[1]+x+A[2]
    assert e.coeff(x) == 1+A[1]
    assert e.coeff(y) == 1
    assert e.coeff(A[1]) == x
    assert e.coeff(A[2]) == 1

def test_identify():
    A = ReferenceFrame('A')
    x = symbols("x")
    assert set(identify(A[1]*A[2])) == set((S(1), A[1], A[2]))
    assert set(identify(x*A[1]*A[2])) == set((x, A[1], A[2]))
    assert set(identify(x*A[1]*A[1])) == set((x, A[1], A[1]))
