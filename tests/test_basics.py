from PyDy import ReferenceFrame, dot, cross, UnitVector
from sympy import symbols
import numpy as N
e1 = N.array([1.0,0.0,0.0])
e2 = N.array([0.0,1.0,0.0])
e3 = N.array([0.0,0.0,1.0])
e1n = N.array([-1.0,0.0,0.0])
e2n = N.array([0.0,-1.0,0.0])
e3n = N.array([0.0,0.0,-1.0])
zero = N.array([0.0,0.0,0.0])

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
    assert vec_eq(cross(A[1], A[1]), zero)
    assert vec_eq(cross(A[1], A[2]), e3)
    assert vec_eq(cross(A[1], A[3]), e2n)
    assert vec_eq(cross(A[2], A[1]), e3n)
    assert vec_eq(cross(A[2], A[2]), zero)
    assert vec_eq(cross(A[2], A[3]), e1)
    assert vec_eq(cross(A[3], A[1]), e2)
    assert vec_eq(cross(A[3], A[2]), e1n)
    assert vec_eq(cross(A[3], A[3]), zero)

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

def test_vector_init():
    x, y = symbols('x y')
    A=ReferenceFrame('A')
    e = y+x*A[1]+x+A[2]
    f = y+x*A[1]-A[2]
    a = UnitVector(e)
    b = UnitVector(f)
    dot(a, b)
