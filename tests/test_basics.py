from PyDy import ReferenceFrame, dot

def eq(a, b, eps=1e-8):
    return abs(a-b) < eps

def test_UnitVector():
	A = ReferenceFrame('A')
	test1 = dot(A[1],A[2])
	assert eq(test1, 0.0)
