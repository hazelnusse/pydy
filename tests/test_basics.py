from pydy import ReferenceFrame, dot, cross, UnitVector, express, \
        Vector, dt

from sympy import symbols, S, Symbol, Function, sin, cos, tan, Matrix, eye, \
    Rational, pprint, trigsimp, expand

zero = Vector({})

def test_Mul_order():
    A = ReferenceFrame('A')
    e1 = UnitVector(A, 1)
    e2 = UnitVector(A, 2)
    e3 = UnitVector(A, 3)
    assert e1*e2 == e1*e2
    assert e1*e2 != e2*e1

    assert e2*e1*e3 == e2*e1*e3
    assert e2*e1*e3 != e3*e2*e1

def test_UnitVector():
    A = ReferenceFrame('A')
    a1 = UnitVector(A, 1)
    a2 = UnitVector(A, 2)
    a3 = UnitVector(A, 3)

def test_dot_cross():
    A = ReferenceFrame('A')
    assert dot(A[1], A[1]) == 1
    assert dot(A[1], A[2]) == 0
    assert dot(A[1], A[3]) == 0

    assert dot(A[2], A[1]) == 0
    assert dot(A[2], A[2]) == 1
    assert dot(A[2], A[3]) == 0

    assert dot(A[3], A[1]) == 0
    assert dot(A[3], A[2]) == 0
    assert dot(A[3], A[3]) == 1

    assert cross(A[1], A[1]) == zero
    assert cross(A[1], A[2]) == A[3]
    assert cross(A[1], A[3]) == -A[2]

    assert cross(A[2], A[1]) == -A[3]
    assert cross(A[2], A[2]) == zero
    assert cross(A[2], A[3]) == A[1]

    assert cross(A[3], A[1]) == A[2]
    assert cross(A[3], A[2]) == -A[1]
    assert cross(A[3], A[3]) == zero

def test_expressions():
    A = ReferenceFrame('A')
    x, y = symbols("x y")
    e = x+x*A[1]+y+A[2]
    assert e == x+x*A[1]+y+A[2]
    assert e != x+x*A[1]+x+A[2]

def test_ReferenceFrame():
    A = ReferenceFrame('A')
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

def test_cross_method():
    q1, q2, q3, t = symbols('q1, q2, q3 t')
    N = ReferenceFrame('N')
    A = N.rotate('A', 1, q1)
    B = N.rotate('B', 2, q2)
    C = N.rotate('C', 3, q3)
    assert cross(N[1], N[1]) == Vector(0) == 0
    assert cross(N[1], N[2]) == N[3]
    assert N[1].cross(N[3]) == Vector({N[2]: -1})

    assert N[2].cross(N[1]) == Vector({N[3]: -1})
    assert N[2].cross(N[2]) == Vector(0)
    assert N[2].cross(N[3]) == N[1]

    assert N[3].cross(N[1]) == N[2]
    assert N[3].cross(N[2]) == Vector({N[1]: -1})
    assert N[3].cross(N[3]) == Vector(0)

    assert N[1].cross(A[1]) == Vector(0)
    assert N[1].cross(A[2]) == A[3]
    assert N[1].cross(A[3]) == Vector(-A[2])

    assert N[2].cross(A[1]) == Vector(-N[3])
    assert N[2].cross(A[2]) == Vector(sin(q1)*N[1])
    assert N[2].cross(A[3]) == Vector(cos(q1)*N[1])

    assert N[1].cross(B[1]) == Vector(sin(q2)*N[2])
    assert N[1].cross(B[2]) == N[3]
    assert N[1].cross(B[3]) == Vector(-cos(q2)*N[2])

def test_cross_different_frames3():
    q1, q2, q3 = symbols('q1 q2 q3')
    N = ReferenceFrame('N')
    A = N.rotate('A', 3, q1)
    B = A.rotate('B', 1, q2)
    C = B.rotate('C', 2, q3)
    assert cross(A[1], C[1]) == sin(q3)*C[2]
    #print "cross(A[1], C[2])" ,cross(A[1], C[2])
    #print "-sin(q3)*C[1]+cos(q3)*C[3]"
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
    #print A[1].express(C)
    #print cos(q3)*C[1] + sin(q3)*C[3]
    assert A[1].express(N) == Vector(cos(q1)*N[1] + sin(q1)*N[2])
    assert A[2].express(N) == Vector(-sin(q1)*N[1] + cos(q1)*N[2])
    #print "A[3].express(N) = ", A[3].express(N), "type(A[3].express(N)) = ",\
    #        type(A[3].express(N))
    assert A[3].express(N) == N[3]
    assert A[1].express(A) == A[1]
    assert A[2].express(A) == A[2]
    assert A[3].express(A) == A[3]
    assert A[1].express(B) == B[1]
    assert A[2].express(B) == Vector(cos(q2)*B[2] - sin(q2)*B[3])
    assert A[3].express(B) == Vector(sin(q2)*B[2] + cos(q2)*B[3])
    assert A[1].express(C) == Vector(cos(q3)*C[1] + sin(q3)*C[3])
    assert A[2].express(C) == Vector(sin(q2)*sin(q3)*C[1] + cos(q2)*C[2] - \
            sin(q2)*cos(q3)*C[3])
    assert A[3].express(C) == Vector(-sin(q3)*cos(q2)*C[1] + sin(q2)*C[2] + \
            cos(q2)*cos(q3)*C[3])

def test_express3():
    q3, q4, q5 = symbols('q3 q4 q5')
    N = ReferenceFrame('N')
    A = N.rotate('A', 3, q3)
    B = A.rotate('B', 1, q4)
    C = B.rotate('C', 2, q5)
    # Check to make sure UnitVectors get converted properly
    assert express(N[1], N) == N[1]
    assert express(N[2], N) == N[2]
    assert express(N[3], N) == N[3]
    assert express(N[1], A) == Vector(cos(q3)*A[1] - sin(q3)*A[2])
    assert express(N[2], A) == Vector(sin(q3)*A[1] + cos(q3)*A[2])
    assert express(N[3], A) == A[3]
    assert express(N[1], B) == Vector(cos(q3)*B[1] - sin(q3)*cos(q4)*B[2] + \
            sin(q3)*sin(q4)*B[3])
    assert express(N[2], B) == Vector(sin(q3)*B[1] + cos(q3)*cos(q4)*B[2] - \
            sin(q4)*cos(q3)*B[3])
    assert express(N[3], B) == Vector(sin(q4)*B[2] + cos(q4)*B[3])
    assert express(N[1], C) == Vector(
            (cos(q3)*cos(q5)-sin(q3)*sin(q4)*sin(q5))*C[1] -
            sin(q3)*cos(q4)*C[2] +
            (sin(q5)*cos(q3)+sin(q3)*sin(q4)*cos(q5))*C[3])
    assert express(N[2], C) == Vector(
            (sin(q3)*cos(q5) + sin(q4)*sin(q5)*cos(q3))*C[1] +
            cos(q3)*cos(q4)*C[2] +
            (sin(q3)*sin(q5) - sin(q4)*cos(q3)*cos(q5))*C[3])
    assert express(N[3], C) == Vector(-sin(q5)*cos(q4)*C[1] + sin(q4)*C[2] +
            cos(q4)*cos(q5)*C[3])

    assert express(A[1], N) == Vector(cos(q3)*N[1] + sin(q3)*N[2])
    assert express(A[2], N) == Vector(-sin(q3)*N[1] + cos(q3)*N[2])
    assert express(A[3], N) == N[3]
    assert express(A[1], A) == A[1]
    assert express(A[2], A) == A[2]
    assert express(A[3], A) == A[3]
    assert express(A[1], B) == B[1]
    assert express(A[2], B) == Vector(cos(q4)*B[2] - sin(q4)*B[3])
    assert express(A[3], B) == Vector(sin(q4)*B[2] + cos(q4)*B[3])
    assert express(A[1], C) == Vector(cos(q5)*C[1] + sin(q5)*C[3])
    assert express(A[2], C) == Vector(sin(q4)*sin(q5)*C[1] + cos(q4)*C[2] -
            sin(q4)*cos(q5)*C[3])
    assert express(A[3], C) == Vector(-sin(q5)*cos(q4)*C[1] + sin(q4)*C[2] +
            cos(q4)*cos(q5)*C[3])

    assert express(B[1], N) == Vector(cos(q3)*N[1] + sin(q3)*N[2])
    assert express(B[2], N) == Vector(-sin(q3)*cos(q4)*N[1] +
            cos(q3)*cos(q4)*N[2] + sin(q4)*N[3])
    assert express(B[3], N) == Vector(sin(q3)*sin(q4)*N[1] -
            sin(q4)*cos(q3)*N[2] + cos(q4)*N[3])
    assert express(B[1], A) == A[1]
    assert express(B[2], A) == Vector(cos(q4)*A[2] + sin(q4)*A[3])
    assert express(B[3], A) == Vector(-sin(q4)*A[2] + cos(q4)*A[3])
    assert express(B[1], B) == B[1]
    assert express(B[2], B) == B[2]
    assert express(B[3], B) == B[3]
    assert express(B[1], C) == Vector(cos(q5)*C[1] + sin(q5)*C[3])
    assert express(B[2], C) == C[2]
    assert express(B[3], C) == Vector(-sin(q5)*C[1] + cos(q5)*C[3])

    assert express(C[1], N) == Vector(
            (cos(q3)*cos(q5)-sin(q3)*sin(q4)*sin(q5))*N[1] +
            (sin(q3)*cos(q5)+sin(q4)*sin(q5)*cos(q3))*N[2] -
                sin(q5)*cos(q4)*N[3])
    assert express(C[2], N) == Vector(
            -sin(q3)*cos(q4)*N[1] + cos(q3)*cos(q4)*N[2] + sin(q4)*N[3])
    assert express(C[3], N) == Vector(
            (sin(q5)*cos(q3)+sin(q3)*sin(q4)*cos(q5))*N[1] +
            (sin(q3)*sin(q5)-sin(q4)*cos(q3)*cos(q5))*N[2] +
            cos(q4)*cos(q5)*N[3])
    assert express(C[1], A) == Vector(cos(q5)*A[1] + sin(q4)*sin(q5)*A[2] -
            sin(q5)*cos(q4)*A[3])
    assert express(C[2], A) == Vector(cos(q4)*A[2] + sin(q4)*A[3])
    assert express(C[3], A) == Vector(sin(q5)*A[1] - sin(q4)*cos(q5)*A[2] +
            cos(q4)*cos(q5)*A[3])
    assert express(C[1], B) == Vector(cos(q5)*B[1] - sin(q5)*B[3])
    assert express(C[2], B) == B[2]
    assert express(C[3], B) == Vector(sin(q5)*B[1] + cos(q5)*B[3])
    assert express(C[1], C) == C[1]
    assert express(C[2], C) == C[2]
    assert express(C[3], C) == C[3] == Vector(C[3])

    #  Check to make sure Vectors get converted back to UnitVectors
    assert N[1] == express(Vector(cos(q3)*A[1] - sin(q3)*A[2]), N)
    assert N[2] == express(Vector(sin(q3)*A[1] + cos(q3)*A[2]), N)
    # Trigsimp doesn't work properly here so the test fails
    #print express(Vector(cos(q3)*B[1] - sin(q3)*cos(q4)*B[2] + sin(q3)*sin(q4)*B[3]), N)
    #print trigsimp(cos(q4)**2 + sin(q3)**2*sin(q4)**2 + cos(q3)**2*cos(q4)**2*tan(q4)**2)
    #print trigsimp(cos(q4)**2 + cos(q3)**2*sin(q4)**2 + sin(q3)**2*sin(q4)**2)
    #assert N[1] == express(Vector(cos(q3)*B[1] - sin(q3)*cos(q4)*B[2] +
    #        sin(q3)*sin(q4)*B[3]), N)
    assert N[2] == express(Vector(sin(q3)*B[1] + cos(q3)*cos(q4)*B[2] -
        sin(q4)*cos(q3)*B[3]), N)
    assert N[3] == express(Vector(sin(q4)*B[2] + cos(q4)*B[3]), N)

    # Trigsimp also doesn't work on this one.
    #print express(Vector(
    #        (cos(q3)*cos(q5)-sin(q3)*sin(q4)*sin(q5))*C[1] -
    #        sin(q3)*cos(q4)*C[2] +
    #        (sin(q5)*cos(q3)+sin(q3)*sin(q4)*cos(q5))*C[3]), N)

    #assert N[1] == express(Vector(
    #        (cos(q3)*cos(q5)-sin(q3)*sin(q4)*sin(q5))*C[1] -
    #        sin(q3)*cos(q4)*C[2] +
    #        (sin(q5)*cos(q3)+sin(q3)*sin(q4)*cos(q5))*C[3]), N)
    # Trigsimp doesn't like this one either.
    #assert N[2] == express(Vector(
    #        (sin(q3)*cos(q5) + sin(q4)*sin(q5)*cos(q3))*C[1] +
    #        cos(q3)*cos(q4)*C[2] +
    #        (sin(q3)*sin(q5) - sin(q4)*cos(q3)*cos(q5))*C[3]), N)
    #print express(Vector(-sin(q5)*cos(q4)*C[1] + sin(q4)*C[2] + cos(q4)*cos(q5)*C[3]), N)
    #assert N[3] == express(Vector(-sin(q5)*cos(q4)*C[1] + sin(q4)*C[2] +
    #        cos(q4)*cos(q5)*C[3]), N)

    assert A[1] == express(Vector(cos(q3)*N[1] + sin(q3)*N[2]), A)
    assert A[2] == express(Vector(-sin(q3)*N[1] + cos(q3)*N[2]), A)

    assert A[2] == express(Vector(cos(q4)*B[2] - sin(q4)*B[3]), A)
    assert A[3] == express(Vector(sin(q4)*B[2] + cos(q4)*B[3]), A)

    assert A[1] == express(Vector(cos(q5)*C[1] + sin(q5)*C[3]), A)

    # Tripsimp messes up here too.
    #print express(Vector(sin(q4)*sin(q5)*C[1] + cos(q4)*C[2] -
    #        sin(q4)*cos(q5)*C[3]), A)
    #assert A[2] == express(Vector(sin(q4)*sin(q5)*C[1] + cos(q4)*C[2] -
    #        sin(q4)*cos(q5)*C[3]), A)

    assert A[3] == express(Vector(-sin(q5)*cos(q4)*C[1] + sin(q4)*C[2] +
            cos(q4)*cos(q5)*C[3]), A)
    assert B[1] == express(Vector(cos(q3)*N[1] + sin(q3)*N[2]), B)
    assert B[2] == express(Vector(-sin(q3)*cos(q4)*N[1] +
            cos(q3)*cos(q4)*N[2] + sin(q4)*N[3]), B)

    # Trigsimp messes up here too.
    #print express(Vector(sin(q3)*sin(q4)*N[1] - sin(q4)*cos(q3)*N[2]
    #    + cos(q4)*N[3]), B)
    #assert B[3] == express(Vector(sin(q3)*sin(q4)*N[1] -
    #        sin(q4)*cos(q3)*N[2] + cos(q4)*N[3]), B)

    assert B[2] == express(Vector(cos(q4)*A[2] + sin(q4)*A[3]), B)
    assert B[3] == express(Vector(-sin(q4)*A[2] + cos(q4)*A[3]), B)
    assert B[1] == express(Vector(cos(q5)*C[1] + sin(q5)*C[3]), B)
    #print express(Vector(-sin(q5)*C[1] + cos(q5)*C[3]), B)
    assert B[3] == express(Vector(-sin(q5)*C[1] + cos(q5)*C[3]), B)

    # Trigsimp fails here too.
    #print express(Vector(
    #        (cos(q3)*cos(q5)-sin(q3)*sin(q4)*sin(q5))*N[1] +
    #        (sin(q3)*cos(q5)+sin(q4)*sin(q5)*cos(q3))*N[2] -
    #            sin(q5)*cos(q4)*N[3]), C)
    #assert C[1] == express(Vector(
    #        (cos(q3)*cos(q5)-sin(q3)*sin(q4)*sin(q5))*N[1] +
    #        (sin(q3)*cos(q5)+sin(q4)*sin(q5)*cos(q3))*N[2] -
    #            sin(q5)*cos(q4)*N[3]), C)
    # Trigsimp fails here too.
    #print trigsimp(dot(express(Vector(
    #        -sin(q3)*cos(q4)*N[1] + cos(q3)*cos(q4)*N[2] + sin(q4)*N[3]), C),
    #        C[3]))
    #assert C[2] == express(Vector(
    #        -sin(q3)*cos(q4)*N[1] + cos(q3)*cos(q4)*N[2] + sin(q4)*N[3]), C)
    #print express(Vector(
    #        (sin(q5)*cos(q3)+sin(q3)*sin(q4)*cos(q5))*N[1] +
    #        (sin(q3)*sin(q5)-sin(q4)*cos(q3)*cos(q5))*N[2] +
    #        cos(q4)*cos(q5)*N[3]), C)
    #assert C[3] == express(Vector(
    #        (sin(q5)*cos(q3)+sin(q3)*sin(q4)*cos(q5))*N[1] +
    #        (sin(q3)*sin(q5)-sin(q4)*cos(q3)*cos(q5))*N[2] +
    #        cos(q4)*cos(q5)*N[3]), C)
    #print express(Vector(cos(q5)*A[1] + sin(q4)*sin(q5)*A[2] -
    #        sin(q5)*cos(q4)*A[3]), C)
    #assert C[1] == express(Vector(cos(q5)*A[1] + sin(q4)*sin(q5)*A[2] -
    #        sin(q5)*cos(q4)*A[3]), C)
    assert C[2] == express(Vector(cos(q4)*A[2] + sin(q4)*A[3]), C)
    assert C[3] == express(Vector(sin(q5)*A[1] - sin(q4)*cos(q5)*A[2] +
            cos(q4)*cos(q5)*A[3]), C)
    assert C[1] == express(Vector(cos(q5)*B[1] - sin(q5)*B[3]), C)
    assert C[3] == express(Vector(sin(q5)*B[1] + cos(q5)*B[3]), C)


def test_dt():
    t = Symbol('t')

    q1 = Function('q1')(t)
    q2 = Function('q2')(t)
    q3 = Function('q3')(t)
    q4 = Function('q4')(t)
    q5 = Function('q5')(t)

    N = ReferenceFrame('N')
    A = N.rotate('A', 3, q3)
    B = A.rotate('B', 1, q4)
    C = B.rotate('C', 2, q5)
    assert dt(N[1], N) == Vector({})
    assert dt(N[2], N) == Vector({})
    assert dt(N[3], N) == Vector({})
    assert dt(N[1], A) == Vector(-q3.diff(t)*N[2])
    assert dt(N[2], A) == Vector(q3.diff(t)*N[1])
    assert dt(N[3], A) == Vector({})
    assert dt(N[1], B) == Vector(-q3.diff(t)*N[2] + sin(q3)*q4.diff(t)*N[3])
    assert dt(N[2], B) == Vector(q3.diff(t)*N[1] - cos(q3)*q4.diff(t)*N[3])
    assert dt(N[3], B) == Vector(q4.diff(t)*A[2])
    assert express(dt(N[1], C), N) == Vector((-q3.diff(t) -
        sin(q4)*q5.diff(t))*N[2] + (sin(q3)*q4.diff(t) +
            cos(q3)*cos(q4)*q5.diff(t))*N[3])

    #print express(dt(N[2], C), N)
    #print Vector((q3.diff(t) +
    #    sin(q4)*q5.diff(t))*N[1] + (sin(q3)*cos(q4)*q5.diff(t) -
    #        cos(q3)*q4.diff(t))*N[3])
    assert express(dt(N[2], C), N) == Vector((q3.diff(t) +
        sin(q4)*q5.diff(t))*N[1] + (sin(q3)*cos(q4)*q5.diff(t) -
            cos(q3)*q4.diff(t))*N[3])

    assert dt(N[3], C) == Vector(q4.diff(t)*A[2] - cos(q4)*q5.diff(t)*B[1])

    #print dt(A[1],N)
    #print Vector(q3.diff(t)*A[2])
    assert dt(A[1], N) == Vector(q3.diff(t)*A[2]) == q3.diff(t)*A[2]
    assert dt(A[2], N) == Vector(-q3.diff(t)*A[1]) == -q3.diff(t)*A[1]
    assert dt(A[3], N) == Vector({}) == 0
    assert dt(A[1], A) == Vector({}) == 0
    assert dt(A[2], A) == Vector({}) == 0
    assert dt(A[3], A) == Vector({}) == 0
    assert dt(A[1], B) == Vector({}) == 0
    assert dt(A[2], B) == Vector(-q4.diff(t)*A[3]) == -q4.diff(t)*A[3]
    assert dt(A[3], B) == Vector(q4.diff(t)*A[2]) == q4.diff(t)*A[2]
    assert dt(A[1], C) == Vector(q5.diff(t)*B[3]) == q5.diff(t)*B[3]
    assert dt(A[2], C) == Vector(sin(q4)*q5.diff(t)*A[1] - q4.diff(t)*A[3]) ==\
            sin(q4)*q5.diff(t)*A[1] - q4.diff(t)*A[3]
    assert dt(A[3], C) == Vector(-cos(q4)*q5.diff(t)*A[1] + q4.diff(t)*A[2]) \
            == -cos(q4)*q5.diff(t)*A[1] + q4.diff(t)*A[2]

    assert dt(B[1], N) == Vector(cos(q4)*q3.diff(t)*B[2] -
            sin(q4)*q3.diff(t)*B[3]) == cos(q4)*q3.diff(t)*B[2] - \
                    sin(q4)*q3.diff(t)*B[3]
    assert dt(B[2], N) == Vector(-cos(q4)*q3.diff(t)*B[1] + q4.diff(t)*B[3]) \
            == -cos(q4)*q3.diff(t)*B[1] + q4.diff(t)*B[3]
    assert dt(B[3], N) == Vector(sin(q4)*q3.diff(t)*B[1] - q4.diff(t)*B[2]) ==\
            sin(q4)*q3.diff(t)*B[1] - q4.diff(t)*B[2]
    assert dt(B[1], A) == Vector({}) == 0
    assert dt(B[2], A) == Vector(q4.diff(t)*B[3]) == q4.diff(t)*B[3]
    assert dt(B[3], A) == Vector(-q4.diff(t)*B[2]) == -q4.diff(t)*B[2]


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
    G = E.rotate('G', 3, q1)
    H = G.rotate('H', 2, q3)
    I = N.rotate('I', 2, q2)

    assert N.get_frames_list(N) == [N]
    assert N.get_frames_list(A) == [N, A]
    assert N.get_frames_list(B) == [N, A, B]
    assert N.get_frames_list(C) == [N, A, B, C]
    assert N.get_frames_list(D) == [N, A, D]
    assert N.get_frames_list(E) == [N, A, D, E]
    assert N.get_frames_list(F) == [N, A, D, E, F]
    assert N.get_frames_list(G) == [N, A, D, E, G]
    assert N.get_frames_list(H) == [N, A, D, E, G, H]
    assert N.get_frames_list(I) == [N, I]

    assert A.get_frames_list(N) == [A, N]
    assert A.get_frames_list(A) == [A]
    assert A.get_frames_list(B) == [A, B]
    assert A.get_frames_list(C) == [A, B, C]
    assert A.get_frames_list(D) == [A, D]
    assert A.get_frames_list(E) == [A, D, E]
    assert A.get_frames_list(F) == [A, D, E, F]
    assert A.get_frames_list(G) == [A, D, E, G]
    assert A.get_frames_list(H) == [A, D, E, G, H]
    assert A.get_frames_list(I) == [A, N, I]

    assert B.get_frames_list(N) == [B, A, N]
    assert B.get_frames_list(A) == [B, A]
    assert B.get_frames_list(B) == [B]
    assert B.get_frames_list(C) == [B, C]
    assert B.get_frames_list(D) == [B, A, D]
    assert B.get_frames_list(E) == [B, A, D, E]
    assert B.get_frames_list(F) == [B, A, D, E, F]
    assert B.get_frames_list(G) == [B, A, D, E, G]
    assert B.get_frames_list(H) == [B, A, D, E, G, H]
    assert B.get_frames_list(I) == [B, A, N, I]

    assert C.get_frames_list(N) == [C, B, A, N]
    assert C.get_frames_list(A) == [C, B, A]
    assert C.get_frames_list(B) == [C, B]
    assert C.get_frames_list(C) == [C]
    assert C.get_frames_list(D) == [C, B, A, D]
    assert C.get_frames_list(E) == [C, B, A, D, E]
    assert C.get_frames_list(F) == [C, B, A, D, E, F]
    assert C.get_frames_list(G) == [C, B, A, D, E, G]
    assert C.get_frames_list(H) == [C, B, A, D, E, G, H]
    assert C.get_frames_list(I) == [C, B, A, N, I]

    assert D.get_frames_list(N) == [D, A, N]
    assert D.get_frames_list(A) == [D, A]
    assert D.get_frames_list(B) == [D, A, B]
    assert D.get_frames_list(C) == [D, A, B, C]
    assert D.get_frames_list(D) == [D]
    assert D.get_frames_list(E) == [D, E]
    assert D.get_frames_list(F) == [D, E, F]
    assert D.get_frames_list(G) == [D, E, G]
    assert D.get_frames_list(H) == [D, E, G, H]
    assert D.get_frames_list(I) == [D, A, N, I]

    assert E.get_frames_list(N) == [E, D, A, N]
    assert E.get_frames_list(A) == [E, D, A]
    assert E.get_frames_list(B) == [E, D, A, B]
    assert E.get_frames_list(C) == [E, D, A, B, C]
    assert E.get_frames_list(D) == [E, D]
    assert E.get_frames_list(E) == [E]
    assert E.get_frames_list(F) == [E, F]
    assert E.get_frames_list(G) == [E, G]
    assert E.get_frames_list(H) == [E, G, H]
    assert E.get_frames_list(I) == [E, D, A, N, I]

    assert F.get_frames_list(N) == [F, E, D, A, N]
    assert F.get_frames_list(A) == [F, E, D, A]
    assert F.get_frames_list(B) == [F, E, D, A, B]
    assert F.get_frames_list(C) == [F, E, D, A, B, C]
    assert F.get_frames_list(D) == [F, E, D]
    assert F.get_frames_list(E) == [F, E]
    assert F.get_frames_list(F) == [F]
    assert F.get_frames_list(G) == [F, E, G]
    assert F.get_frames_list(H) == [F, E, G, H]
    assert F.get_frames_list(I) == [F, E, D, A, N, I]

    assert G.get_frames_list(N) == [G, E, D, A, N]
    assert G.get_frames_list(A) == [G, E, D, A]
    assert G.get_frames_list(B) == [G, E, D, A, B]
    assert G.get_frames_list(C) == [G, E, D, A, B, C]
    assert G.get_frames_list(D) == [G, E, D]
    assert G.get_frames_list(E) == [G, E]
    assert G.get_frames_list(F) == [G, E, F]
    assert G.get_frames_list(G) == [G]
    assert G.get_frames_list(H) == [G, H]
    assert G.get_frames_list(I) == [G, E, D, A, N, I]

    assert H.get_frames_list(N) == [H, G, E, D, A, N]
    assert H.get_frames_list(A) == [H, G, E, D, A]
    assert H.get_frames_list(B) == [H, G, E, D, A, B]
    assert H.get_frames_list(C) == [H, G, E, D, A, B, C]
    assert H.get_frames_list(D) == [H, G, E, D]
    assert H.get_frames_list(E) == [H, G, E]
    assert H.get_frames_list(F) == [H, G, E, F]
    assert H.get_frames_list(G) == [H, G]
    assert H.get_frames_list(H) == [H]
    assert H.get_frames_list(I) == [H, G, E, D, A, N, I]

    assert I.get_frames_list(N) == [I, N]
    assert I.get_frames_list(A) == [I, N, A]
    assert I.get_frames_list(B) == [I, N, A, B]
    assert I.get_frames_list(C) == [I, N, A, B, C]
    assert I.get_frames_list(D) == [I, N, A, D]
    assert I.get_frames_list(E) == [I, N, A, D, E]
    assert I.get_frames_list(F) == [I, N, A, D, E, F]
    assert I.get_frames_list(G) == [I, N, A, D, E, G]
    assert I.get_frames_list(H) == [I, N, A, D, E, G, H]
    assert I.get_frames_list(I) == [I]

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
    B = A.rotate('B', 1, q2)
    for i in (1, 2, 3):
        for j in (1, 2, 3):
            a = cross(N[i], A[j])
            b = express(cross(N[i], A[j]), A)
            assert a == b

    for i in range(1, 4):
        for j in range(1, 4):
            a = cross(N[i], B[j])
            b = express(cross(B[j], N[i]), N)
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
    #print "v2: ", v2, "type(v2): ", type(v2), "v2.dict: ", v2.dict
    #print "v3: ", v3, "type(v3): ", type(v3), "v3.dict: ", v3.dict
    v4 = v2 + v3
    assert v4 == v2 + v3
    v3 = Vector({B[1]: q1*sin(q2), B[2]: t*u1*q1*sin(q2)})
    v5 = Vector({B[1]: q1*sin(q2), B[2]: t*u1*q1*sin(q2)})
    assert v3 == v5
    # Another way to generate the same vector
    v5 = Vector(q1*sin(q2)*B[1] + t*u1*q1*sin(q2)*B[2])
    assert v3 == v5
    assert v5.dict == {B[1]: q1*sin(q2), B[2]: t*u1*q1*sin(q2)}

def test_mag():
    A = ReferenceFrame('A')
    v1 = Vector(A[1])
    v2 = Vector(A[1] + A[2])
    v3 = Vector(A[1] + A[2] + A[3])
    v4 = -A[1]
    assert v1.mag() == 1
    #print v2.mag()
    assert v2.mag() == 2**Rational(1,2)
    assert v3.mag() == 3**Rational(1,2)
    assert v4.mag() == 1


def test_rotate_Euler_Space():
    t = Symbol('t')
    q1 = Function('q1')(t)
    q2 = Function('q2')(t)
    q3 = Function('q3')(t)
    A = ReferenceFrame('A')

    c1 = cos(q1)
    c2 = cos(q2)
    c3 = cos(q3)
    s1 = sin(q1)
    s2 = sin(q2)
    s3 = sin(q3)

    q1p = q1.diff(t)
    q2p = q2.diff(t)
    q3p = q3.diff(t)

    #### Rotation matrices from Spacecraft Dynamics, by Kane, Likins, Levinson,
    #### 1982, Appendix I, pg. 422


    ###########  DIRECTION COSINE MATRICES AS FUNCTIONS OF ORIENTATION ANGLES
    #### Euler Angles (Body Fixed rotations) ####
    #### Body 1-2-3
    B = A.rotate('B', 'BODY123', (q1, q2, q3))
    R123_Body = Matrix([ [ c2*c3, -c2*s3, s2],
                         [s1*s2*c3 + s3*c1, -s1*s2*s3 + c3*c1, -s1*c2],
                         [-c1*s2*c3 + s3*s1, c1*s2*s3 + c3*s1, c1*c2]])
    W_B_A = Vector((q1p*c2*c3 + q2p*s3)*B[1] + (-q1p*c2*s3+q2p*c3)*B[2] +
            (q1p*s2 + q3p)*B[3])
    #print B.get_rot_matrices(A)[0]
    #print R123_Body
    assert B.get_rot_matrices(A)[0] == R123_Body
    #print W_B_A
    #print 'b1>', dot(B.ang_vel(A), B[1])
    #print 'b2>', dot(B.ang_vel(A), B[2])
    #print 'b3>', dot(B.ang_vel(A), B[3])
    # Currently trigsimp can't handle this kind of thing...
    #assert B.ang_vel(A) == W_B_A

    #### Body 1-3-2
    B = A.rotate('B', 'BODY132', (q1, q2, q3))
    R132_Body = Matrix([ [c2*c3, -s2, c2*s3],
                         [c1*s2*c3 + s3*s1, c1*c2, c1*s2*s3 - c3*s1],
                         [s1*s2*c3 - s3*c1, s1*c2, s1*s2*s3 + c3*c1]])
    W_B_A = Vector((q1p*c2*c3 - q2p*s3)*B[1] + (-q1p*s2+q3p)*B[2] + (q1p*c2*s3
        + q2p*c3)*B[3])
    #print W_B_A
    #print B.ang_vel(A)
    #print B.get_rot_matrices(A)[0]
    #print R132_Body
    assert B.get_rot_matrices(A)[0] == R132_Body
    #assert B.ang_vel(A) == W_B_A

    #### Space Fixed rotations ####
    #### Space 1-2-3
    B = A.rotate('B', 'SPACE123', (q1, q2, q3))
    R123_Space = Matrix([ [c2*c3, s1*s2*c3 - s3*c1, c1*s2*c3 + s3*s1],
                          [c2*s3, s1*s2*s3 + c3*c1, c1*s2*s3 - c3*s1],
                          [-s2, s1*c2, c1*c2]])
    #print B.get_rot_matrices(A)[0]
    #print R123_Space
    assert B.get_rot_matrices(A)[0] == R123_Space

    #### Space 1-3-2
    B = A.rotate('B', 'SPACE132', (q1, q2, q3))
    R132_Space = Matrix([ [c2*c3, -c1*s2*c3 + s3*s1, s1*s2*c3 + s3*c1],
                          [s2, c1*c2, -s1*c2],
                          [-c2*s3, c1*s2*s3 + c3*s1, -s1*s2*s3 + c3*c1]])
    print B.get_rot_matrices(A)[0]
    print R132_Space
    assert B.get_rot_matrices(A)[0] == R132_Space


def test_Point_get_point_list():
    t = Symbol('t')
    x = Function('x')(t)
    l1, l2, l3 = symbols('l1 l2 l3')
    N = ReferenceFrame("N")
    A = N.rotate('A', 1, x)
    P1 = N.O.locate('P1', l1*N[1])
    P2 = N.O.locate('P2', l2*A[2])
    P3 = P2.locate('P3', l3*A[3])
    P4 = P3.locate('P4', 2*A[2] + 3*A[1])
    P5 = P3.locate('P5', 6*N[1] + 4*A[2])
    assert P1.get_point_list(N.O) == [P1, N.O]
    assert P1.get_point_list(P1) == [P1]
    assert P1.get_point_list(P2) == [P1, N.O, P2]
    assert P1.get_point_list(P3) == [P1, N.O, P2, P3]
    assert P1.get_point_list(P4) == [P1, N.O, P2, P3, P4]
    assert P1.get_point_list(P5) == [P1, N.O, P2, P3, P5]
    assert P2.get_point_list(N.O) == [P2, N.O]
    assert P2.get_point_list(P1) == [P2, N.O, P1]
    assert P2.get_point_list(P2) == [P2]
    assert P2.get_point_list(P3) == [P2, P3]
    assert P2.get_point_list(P4) == [P2, P3, P4]
    assert P2.get_point_list(P5) == [P2, P3, P5]
    assert P3.get_point_list(N.O) == [P3, P2, N.O]
    assert P3.get_point_list(P1) == [P3, P2, N.O, P1]
    assert P3.get_point_list(P2) == [P3, P2]
    assert P3.get_point_list(P3) == [P3]
    assert P3.get_point_list(P4) == [P3, P4]
    assert P3.get_point_list(P5) == [P3, P5]
    assert P4.get_point_list(N.O) == [P4, P3, P2, N.O]
    assert P4.get_point_list(P1) == [P4, P3, P2, N.O, P1]
    assert P4.get_point_list(P2) == [P4, P3, P2]
    assert P4.get_point_list(P3) == [P4, P3]
    assert P4.get_point_list(P4) == [P4]
    assert P4.get_point_list(P5) == [P4, P3, P5]
    assert P5.get_point_list(N.O) == [P5, P3, P2, N.O]
    assert P5.get_point_list(P1) == [P5, P3, P2, N.O, P1]
    assert P5.get_point_list(P2) == [P5, P3, P2]
    assert P5.get_point_list(P3) == [P5, P3]
    assert P5.get_point_list(P4) == [P5, P3, P4]
    assert P5.get_point_list(P5) == [P5]
