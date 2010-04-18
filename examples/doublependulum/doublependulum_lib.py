from __future__ import division
from math import sin, cos

def eoms(_x, t, _params):
    """Double Pendulum equations of motion.

    _x is an array/list in the following order:
        q1:  Angle of first pendulum link relative to vertical (0 downwards)
        q2:  Angle of second pendulum link relative to first (0 is fully extended)
        u1:  A[1] measure number of the inertial angular velocity of the first link.
        u2:  A[1] measure number of the inertial angular velocity of the second link.

    _params is an array/list in the following order:
        m1:  Mass of first pendulum point mass.
        m2:  Mass of second pendulum point mass.
        l1:  Length of first pendulum link.
        l2:  Length of second pendulum link.
         g:  Gravitational constant.

    """
    # Unpack function arguments
    q1, q2, u1, u2 = _x

    # Unpack function parameters
    m1, m2, g, l1, l2, b1, b2 = _params

    # Trigonometric functions
    c1 = cos(q1)
    c2 = cos(q2)
    s1 = sin(q1)
    s2 = sin(q2)

    # Nested terms
    _M00 = -m1*l1**2 - m2*l1**2
    _M01 = -l1*l2*m2*c2
    _M11 = -m2*l2**2
    _rhs0 = b1*u1 + b2*u1 - b2*u2 + g*l1*m1*s1 + g*l1*m2*s1 - l1*l2*m2*u2**2*s2
    _rhs1 = b2*u2 - b2*u1 + g*l2*m2*c1*s2 + g*l2*m2*c2*s1 + l1*l2*m2*u1**2*s2

    # Calculate return values
    q1d = u1
    q2d = -u1 + u2
    u1d = (_M11*_rhs0 - _M01*_rhs1)/(_M00*_M11 - _M01**2)
    u2d = (_M00*_rhs1 - _M01*_rhs0)/(_M00*_M11 - _M01**2)

    # Return calculated values
    return [q1d, q2d, u1d, u2d]

def energy(_x, _params):
    """Kinetic and Potential Energy of double pendulum.

    _x is an array/list in the following order:
        q1:  Angle of first pendulum link relative to vertical (0 downwards)
        q2:  Angle of second pendulum link relative to first (0 is fully extended)
        u1:  A[1] measure number of the inertial angular velocity of the first link.
        u2:  A[1] measure number of the inertial angular velocity of the second link.

    _params is an array/list in the following order:
        m1:  Mass of first pendulum point mass.
        m2:  Mass of second pendulum point mass.
        l1:  Length of first pendulum link.
        l2:  Length of second pendulum link.
         g:  Gravitational constant.
    Returns a list/array of kinetic energy and potential energy, respectively.

    """
    # Unpack function arguments
    q1, q2, u1, u2 = _x

    # Unpack function parameters
    m1, m2, g, l1, l2, b1, b2 = _params

    # Trigonometric functions
    c1 = cos(q1)
    c2 = cos(q2)
    s1 = sin(q1)
    s2 = sin(q2)

    # Calculate return values
    ke = m2*(2*l1*l2*c2*u1*u2 + l1**2*u1**2 + l2**2*u2**2)/2 + m1*l1**2*u1**2/2
    pe = g*(m2*(-l1*c1 - l2*(c1*c2 - s1*s2)) - l1*m1*c1)

    # Return calculated values
    return [ke, pe]

def anim(_x, _params):
    """Calculate configuration of double pendulum for purposes of animation.

    _x is an array/list of the configuration coordinates of the disc:
        q1:  Angle of first pendulum link relative to vertical (0 downwards)
        q2:  Angle of second pendulum link relative to first (0 is fully extended)
        u1:  A[1] measure number of the inertial angular velocity of the first link.
        u2:  A[1] measure number of the inertial angular velocity of the second link.

    _params is the radius of the disc.
        l1:  Length of first pendulum link.
        l2:  Length of second pendulum link.

    Output is four 3-tuples in the following order:
          P1:  Position of first point mass.
          P2:  Position of second point mass.
        B[3]:  Axis along second link of pendulum.

    """
    # Unpack function arguments
    q1, q2 = _x

    # Unpack function parameters
    m1, m2, g, l1, l2, b1, b2 = _params

    # Trigonometric functions
    c1 = cos(q1)
    c2 = cos(q2)
    s1 = sin(q1)
    s2 = sin(q2)

    # Calculate return values
    P1_1 = 0
    P1_2 = l1*s1
    P1_3 = -l1*c1
    P2_1 = 0
    P2_2 = l1*s1 - l2*(-c1*s2 - c2*s1)
    P2_3 = -l1*c1 - l2*(c1*c2 - s1*s2)
    B3_1 = 0
    B3_2 = -c1*s2 - c2*s1
    B3_3 = c1*c2 - s1*s2

    # Return calculated values
    return [[P1_1, P1_2, P1_3], [P2_1, P2_2, P2_3], [B3_1, B3_2, B3_3]]

