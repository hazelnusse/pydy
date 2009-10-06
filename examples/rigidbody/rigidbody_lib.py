from __future__ import division
from math import sin, cos, tan

def eoms(_x, t, _params):
    """Rigidy body equations of motion.

    _x is an array/list in the following order:
        q1:  Yaw           q2:  Lean   |-(Euler 3-1-2 angles used to orient A
        q3:  Pitch /
        q4:  N[1] displacement of mass center.
        q5:  N[2] displacement of mass center.
        q6:  N[3] displacement of mass center.
        u1:  A[1] measure number of angular velocity
        u2:  A[2] measure number of angular velocity
        u3:  A[3] measure number of angular velocity
        u4:  N[1] velocity of mass center.
        u5:  N[2] velocity of mass center.
        u6:  N[3] velocity of mass center.

    _params is an array/list in the following order:
        m:  Mass of first pendulum point mass.
        g:  Gravitational constant.
      I11:  Principal moment of inertia about A[1]
      I22:  Principal moment of inertia about A[2]
      I33:  Principal moment of inertia about A[3]

    """
    # Unpack function arguments
    q1, q2, q3, q4, q5, q6, u1, u2, u3, u4, u5, u6 = _x

    # Unpack function parameters
    m, g, I11, I22, I33 = _params

    # Trigonometric functions
    c2 = cos(q2)
    c3 = cos(q3)
    s3 = sin(q3)
    t2 = tan(q2)

    # Calculate return values
    q1d = c3*u3/c2 - s3*u1/c2
    q2d = c3*u1 + s3*u3
    q3d = s3*t2*u1 - c3*t2*u3 + u2
    q4d = u4
    q5d = u5
    q6d = u6
    u1d = (I22 - I33)*u2*u3/I11
    u2d = (I33 - I11)*u1*u3/I22
    u3d = -(I22 - I11)*u1*u2/I33
    u4d = 0
    u5d = 0
    u6d = g

    # Return calculated values
    return [q1d, q2d, q3d, q4d, q5d, q6d, u1d, u2d, u3d, u4d, u5d, u6d]

def energy(_x, _params):
    """Kinetic and Potential Energy of rigid body.

    _x is an array/list in the following order:
        q1:  Yaw           q2:  Lean   |-(Euler 3-1-2 angles used to orient A
        q3:  Pitch /
        q4:  N[1] displacement of mass center.
        q5:  N[2] displacement of mass center.
        q6:  N[3] displacement of mass center.
        u1:  A[1] measure number of angular velocity
        u2:  A[2] measure number of angular velocity
        u3:  A[3] measure number of angular velocity
        u4:  N[1] velocity of mass center.
        u5:  N[2] velocity of mass center.
        u6:  N[3] velocity of mass center.

    _params is an array/list of:
        m:  Mass of first pendulum point mass.
        g:  Gravitational constant.
      I11:  Principal moment of inertia about A[1]
      I22:  Principal moment of inertia about A[2]
      I33:  Principal moment of inertia about A[3]

    Returns a list/array of kinetic energy and potential energy, respectively.

    """
    # Unpack function arguments
    q1, q2, q3, q4, q5, q6, u1, u2, u3, u4, u5, u6 = _x

    # Unpack function parameters
    m, g, I11, I22, I33 = _params

    # Calculate return values
    ke = I11*u1**2/2 + I22*u2**2/2 + I33*u3**2/2 + m*u4**2/2 + m*u5**2/2 + m*u6**2/2
    pe = -g*m*q6

    # Return calculated values
    return [ke, pe]

def anim(_x, _params):
    """Calculate configuration of pendulum for purposes of animation.

    _x is an array/list of the configuration coordinates of the disc:

    _params is an array/list of:
        m:  Mass of first pendulum point mass.
        g:  Gravitational constant.
      I11:  Principal moment of inertia about A[1]
      I22:  Principal moment of inertia about A[2]
      I33:  Principal moment of inertia about A[3]

    Output is:
          AO:  Position of mass center.
        A[1]:  Body fixed A[1] unit vector
        A[3]:  Body fixed A[3] unit vector

    """
    # Unpack function arguments
    q1, q2, q3, q4, q5, q6 = _x

    # Unpack function parameters
    m, g, I11, I22, I33 = _params

    # Trigonometric functions
    c1 = cos(q1)
    c2 = cos(q2)
    c3 = cos(q3)
    s1 = sin(q1)
    s2 = sin(q2)
    s3 = sin(q3)

    # Calculate return values
    AO_1 = q4
    AO_2 = q5
    AO_3 = q6
    A1_1 = c1*c3 - s1*s2*s3
    A1_2 = c3*s1 + c1*s2*s3
    A1_3 = -c2*s3
    A3_1 = c1*s3 + c3*s1*s2
    A3_2 = s1*s3 - c1*c3*s2
    A3_3 = c2*c3

    # Return calculated values
    return [[AO_1, AO_2, AO_3], [A1_1, A1_2, A1_3], [A3_1, A3_2, A3_3]]

