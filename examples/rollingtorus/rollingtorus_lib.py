from __future__ import division
from numpy import sin, cos, tan

def eoms(_x, t, _params):
    """Equations of motion for rolling torus.

    _x is an array/list with the following ordering:
        q1:  Yaw angle (0 is aligned with N[1])
        q2:  Lean angle (0 is upright)
        q3:  Spin angle
        q4:  N[1] measure number of contact point position relative to origin
        q5:  N[2] measure number of contact point position relative to origin
        u1:  B[1] measure number of angular velocity of C in N
        u2:  B[2] measure number of angular velocity of C in N
        u3:  B[3] measure number of angular velocity of C in N

    _params is a array/list with the following ordering:
         m:  Mass of disc
         g:  Gravitational constant
         r1: Major radius of torus
         r2: Minor radius of torus

    """
    # Unpack function arguments
    q1, q2, q3, q4, q5, u1, u2, u3 = _x

    # Unpack function parameters
    m, g, r1, r2, I, J = _params

    # Trigonometric functions
    c1 = cos(q1)
    c2 = cos(q2)
    s1 = sin(q1)
    s2 = sin(q2)
    t2 = tan(q2)

    # Nested terms
    _M00 = -2*r1*r2*c2 - 13*r2**2/8 - 3*r1**2/2
    _M11 = -2*r1*r2*c2 - 2*r1**2 - 3*r2**2/4 - r2**2*c2**2
    _M12 = r1*r2*s2 + r2**2*c2*s2
    _M22 = -5*r2**2/8 - r1**2/2 - r2**2*s2**2
    _rhs0 = -g*r1*s2 - 2*r1**2*u2*u3 - 3*r2**2*u2*u3/4 + r1*r2*u3**2*s2 + r2**2*u3**2*c2*s2 - r1*r2*u1**2*s2 - r2**2*c2**2*u2*u3 - 2*r1*r2*c2*u2*u3
    _rhs1 = r1**2*u1*u3 + r2**2*u1*u3 - r2**2*c2**2*u1*u3 + r1*r2*u1*u3/c2 - r1*r2*s2*u1*u2 - r2**2*c2*s2*u1*u2
    _rhs2 = r1**2*u1*u2 + 3*r2**2*u1*u2/4 + r2**2*s2**2*u1*u2 + r2**2*c2*s2*u1*u3 - r1*r2*s2*u1*u3 - r2**2*s2*u1*u3/c2

    # Calculate return values
    q1d = u3/c2
    q2d = u1
    q3d = -t2*u3 + u2
    q4d = (-r1*c1 - r2*c1*c2)*q3d - r2*q2d*s1
    q5d = (-r1*s1 - r2*c2*s1)*q3d + r2*q2d*c1
    u1d = _rhs0/_M00
    u2d = (_M22*_rhs1 - _M12*_rhs2)/(_M11*_M22 - _M12**2)
    u3d = (_M11*_rhs2 - _M12*_rhs1)/(_M11*_M22 - _M12**2)

    # Return calculated values
    return [q1d, q2d, q3d, q4d, q5d, u1d, u2d, u3d]

def anim(_x, _params):
    """Calculate position and orientation of torus for purposes of animation.

    _x is an array/list of the configuration coordinates of the disc:
        q1:  Yaw angle (0 is aligned with N[1])
        q2:  Lean angle (0 is upright)
        q3:  Spin angle
        q4:  N[1] measure number of contact point position relative to origin
        q5:  N[2] measure number of contact point position relative to origin

    _params is the radius of the torus.

    Output is four 3-tuples in the following order:
          CO:  Postion of the center of the disc in inertial coordinates
        B[2]:  Vector normal to the plane of the disc (axis of rotation)
        C[1]:  1st body fixed unit vector in the plane of the disc
        C[3]:  3rd body fixed unit vector in the plane of the disc

    """
    # Unpack function arguments
    q1, q2, q3, q4, q5 = _x

    # Unpack function parameters
    m, g, r1, r2, I, J = _params

    # Trigonometric functions
    c1 = cos(q1)
    c2 = cos(q2)
    c3 = cos(q3)
    s1 = sin(q1)
    s2 = sin(q2)
    s3 = sin(q3)

    # Calculate return values
    CO_1 = -r1*s1*s2 + q4
    CO_2 = r1*c1*s2 + q5
    CO_3 = -r2 - r1*c2
    B2_1 = -c2*s1
    B2_2 = c1*c2
    B2_3 = s2
    C1_1 = c1*c3 - s1*s2*s3
    C1_2 = c3*s1 + c1*s2*s3
    C1_3 = -c2*s3
    C3_1 = c1*s3 + c3*s1*s2
    C3_2 = s1*s3 - c1*c3*s2
    C3_3 = c2*c3

    # Return calculated values
    return [[CO_1, CO_2, CO_3], [B2_1, B2_2, B2_3], [C1_1, C1_2, C1_3], [C3_1, C3_2, C3_3]]

def energy(_x, _params):
    """Kinetic and Potential Energy of rolling torus.

    _x is an array/list in the following order:
        q2:  Lean angle
        u1:  B[1] measure number of angular velocity (roll rate)
        u2:  B[2] measure number of angular velocity (spin rate)
        u3:  B[3] measure number of angular velocity (yaw-like rate)

    _params is an array/list in the following order:
         m:  Disc mass
         g:  Gravitational constant
         r1: Major radius of torus
         r2: Minor radius of torus

    Returns a list/array of kinetic energy and potential energy, respectively.

    """
    # Unpack function arguments
    q1, q2, q3, q4, q5, u1, u2, u3 = _x

    # Unpack function parameters
    m, g, r1, r2, I, J = _params

    # Trigonometric functions
    c2 = cos(q2)
    s2 = sin(q2)

    # Calculate return values
    ke = m*r1*r2*u1**2*c2 + m*r1*r2*u2**2*c2 - 1.0*m*r1*r2*s2*u2*u3 - 1.0*m*r2**2*c2*s2*u2*u3 + m*r1**2*u2**2 + 0.25*m*r1**2*u3**2 + 0.375*m*r2**2*u2**2 + 0.75*m*r1**2*u1**2 + 0.3125*m*r2**2*u3**2 + 0.8125*m*r2**2*u1**2 + 0.5*m*r2**2*c2**2*u2**2 + 0.5*m*r2**2*s2**2*u3**2
    pe = -g*m*(r1 + r2) - g*m*(-r2 - r1*c2)

    # Return calculated values
    return [ke, pe]

