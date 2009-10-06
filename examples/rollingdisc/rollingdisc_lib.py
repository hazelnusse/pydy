from __future__ import division
from numpy import sin, cos, tan

def eoms(_x, t, _params):
    """Equations of motion for rolling disc.

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
         r:  Radius of disc

    """
    # Unpack function arguments
    q1, q2, q3, q4, q5, u1, u2, u3 = _x

    # Unpack function parameters
    m, g, r = _params

    # Trigonometric functions
    c1 = cos(q1)
    c2 = cos(q2)
    s1 = sin(q1)
    s2 = sin(q2)
    t2 = tan(q2)

    # Calculate return values
    q1d = u3/c2
    q2d = u1
    q3d = -t2*u3 + u2
    q4d = -r*q3d*c1
    q5d = -r*q3d*s1
    u1d = -u3**2*t2/5 + 6*u2*u3/5 + 4*g*s2/(5*r)
    u2d = -2*u1*u3/3
    u3d = -2*u1*u2 + t2*u1*u3

    # Return calculated values
    return [q1d, q2d, q3d, q4d, q5d, u1d, u2d, u3d]

def anim(_x, _params):
    """Calculate position and orientation of disc for purposes of animation.

    _x is an array/list of the configuration coordinates of the disc:
        q1:  Yaw angle (0 is aligned with N[1])
        q2:  Lean angle (0 is upright)
        q3:  Spin angle
        q4:  N[1] measure number of contact point position relative to origin
        q5:  N[2] measure number of contact point position relative to origin

    _params is the radius of the disc.

    Output is four 3-tuples in the following order:
          CO:  Postion of the center of the disc in inertial coordinates
        B[2]:  Vector normal to the plane of the disc (axis of rotation)
        C[1]:  1st body fixed unit vector in the plane of the disc
        C[3]:  3rd body fixed unit vector in the plane of the disc

    """
    # Unpack function arguments
    q1, q2, q3, q4, q5 = _x

    # Unpack function parameters
    r = _params

    # Trigonometric functions
    c1 = cos(q1)
    c2 = cos(q2)
    c3 = cos(q3)
    s1 = sin(q1)
    s2 = sin(q2)
    s3 = sin(q3)

    # Calculate return values
    CO_1 = -r*s1*s2 + q4
    CO_2 = r*c1*s2 + q5
    CO_3 = -r*c2
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

def qdot_to_u(_x):
    """Mapping from time derivates of coordinates to generalized speeds.

    _x is an array/list in the following order:
         q2:  Lean angle
        q1d:  Yaw rate
        q2d:  Lean rate
        q3d:  Spin rate

    """
    # Unpack function arguments
    q2, q1d, q2d, q3d = _x

    # Trigonometric functions
    c2 = cos(q2)
    s2 = sin(q2)

    # Calculate return values
    u1 = q2d
    u2 = q1d*s2 + q3d
    u3 = q1d*c2

    # Return calculated values
    return [u1, u2, u3]

def evals(_x, _params):
    """Eigenvalues of the linearized equations of motion about the upright
configuration with spin rate treated as a parameter.

    """
    # Unpack function arguments
    u2 = _x

    # Unpack function parameters
    g, r = _params

    # Calculate return values
    ev1 = (4*g/(5*r) - 12*u2**2/5)**(1/2)
    ev2 = 0
    ev3 = -(4*g/(5*r) - 12*u2**2/5)**(1/2)

    # Return calculated values
    return [ev1, ev2, ev3]

def critical_speed(_x):
    """Critical speed of rolling disc.

    _x is an array/list in the following order:
        g:  Gravitational constant
        r:  Radius of disc

    """
    # Unpack function arguments
    g, r = _x

    # Calculate return values
    cs = -3**(1/2)*g**(1/2)*(1/r)**(1/2)/3

    # Return calculated values
    return [cs]

def energy(_x, _params):
    """Kinetic and Potential Energy of rolling disc.

    _x is an array/list in the following order:
        q2:  Lean angle
        u1:  B[1] measure number of angular velocity (roll rate)
        u2:  B[2] measure number of angular velocity (spin rate)
        u3:  B[3] measure number of angular velocity (yaw-like rate)

    _params is an array/list in the following order:
         m:  Disc mass
         g:  Gravitational constant
         r:  Radius of disc

    Returns a list/array of kinetic energy and potential energy, respectively.

    """
    # Unpack function arguments
    q1, q2, q3, q4, q5, u1, u2, u3 = _x

    # Unpack function parameters
    m, g, r = _params

    # Trigonometric functions
    c2 = cos(q2)

    # Calculate return values
    ke = 0.125*m*r**2*u3**2 + 0.75*m*r**2*u2**2 + 0.625*m*r**2*u1**2
    pe = -g*m*r + g*m*r*c2

    # Return calculated values
    return [ke, pe]

