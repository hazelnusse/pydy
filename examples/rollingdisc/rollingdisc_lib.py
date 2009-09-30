from __future__ import division
from numpy import sin, cos, tan

def eoms(_x, t, _params):
    """Equations of motion for rolling disc.

    _x is an array/list with the following ordering:
        q0:  Yaw angle (0 is aligned with N[1])
        q1:  Lean angle (0 is upright)
        q2:  Spin angle
        q3:  N[1] measure number of contact point position relative to origin
        q4:  N[2] measure number of contact point position relative to origin
        u0:  B[1] measure number of angular velocity of C in N
        u1:  B[2] measure number of angular velocity of C in N
        u2:  B[3] measure number of angular velocity of C in N

    _params is a array/list with the following ordering:
         m:  Mass of disc
         g:  Gravitational constant
         r:  Radius of disc

    """
    # Unpack function arguments
    q0, q1, q2, q3, q4, u0, u1, u2 = _x

    # Unpack function parameters
    m, g, r = _params

    # Trigonometric functions
    c0 = cos(q0)
    c1 = cos(q1)
    s0 = sin(q0)
    s1 = sin(q1)
    t1 = tan(q1)

    # Calculate return values
    q0d = u2/c1
    q1d = u0
    q2d = -t1*u2 + u1
    q3d = -r*q2d*c0
    q4d = -r*q2d*s0
    u0d = -4*(-g*m*r*s1 - 3*m*r**2*u1*u2/2 + m*r**2*u2**2*t1/4)/(5*m*r**2)
    u1d = -2*u0*u2/3
    u2d = -4*(m*r**2*u0*u1/2 - m*r**2*t1*u0*u2/4)/(m*r**2)

    # Return calculated values
    return [q0d, q1d, q2d, q3d, q4d, u0d, u1d, u2d]

def anim(_x, _params):
    """Calculate position and orientation of disc for purposes of animation.

    _x is an array/list of the configuration coordinates of the disc:
        q0:  Yaw angle (0 is aligned with N[1])
        q1:  Lean angle (0 is upright)
        q2:  Spin angle
        q3:  N[1] measure number of contact point position relative to origin
        q4:  N[2] measure number of contact point position relative to origin

    _params is the radius of the disc.

    Output is four 3-tuples in the following order:
          CO:  Postion of the center of the disc in inertial coordinates
        B[2]:  Vector normal to the plane of the disc (axis of rotation)
        C[1]:  1st body fixed unit vector in the plane of the disc
        C[3]:  3rd body fixed unit vector in the plane of the disc

    """
    # Unpack function arguments
    q0, q1, q2, q3, q4 = _x

    # Unpack function parameters
    r = _params

    # Trigonometric functions
    c0 = cos(q0)
    c1 = cos(q1)
    c2 = cos(q2)
    s0 = sin(q0)
    s1 = sin(q1)
    s2 = sin(q2)

    # Calculate return values
    CO_1 = -r*s0*s1 + q3
    CO_2 = r*c0*s1 + q4
    CO_3 = -r*c1
    B2_1 = -c1*s0
    B2_2 = c0*c1
    B2_3 = s1
    C1_1 = c0*c2 - s0*s1*s2
    C1_2 = c2*s0 + c0*s1*s2
    C1_3 = -c1*s2
    C3_1 = c0*s2 + c2*s0*s1
    C3_2 = s0*s2 - c0*c2*s1
    C3_3 = c1*c2

    # Return calculated values
    return [[CO_1, CO_2, CO_3], [B2_1, B2_2, B2_3], [C1_1, C1_2, C1_3], [C3_1, C3_2, C3_3]]

def qdot_to_u(_x):
    """Mapping from time derivates of coordinates to generalized speeds.

    _x is an array/list in the following order:
         q1:  Lean angle
        q0d:  Yaw rate
        q1d:  Lean rate
        q2d:  Spin rate

    """
    # Unpack function arguments
    q1, q0d, q1d, q2d = _x

    # Trigonometric functions
    c1 = cos(q1)
    s1 = sin(q1)

    # Calculate return values
    u0 = q1d
    u1 = q0d*s1 + q2d
    u2 = q0d*c1

    # Return calculated values
    return [u0, u1, u2]

def evals(_x, _params):
    """Eigenvalues of the linearized equations of motion about the upright
configuration with spin rate treated as a parameter.

    """
    # Unpack function arguments
    u1 = _x

    # Unpack function parameters
    g, r = _params

    # Calculate return values
    ev0 = (4*g/(5*r) - 12*u1**2/5)**(1/2)
    ev1 = -(4*g/(5*r) - 12*u1**2/5)**(1/2)
    ev2 = 0

    # Return calculated values
    return [ev0, ev1, ev2]

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
        q1:  Lean angle
        u0:  B[1] measure number of angular velocity (roll rate)
        u1:  B[2] measure number of angular velocity (spin rate)
        u3:  B[3] measure number of angular velocity (yaw-like rate)

    _params is an array/list in the following order:
         m:  Disc mass
         g:  Gravitational constant
         r:  Radius of disc

    Returns a list/array of kinetic energy and potential energy, respectively.

    """
    # Unpack function arguments
    q0, q1, q2, q3, q4, u0, u1, u2 = _x

    # Unpack function parameters
    m, g, r = _params

    # Trigonometric functions
    c1 = cos(q1)

    # Calculate return values
    ke = 0.125*m*r**2*u2**2 + 0.75*m*r**2*u1**2 + 0.625*m*r**2*u0**2
    pe = -g*m*r + g*m*r*c1

    # Return calculated values
    return [ke, pe]

