from __future__ import division
from math import sin, cos

def eoms(_x, t, _params):
    """Point mass pendulum equations of motion.

    _x is an array/list in the following order:
        q1:  Angle of pendulum link relative to vertical (0 downwards)
        u1:  A[1] measure number of the inertial angular velocity of the first link.

    _params is an array/list in the following order:
        m:  Mass of first pendulum point mass.
        l:  Length of first pendulum link.
        g:  Gravitational constant.
        b:  Damping coefficient at hinge.

    """
    # Unpack function arguments
    q1, u1 = _x

    # Unpack function parameters
    m, g, l, b = _params

    # Trigonometric functions
    s1 = sin(q1)

    # Calculate return values
    q1d = u1
    u1d = -g*s1/l - b*u1/(l**2*m)

    # Return calculated values
    return [q1d, u1d]

def energy(_x, _params):
    """Kinetic and Potential Energy of point mass pendulum.

    _x is an array/list in the following order:
        q1:  Angle of first pendulum link relative to vertical (0 downwards)
        u1:  A[1] measure number of the inertial angular velocity of the first link.

    _params is an array/list in the following order:
        m:  Mass of first pendulum point mass.
        l:  Length of first pendulum link.
        g:  Gravitational constant.
    Returns a list/array of kinetic energy and potential energy, respectively.

    """
    # Unpack function arguments
    q1, u1 = _x

    # Unpack function parameters
    m, g, l, b = _params

    # Trigonometric functions
    c1 = cos(q1)

    # Calculate return values
    ke = m*l**2*u1**2/2
    pe = g*l*m*(1 - c1)

    # Return calculated values
    return [ke, pe]

def anim(_x, _params):
    """Calculate configuration of pendulum for purposes of animation.

    _x is an array/list of the configuration coordinates of the disc:
        q1:  Angle of first pendulum link relative to vertical (0 downwards)
        u1:  A[1] measure number of the inertial angular velocity of the first link.

    _params is the radius of the disc.
        m:  Mass of first pendulum point mass.
        l:  Length of first pendulum link.
        g:  Gravitational constant.

    Output is:
          P:  Position of first point mass.

    """
    # Unpack function arguments
    q1 = _x

    # Unpack function parameters
    m, g, l, b = _params

    # Trigonometric functions
    c1 = cos(q1)
    s1 = sin(q1)

    # Calculate return values
    P_1 = 0
    P_2 = l*s1
    P_3 = -l*c1

    # Return calculated values
    return [P_1, P_2, P_3]

