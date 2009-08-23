# Sun Aug 23 13:41:31 2009
from numpy import sin, cos, tan, vectorize

def f(x, t, parameter_list):
    # Unpacking the parameters
    m, g, l, b = parameter_list
    # Unpacking the states (q's and u's)
    q1, u1 = x
    s1 = sin(q1)
    # Kinematic differential equations
    q1p = u1
    # Dynamic differential equations
    u1p = -g*s1/l - b*u1/(l**2*m)
    return [q1p, u1p]

def qdot2u(q, qd, parameter_list):
    # Unpacking the parameters
    m, g, l, b = parameter_list
    # Unpacking the q's and qdots
    q1 = q
    q1p = qd
    s1 = sin(q1)
    # Kinematic differential equations
    u1 = q1p
    return [u1]

def animate(q, parameter_list):
    # Unpacking the parameters
    m, g, l, b = parameter_list
    # Unpacking the coordinates
    q1 = q
    # Trigonometric functions needed
    s1 = sin(q1)
    c1 = cos(q1)
    # Position of Points and Axis/Angle Calculations
    A3_1 = 0
    A3_2 = -s1
    A3_3 = c1
    return [A3_1, A3_2, A3_3, 0]