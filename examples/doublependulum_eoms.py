# Sun Aug 23 15:08:51 2009
from numpy import sin, cos, tan, vectorize

def f(x, t, parameter_list):
    # Unpacking the parameters
    m1, m2, g, l1, l2, b1, b2 = parameter_list
    # Unpacking the states (q's and u's)
    q1, q2, u1, u2 = x
    c2 = cos(q2)
    s1 = sin(q1)
    s2 = sin(q2)
    c1 = cos(q1)
    # Kinematic differential equations
    q1p = u1
    q2p = -u1 + u2
    # Dynamic differential equations
    u1p = m2*l2**2*(b1*u1 + b2*u1 - b2*u2 + g*l1*m1*s1 + g*l1*m2*s1 - l1*l2*m2*u2**2*s2)/(-l1**2*l2**2*m2**2 + l1**2*l2**2*m2**2*c2**2 - m1*m2*l1**2*l2**2) + l1*l2*m2*(b2*u2 - b2*u1 + g*l2*m2*(c1*s2 + c2*s1) + l1*l2*m2*u1**2*s2)*c2/(l1**2*l2**2*m2**2 + m1*m2*l1**2*l2**2 - l1**2*l2**2*m2**2*c2**2)
    u2p = (m1*l1**2 + m2*l1**2)*(b2*u2 - b2*u1 + g*l2*m2*(c1*s2 + c2*s1) + l1*l2*m2*u1**2*s2)/(-l1**2*l2**2*m2**2 + l1**2*l2**2*m2**2*c2**2 - m1*m2*l1**2*l2**2) + l1*l2*m2*(b1*u1 + b2*u1 - b2*u2 + g*l1*m1*s1 + g*l1*m2*s1 - l1*l2*m2*u2**2*s2)*c2/(l1**2*l2**2*m2**2 + m1*m2*l1**2*l2**2 - l1**2*l2**2*m2**2*c2**2)
    return [q1p, q2p, u1p, u2p]

def qdot2u(q, qd, parameter_list):
    # Unpacking the parameters
    m1, m2, g, l1, l2, b1, b2 = parameter_list
    # Unpacking the q's and qdots
    q1, q2 = q
    q1p, q2p = qd
    c2 = cos(q2)
    s1 = sin(q1)
    s2 = sin(q2)
    c1 = cos(q1)
    # Kinematic differential equations
    u1 = q1p
    u2 = q1p + q2p
    return [u1, u2]

def animate(q, parameter_list):
    # Unpacking the parameters
    m1, m2, g, l1, l2, b1, b2 = parameter_list
    # Unpacking the coordinates
    q1, q2 = q
    # Trigonometric functions needed
    c2 = cos(q2)
    s1 = sin(q1)
    s2 = sin(q2)
    c1 = cos(q1)
    # Position of Points and Axis/Angle Calculations
    A3_1 = 0
    A3_2 = -s1
    A3_3 = c1
    B3_1 = 0
    B3_2 = -c1*s2 - c2*s1
    B3_3 = c1*c2 - s1*s2
    return [A3_1, A3_2, A3_3, 0], [B3_1, B3_2, B3_3, 0]