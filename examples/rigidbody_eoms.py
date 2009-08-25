# Tue Aug 25 15:16:01 2009
from numpy import sin, cos, tan, vectorize

def f(x, t, parameter_list):
    # Unpacking the parameters
    m, g, I11, I22, I33 = parameter_list
    # Unpacking the states (q's and u's)
    q1, q2, q3, q4, q5, q6, u1, u2, u3, u4, u5, u6 = x
    c3 = cos(q3)
    c2 = cos(q2)
    s2 = sin(q2)
    s3 = sin(q3)
    # Kinematic differential equations
    q1p = c3*u3/c2 - s3*u1/c2
    q2p = c3*u1 + s3*u3
    q3p = s2*s3*u1/c2 - c3*s2*u3/c2 + u2
    q4p = u4
    q5p = u5
    q6p = u6
    # Dynamic differential equations
    u1p = (I22 - I33)*u2*u3/I11
    u2p = (I33 - I11)*u1*u3/I22
    u3p = (I11 - I22)*u1*u2/I33
    u4p = 0
    u5p = 0
    u6p = g
    return [q1p, q2p, q3p, q4p, q5p, q6p, u1p, u2p, u3p, u4p, u5p, u6p]

def qdot2u(q, qd, parameter_list):
    # Unpacking the parameters
    m, g, I11, I22, I33 = parameter_list
    # Unpacking the q's and qdots
    q1, q2, q3, q4, q5, q6 = q
    q1p, q2p, q3p, q4p, q5p, q6p = qd
    c3 = cos(q3)
    c2 = cos(q2)
    s2 = sin(q2)
    s3 = sin(q3)
    # Kinematic differential equations
    u1 = q2p*c3 - q1p*c2*s3
    u2 = q3p + q1p*s2
    u3 = q2p*s3 + q1p*c2*c3
    u4 = q4p
    u5 = q5p
    u6 = q6p
    return [u1, u2, u3, u4, u5, u6]

def animate(q, parameter_list):
    # Unpacking the parameters
    m, g, I11, I22, I33 = parameter_list
    # Unpacking the coordinates
    q1, q2, q3, q4, q5, q6 = q
    # Trigonometric functions needed
    s1 = sin(q1)
    s2 = sin(q2)
    c1 = cos(q1)
    c2 = cos(q2)
    # Position of Points and Axis/Angle Calculations
    p_NO_CO_1 = q4
    p_NO_CO_2 = q5
    p_NO_CO_3 = q6
    C2_1 = -c2*s1
    C2_2 = c1*c2
    C2_3 = s2
    return [p_NO_CO_1, p_NO_CO_2, p_NO_CO_3], [C2_1, C2_2, C2_3, q3]