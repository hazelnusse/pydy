# Wed Aug 12 17:26:54 2009
from numpy import sin, cos, tan, vectorize

def f(x, t, parameter_list):
    # Unpacking the parameters
    m, g, r = parameter_list
    # Unpacking the states (q's and u's)
    q1, q2, q3, q4, q5, u1, u2, u3 = x
    s1 = sin(q1)
    s2 = sin(q2)
    c1 = cos(q1)
    c2 = cos(q2)
    # Kinematic differential equations
    q1p = u3/c2
    q2p = u1
    q3p = -s2*u3/c2 + u2
    # Dependent differential equations
    q4p = -q3p*r*c1
    q5p = -q3p*r*s1
    # Dynamic differential equations
    u1p = 2*u2*u3 - u3**2*s2/c2 + 4*g*s2/(5*r)
    u2p = -2*u1*u3/3
    u3p = -2*u1*u2 + s2*u1*u3/c2
    return [q1p, q2p, q3p, q4p, q5p, u1p, u2p, u3p]

def qdot2u(q, qd, parameter_list):
    # Unpacking the parameters
    m, g, r = parameter_list
    # Unpacking the q's and qdots
    q1, q2, q3, q4, q5 = q
    q1p, q2p, q3p = qd
    s1 = sin(q1)
    s2 = sin(q2)
    c1 = cos(q1)
    c2 = cos(q2)
    # Kinematic differential equations
    u1 = q2p
    u2 = q3p + q1p*s2
    u3 = q1p*c2
    return [u1, u2, u3]

def animate(q, parameter_list):
    # Unpacking the parameters
    m, g, r = parameter_list
    # Unpacking the coordinates
    q1, q2, q3, q4, q5 = q
    # Trigonometric functions needed
    c3 = cos(q3)
    c2 = cos(q2)
    s1 = sin(q1)
    s2 = sin(q2)
    c1 = cos(q1)
    s3 = sin(q3)
    # Position of Points and Axis/Angle Calculations
    p_N1_CO_1 = -r*s1*s2 + q4
    p_N1_CO_2 = r*c1*s2 + q5
    p_N1_CO_3 = -r*c2
    B2_1 = -c2*s1
    B2_2 = c1*c2
    B2_3 = s2
    C1_1 = c1*c3 - s1*s2*s3
    C1_2 = c3*s1 + c1*s2*s3
    C1_3 = -c2*s3
    C3_1 = c1*s3 + c3*s1*s2
    C3_2 = s1*s3 - c1*c3*s2
    C3_3 = c2*c3
    return [p_N1_CO_1, p_N1_CO_2, p_N1_CO_3], [B2_1, B2_2, B2_3, q3], [C1_1, C1_2, C1_3, 0], [C3_1, C3_2, C3_3, 0]