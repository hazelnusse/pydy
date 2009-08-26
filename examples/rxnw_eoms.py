# Tue Aug 25 17:42:58 2009
from numpy import sin, cos, tan, vectorize

def f(x, t, parameter_list):
    # Unpacking the parameters
    m1, m2, I1, I2, J1, J2, l, M = parameter_list
    # Unpacking the states (q's and u's)
    q1, q2, q3, q4, q5, q6, q7, q8, u1, u2, u3, u4, u5, u6, u7, u8 = x
    c3 = cos(q3)
    c2 = cos(q2)
    s1 = sin(q1)
    s2 = sin(q2)
    c1 = cos(q1)
    s3 = sin(q3)
    # Kinematic differential equations
    q1p = c3*u1/c2 - s3*u2/c2
    q2p = c3*u2 + s3*u1
    q3p = s2*s3*u2/c2 - c3*s2*u1/c2 + u3
    q4p = u4
    q5p = u5
    q6p = s2*u8 + c2*c3*u6 - c2*s3*u7
    q7p = (c1*c3 - s1*s2*s3)*u7 + (c1*s3 + c3*s1*s2)*u6 - c2*s1*u8
    q8p = (c3*s1 + c1*s2*s3)*u7 + (s1*s3 - c1*c3*s2)*u6 + c1*c2*u8
    # Dynamic differential equations
    u1p = (-J2*u3*u4 - J2*u3*u5 - (I1 - J1)*u2*u3 - 2*(J2 - I2)*u2*u3 + 2*m2*l**2*u2*u3)/(-I1 - 2*I2 - 2*m2*l**2)
    u2p = (J1 - I1)*u1*u3/I1 - 2*M/I1
    u3p = (J2*u1*u4 + J2*u1*u5 - 2*(I2 - J2)*u1*u2 - 2*m2*l**2*u1*u2)/(-J1 - 2*I2 - 2*m2*l**2)
    u4p = -M*(I1*J1 + J1*J2 + 2*I1*I2 + 2*I2*J2 + 2*I1*m2*l**2 + 2*J2*m2*l**2)/(-I1*J1*J2 - 2*I1*I2*J2 - 2*I1*J2*m2*l**2) - (J1 - I1)*u1*u3/I1 + M/I1
    u5p = -M*(I1*J1 + J1*J2 + 2*I1*I2 + 2*I2*J2 + 2*I1*m2*l**2 + 2*J2*m2*l**2)/(-I1*J1*J2 - 2*I1*I2*J2 - 2*I1*J2*m2*l**2) - (J1 - I1)*u1*u3/I1 + M/I1
    u6p = (m1*u2*u8 - m1*u3*u7 - 2*m2*u3*u7 + 2*m2*u2*u8)/(-m1 - 2*m2)
    u7p = (m1*u3*u6 - m1*u1*u8 - 2*m2*u1*u8 + 2*m2*u3*u6)/(-m1 - 2*m2)
    u8p = (m1*u1*u7 - m1*u2*u6 - 2*m2*u2*u6 + 2*m2*u1*u7)/(-m1 - 2*m2)
    return [q1p, q2p, q3p, q4p, q5p, q6p, q7p, q8p, u1p, u2p, u3p, u4p, u5p, u6p, u7p, u8p]

def qdot2u(q, qd, parameter_list):
    # Unpacking the parameters
    m1, m2, I1, I2, J1, J2, l, M = parameter_list
    # Unpacking the q's and qdots
    q1, q2, q3, q4, q5, q6, q7, q8 = q
    q1p, q2p, q3p, q4p, q5p, q6p, q7p, q8p = qd
    c3 = cos(q3)
    c2 = cos(q2)
    s1 = sin(q1)
    s2 = sin(q2)
    c1 = cos(q1)
    s3 = sin(q3)
    # Kinematic differential equations
    u1 = q2p*s3 + q1p*c2*c3
    u2 = q2p*c3 - q1p*c2*s3
    u3 = q3p + q1p*s2
    u4 = q4p
    u5 = q5p
    u6 = q7p*(c1*s3 + c3*s1*s2) + q8p*(s1*s3 - c1*c3*s2) + q6p*c2*c3
    u7 = q7p*(c1*c3 - s1*s2*s3) + q8p*(c3*s1 + c1*s2*s3) - q6p*c2*s3
    u8 = q6p*s2 + q8p*c1*c2 - q7p*c2*s1
    return [u1, u2, u3, u4, u5, u6, u7, u8]

def animate(q, parameter_list):
    # Unpacking the parameters
    m1, m2, I1, I2, J1, J2, l, M = parameter_list
    # Unpacking the coordinates
    q1, q2, q3, q4, q5, q6, q7, q8 = q
    # Trigonometric functions needed
    c3 = cos(q3)
    s5 = sin(q5)
    c2 = cos(q2)
    c5 = cos(q5)
    s1 = sin(q1)
    s2 = sin(q2)
    c1 = cos(q1)
    c4 = cos(q4)
    s3 = sin(q3)
    s4 = sin(q4)
    # Position of Points and Axis/Angle Calculations
    p_NO_AO_1 = q6
    p_NO_AO_2 = q7
    p_NO_AO_3 = q8
    p_NO_BO_1 = l*c2*s3 + q6
    p_NO_BO_2 = -l*(c1*c3 - s1*s2*s3) + q7
    p_NO_BO_3 = -l*(c3*s1 + c1*s2*s3) + q8
    p_NO_CO_1 = -l*c2*s3 + q6
    p_NO_CO_2 = l*(c1*c3 - s1*s2*s3) + q7
    p_NO_CO_3 = l*(c3*s1 + c1*s2*s3) + q8
    A1_1 = c2*c3
    A1_2 = c1*s3 + c3*s1*s2
    A1_3 = s1*s3 - c1*c3*s2
    A2_1 = -c2*s3
    A2_2 = c1*c3 - s1*s2*s3
    A2_3 = c3*s1 + c1*s2*s3
    A3_1 = s2
    A3_2 = -c2*s1
    A3_3 = c1*c2
    B1_1 = -s2*s4 + c2*c3*c4
    B1_2 = c1*c4*s3 + c2*s1*s4 + c3*c4*s1*s2
    B1_3 = c4*s1*s3 - c1*c2*s4 - c1*c3*c4*s2
    B3_1 = c4*s2 + c2*c3*s4
    B3_2 = c1*s3*s4 - c2*c4*s1 + c3*s1*s2*s4
    B3_3 = c1*c2*c4 + s1*s3*s4 - c1*c3*s2*s4
    C1_1 = -s2*s5 + c2*c3*c5
    C1_2 = c1*c5*s3 + c2*s1*s5 + c3*c5*s1*s2
    C1_3 = c5*s1*s3 - c1*c2*s5 - c1*c3*c5*s2
    C3_1 = c5*s2 + c2*c3*s5
    C3_2 = c1*s3*s5 - c2*c5*s1 + c3*s1*s2*s5
    C3_3 = c1*c2*c5 + s1*s3*s5 - c1*c3*s2*s5
    return [p_NO_AO_1, p_NO_AO_2, p_NO_AO_3], [p_NO_BO_1, p_NO_BO_2, p_NO_BO_3], [p_NO_CO_1, p_NO_CO_2, p_NO_CO_3], [A1_1, A1_2, A1_3, 0], [A2_1, A2_2, A2_3, 0], [A3_1, A3_2, A3_3, 0], [B1_1, B1_2, B1_3, 0], [B3_1, B3_2, B3_3, 0], [C1_1, C1_2, C1_3, 0], [C3_1, C3_2, C3_3, 0]