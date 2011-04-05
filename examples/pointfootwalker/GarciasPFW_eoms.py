# Sun Aug 23 13:12:56 2009
from numpy import sin, cos, tan, vectorize

def f(x, t, parameter_list):
    # Unpacking the parameters
    Mh, mf, g, L, q3 = parameter_list
    # Unpacking the states (q's and u's)
    q1, q2, u1, u2 = x
    s1 = sin(q1)
    cos(q3) = cos(q3)
    c2 = cos(q2)
    sin(q3) = sin(q3)
    s2 = sin(q2)
    c1 = cos(q1)
    # Kinematic differential equations
    q1p = u1
    q2p = -u2 + u1
    # Dynamic differential equations
    u1p = -L**4*mf**2*u2**2*s2/(-L**4*mf**2 + L**4*mf**2*c2**2 - Mh*mf*L**4) + g*L**3*mf**2*cos(q3)*c1/(-L**4*mf**2 + L**4*mf**2*c2**2 - Mh*mf*L**4) + L**4*mf**2*u1**2*c2*s2/(-L**4*mf**2 + L**4*mf**2*c2**2 - Mh*mf*L**4) - g*L**3*mf**2*sin(q3)*s1/(-L**4*mf**2 + L**4*mf**2*c2**2 - Mh*mf*L**4) + Mh*g*mf*L**3*cos(q3)*c1/(-L**4*mf**2 + L**4*mf**2*c2**2 - Mh*mf*L**4) + g*L**3*mf**2*c2**2*sin(q3)*s1/(-L**4*mf**2 + L**4*mf**2*c2**2 - Mh*mf*L**4) - Mh*g*mf*L**3*sin(q3)*s1/(-L**4*mf**2 + L**4*mf**2*c2**2 - Mh*mf*L**4) - g*L**3*mf**2*c2**2*cos(q3)*c1/(-L**4*mf**2 + L**4*mf**2*c2**2 - Mh*mf*L**4) - g*L**3*mf**2*cos(q3)*c2*s1*s2/(-L**4*mf**2 + L**4*mf**2*c2**2 - Mh*mf*L**4) - g*L**3*mf**2*c1*c2*sin(q3)*s2/(-L**4*mf**2 + L**4*mf**2*c2**2 - Mh*mf*L**4)
    u2p = L**4*mf**2*u1**2*s2/(-L**4*mf**2 + L**4*mf**2*c2**2 - Mh*mf*L**4) + Mh*mf*L**4*u1**2*s2/(-L**4*mf**2 + L**4*mf**2*c2**2 - Mh*mf*L**4) - L**4*mf**2*u2**2*c2*s2/(-L**4*mf**2 + L**4*mf**2*c2**2 - Mh*mf*L**4) - g*L**3*mf**2*cos(q3)*s1*s2/(-L**4*mf**2 + L**4*mf**2*c2**2 - Mh*mf*L**4) - g*L**3*mf**2*c1*sin(q3)*s2/(-L**4*mf**2 + L**4*mf**2*c2**2 - Mh*mf*L**4) - Mh*g*mf*L**3*cos(q3)*s1*s2/(-L**4*mf**2 + L**4*mf**2*c2**2 - Mh*mf*L**4) - Mh*g*mf*L**3*c1*sin(q3)*s2/(-L**4*mf**2 + L**4*mf**2*c2**2 - Mh*mf*L**4)
    return [q1p, q2p, u1p, u2p]

def qdot2u(q, qd, parameter_list):
    # Unpacking the parameters
    Mh, mf, g, L, q3 = parameter_list
    # Unpacking the q's and qdots
    q1, q2 = q
    q1p, q2p = qd
    s1 = sin(q1)
    cos(q3) = cos(q3)
    c2 = cos(q2)
    sin(q3) = sin(q3)
    s2 = sin(q2)
    c1 = cos(q1)
    # Kinematic differential equations
    u1 = q1p
    u2 = q1p - q2p
    return [u1, u2]