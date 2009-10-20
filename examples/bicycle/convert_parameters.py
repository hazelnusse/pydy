from __future__ import division
from math import sin, cos

def convert_params(_x):
    # Unpack function arguments
    rr, rrt, rf, rft, w, c, q1, xb, zb, xh, zh, mr, mb, mh, mf, IRxx, IRyy, IBxx, IByy, IBzz, IBxz, IHxx, IHyy, IHzz, IHxz, IFxx, IFyy, g = _x

    # Trigonometric functions
    c1 = cos(q1)
    s1 = sin(q1)

    # Calculate return values
    lr = (c + w)*c1 + (-rr - rrt)*s1
    ls = w*s1 + (rr + rrt - rf - rft)*c1
    lf = (rf + rft)*s1 - c*c1
    l1 = (-rr - rrt - (mb*zb + mr*(-rr - rrt))/(mb + mr))*s1 + mb*xb*c1/(mb + mr)
    l2 = -(-rr - rrt - (mb*zb + mr*(-rr - rrt))/(mb + mr))*c1 + mb*xb*s1/(mb + mr)
    l3 = (-rf - rft - (mf*(-rf - rft) + mh*zh)/(mf + mh))*s1 + mh*(xh - w)*c1/(mf + mh)
    l4 = -(-rf - rft - (mf*(-rf - rft) + mh*zh)/(mf + mh))*c1 + mh*(xh - w)*s1/(mf + mh)

    # Nested terms
    l1c = -l1
    l1d = -l1 + xb*c1 - (rr + rrt + zb)*s1
    l1e = -l3 - lf - lr + xh*c1 - (rr + rrt + zh)*s1
    l1f = -l3
    l3c = -l2
    l3d = -l2 + xb*s1 + (rr + rrt + zb)*c1
    l3e = -l4 - ls + xh*s1 + (rr + rrt + zh)*c1
    l3f = -l4

    mcd = mb + mr
    mef = mf + mh
    IC22 = IRyy
    ICD11 = IRxx - 2*IBxz*c1*s1 + IBxx*c1**2 + IBzz*s1**2 + mb*l3d**2 + mr*l3c**2
    ICD13 = IBxx*c1*s1 - IBzz*c1*s1 - l1c*l3c*mr - l1d*l3d*mb + IBxz*c1**2 - IBxz*s1**2
    ICD22 = IByy + mb*(l1d**2 + l3d**2) + mr*(l1c**2 + l3c**2)
    ICD33 = IRxx + 2*IBxz*c1*s1 + IBxx*s1**2 + IBzz*c1**2 + mb*l1d**2 + mr*l1c**2
    IEF11 = IFxx - 2*IHxz*c1*s1 + IHxx*c1**2 + IHzz*s1**2 + mf*l3f**2 + mh*l3e**2
    IEF22 = IHyy + mf*(l1f**2 + l3f**2) + mh*(l1e**2 + l3e**2)
    IEF33 = IFxx + 2*IHxz*c1*s1 + IHxx*s1**2 + IHzz*c1**2 + mf*l1f**2 + mh*l1e**2
    IEF13 = IHxx*c1*s1 - IHzz*c1*s1 - l1e*l3e*mh - l1f*l3f*mf + IHxz*c1**2 - IHxz*s1**2
    IF22 = IFyy

    # Return calculated values
    return [rr, rrt, rf, rft, lr, ls, lf, l1, l2, l3, l4, mcd, mef, IC22,
            ICD11, ICD22, ICD33, ICD13, IEF11, IEF22, IEF33, IEF13, IF22, g]

