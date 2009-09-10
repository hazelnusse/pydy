from sympy import sin, cos, symbols, Matrix
from pydy import ReferenceFrame, Vector, dot
import numpy as np
import time
from scipy.optimize import fsolve

def cythonit(pyxfile, html=False):
    import os
    os.system('cython ' + pyxfile + '.pyx -a')
    if html:
        os.system('firefox ' + pyxfile + '.html')

def compileit(cfile, outputfile):
    import os
    os.system('gcc -shared -pthread -fPIC -fwrapv -O2 -Wall\
    -fno-strict-aliasing -I/usr/include/python2.6 -o ' + outputfile + ' ' +\
    cfile)

def generate_kf_module():
    # Lean, Pitch, Steer
    l, p, s = symbols('l p s')

    # Frame geometric parameters
    rr, rrt, rf, rft, lr, lf, ls = symbols('rr rrt rf rft lr lf ls')
    A = ReferenceFrame('A')
    B = A.rotate('B', 1, l)
    D = B.rotate('D', 2, p)
    E = D.rotate('E', 3, s)

    # Vector in the plane of the front wheel, pointed towards the ground
    g = Vector(A[3] - dot(E[2], A[3])*E[2]).normalized

    # Holonomic constraint
    f0 = dot(A[3], -rrt*A[3] - rr*B[3] + lr*D[1] + ls*D[3] + lf*E[1] + rf*g + rft*A[3])
    f1 = f0.diff(p)

    # Vector valued function
    F = Matrix([f0, f1])
    X = Matrix([l, p])
    J = F.jacobian(X)

    # Generate string representations of each function
    f0_s = str(f0).replace('sin', 'np.sin').replace('cos', 'np.cos')
    f1_s = str(f1).replace('sin', 'np.sin').replace('cos', 'np.cos')

    J00_s = str(J[0,0]).replace('sin', 'np.sin').replace('cos', 'np.cos')
    J01_s = str(J[0,1]).replace('sin', 'np.sin').replace('cos', 'np.cos')
    J10_s = str(J[1,0]).replace('sin', 'np.sin').replace('cos', 'np.cos')
    J11_s = str(J[1,1]).replace('sin', 'np.sin').replace('cos', 'np.cos')

    # Decide on what type of float to use here
    dtype = "float64"
    fh = open('kinematic_feasibility.pyx', 'w')

    f_s  = "from __future__ import division\n"
    f_s += "import numpy as np\n"
    f_s += "cimport numpy as np\n"
    f_s += "DTYPE = np." + dtype + "\n"
    f_s += "ctypedef np." + dtype + "_t DTYPE_t\n"
    f_s += "cimport cython\n"
    f_s += "@cython.boundscheck(False)\n"
    f_s += "def f(np.ndarray[DTYPE_t, ndim=1] x, np.ndarray[DTYPE_t, ndim=1] params):\n"
    f_s += '''    """Computes the holonomic constraint and its partial derivative
                  with respect to pitch.

                  x:  Numpy array of lean and steer, length 2

                  params:  Numpy array of parameters, length 8
                          in the following order:
                         rr:  Rear wheel radius.
                        rrt:  Rear wheel tire radius.
                         rf:  Front wheel radius.
                        rft:  Front wheel tire radius.
                         lr:  Rear wheel center perpendicular distance from steer axis.
                         lf:  Front wheel center perpendicular distance from steer axis.
                         ls:  Steer axis offset.
                          s:  Steer angle.  (treated as a parameter)

                 Returns a numpy array of the value of the holonomic constraint
                 in the first entry, and the partial derivative of the holonomic
                 constraint with respect to pitch in the second entry.  The
                 zeros of this function occur on the boundary of the
                 kinematically feasibly region in the lean/steer plane.

                """'''
    f_s += "    # Generated " + time.asctime() + "\n"
    f_s += "    cdef np.float64_t l = x[0]\n"
    f_s += "    cdef np.float64_t p = x[1]\n"
    f_s += "    cdef np.float64_t rr = params[0]\n"
    f_s += "    cdef np.float64_t rrt = params[1]\n"
    f_s += "    cdef np.float64_t rf = params[2]\n"
    f_s += "    cdef np.float64_t rft = params[3]\n"
    f_s += "    cdef np.float64_t lr = params[4]\n"
    f_s += "    cdef np.float64_t lf = params[5]\n"
    f_s += "    cdef np.float64_t ls = params[6]\n"
    f_s += "    cdef np.float64_t s = params[7]\n"
    f_s += "    cdef np.ndarray[DTYPE_t, ndim=1] F = np.zeros([2], dtype=DTYPE)\n"
    f_s += "    F[0] = " + f0_s + "\n"
    f_s += "    F[1] = " + f1_s + "\n"
    f_s += "    return F\n\n\n"


    f_s += "@cython.boundscheck(False)\n"
    f_s += "def df(np.ndarray[DTYPE_t, ndim=1] x, np.ndarray[DTYPE_t, ndim=1] params):\n"
    f_s += '''    """Evaluates holonomic constraint and its partial derivative
                 with respect to pitch, with the steer angle treated as parameter.

                 x:  Numpy array of lean and steer, length 2

                 params:  Numpy array of parameters, length 8
                          in the following order:
                         rr:  Rear wheel radius.
                        rrt:  Rear wheel tire radius.
                         rf:  Front wheel radius.
                        rft:  Front wheel tire radius.
                         lr:  Rear wheel center perpendicular distance from steer axis.
                         lf:  Front wheel center perpendicular distance from steer axis.
                         ls:  Steer axis offset.
                          s:  Steer angle
                 """'''
    f_s += "    # Generated " + time.asctime() + "\n"
    f_s += "    cdef np.float64_t l = x[0]\n"
    f_s += "    cdef np.float64_t p = x[1]\n"
    f_s += "    cdef np.float64_t rr = params[0]\n"
    f_s += "    cdef np.float64_t rrt = params[1]\n"
    f_s += "    cdef np.float64_t rf = params[2]\n"
    f_s += "    cdef np.float64_t rft = params[3]\n"
    f_s += "    cdef np.float64_t lr = params[4]\n"
    f_s += "    cdef np.float64_t lf = params[5]\n"
    f_s += "    cdef np.float64_t ls = params[6]\n"
    f_s += "    cdef np.float64_t s = params[7]\n"
    f_s += "    cdef np.ndarray[DTYPE_t, ndim=2] dF = np.zeros([2,2], dtype=DTYPE)\n"
    f_s += "    dF[0, 0] = " + J00_s + "\n"
    f_s += "    dF[0, 1] = " + J01_s + "\n"
    f_s += "    dF[1, 0] = " + J10_s + "\n"
    f_s += "    dF[1, 1] = " + J11_s + "\n"
    f_s += "    return dF\n"

    f_s += "@cython.boundscheck(False)\n"
    f_s += "def convert_geometric(np.ndarray[DTYPE_t, ndim=1] x):\n"
    f_s += "    cdef np.ndarray[DTYPE_t, ndim=1] lengths = np.zeros([3], dtype=DTYPE)\n"
    f_s += "    cdef np.float64_t w = x[0]\n"
    f_s += "    cdef np.float64_t c = x[1]\n"
    f_s += "    cdef np.float64_t lmbda = x[2]\n"
    f_s += "    cdef np.float64_t rr = x[3]\n"
    f_s += "    cdef np.float64_t rrt = x[4]\n"
    f_s += "    cdef np.float64_t rf = x[5]\n"
    f_s += "    cdef np.float64_t rft = x[6]\n"
    f_s += "    lengths[0] = (w+c)*np.cos(lmbda)-(rr+rrt)*np.sin(lmbda)\n"
    f_s += "    lengths[1]  = (rf+rft)*np.sin(lmbda)-c*np.cos(lmbda)\n"
    f_s += "    lengths[2] = w*np.sin(lmbda) + (rr+rrt-rf-rft)*np.cos(lmbda)\n"
    f_s += "    return lengths\n"

    fh.write(f_s)
    fh.close()
    cythonit('kinematic_feasibility')
    compileit('kinematic_feasibility.c', 'kinematic_feasibility.so')

# Uncomment if you want to regenerate the cython generated kinematic
# feasibility module
#generate_kf_module()

import kinematic_feasibility as kf
w = 1.02
c = 0.08
lmbda = np.pi/10.
rr = 0.3
rrt = 0.0
rf = 0.35
rft = 0.0
benchmark_parameters = np.array([w, c, lmbda, rr, rrt, rf, rft],\
        dtype=np.float64)

lr,lf,ls = kf.convert_geometric(benchmark_parameters)

x = np.array([0.0, lmbda], dtype=np.float64)
params = np.array([rr, rrt, rf, rft, lr, lf, ls, 0.0], dtype=np.float64)

step = 0.01
steer_range_12 = np.arange(0.01, np.pi, step)
steer_range_34 = np.arange(-0.01, -np.pi, -step)

# First quadrant
# Initial guess on lean and pitch
x0 = [np.pi/2. + 0.01, np.pi/2 - lmbda]
# List for storing the lean and pitch
lp_range1 = []
for steer in steer_range_12:
    params[-1] = steer
    lp_range1.append(fsolve(kf.f, x0, args=params))
    x0 = lp_range1[-1]
lp1 = np.array(lp_range1)

# Second quadrant
# Initial guess on lean and pitch
x0 = [-np.pi/2. + 0.01, np.pi/2 - lmbda]
# List for storing the lean and pitch
lp_range2 = []
for steer in steer_range_12:
    params[-1] = steer
    lp_range2.append(fsolve(kf.f, x0, args=params))
    x0 = lp_range2[-1]
lp2 = np.array(lp_range2)

# Third quadrant
# Initial guess on lean and pitch
x0 = [-np.pi/2. - 0.01, np.pi/2 - lmbda]
# List for storing the lean and pitch
lp_range3 = []
for steer in steer_range_34:
    params[-1] = steer
    lp_range3.append(fsolve(kf.f, x0, args=params))
    x0 = lp_range3[-1]
lp3 = np.array(lp_range3)

# Fourth quadrant
# Initial guess on lean and pitch
x0 = [np.pi/2. + 0.01, np.pi/2 - lmbda]
# List for storing the lean and pitch
lp_range4 = []
for steer in steer_range_34:
    params[-1] = steer
    lp_range4.append(fsolve(kf.f, x0, args=params))
    x0 = lp_range4[-1]
lp4 = np.array(lp_range4)


import matplotlib.pyplot as plt
plt.plot(lp1[:,0], steer_range_12)
plt.plot(lp2[:,0], steer_range_12)
plt.plot(lp3[:,0], steer_range_34)
plt.plot(lp4[:,0], steer_range_34)
plt.axis([-np.pi/2, np.pi/2., -np.pi, np.pi])
plt.xlabel('Lean')
plt.ylabel('Steer')
plt.show()


