from __future__ import division
import numpy as np
cimport numpy as np
DTYPE = np.float64
ctypedef np.float64_t DTYPE_t
cimport cython
@cython.boundscheck(False)
def f(np.ndarray[DTYPE_t, ndim=1] x, np.ndarray[DTYPE_t, ndim=1] params):
    """Computes the holonomic constraint and its partial derivative
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

                """    # Generated Thu Sep 10 13:14:27 2009
    cdef np.float64_t l = x[0]
    cdef np.float64_t p = x[1]
    cdef np.float64_t rr = params[0]
    cdef np.float64_t rrt = params[1]
    cdef np.float64_t rf = params[2]
    cdef np.float64_t rft = params[3]
    cdef np.float64_t lr = params[4]
    cdef np.float64_t lf = params[5]
    cdef np.float64_t ls = params[6]
    cdef np.float64_t s = params[7]
    cdef np.ndarray[DTYPE_t, ndim=1] F = np.zeros([2], dtype=DTYPE)
    F[0] = rft - rrt + lf*(np.sin(l)*np.sin(s) - np.cos(l)*np.cos(s)*np.sin(p)) - rr*np.cos(l) + ls*np.cos(l)*np.cos(p) - lr*np.cos(l)*np.sin(p) + rf*(np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))*(-np.cos(s)*np.sin(l) - np.cos(l)*np.sin(p)*np.sin(s))/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(1/2) + rf/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(1/2)
    F[1] = -lr*np.cos(l)*np.cos(p) - ls*np.cos(l)*np.sin(p) - lf*np.cos(l)*np.cos(p)*np.cos(s) + rf*(np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))*np.cos(l)*np.cos(p)*np.sin(s)/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(3/2) + rf*(-np.cos(s)*np.sin(l) - np.cos(l)*np.sin(p)*np.sin(s))*np.cos(l)*np.cos(p)*np.sin(s)/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(1/2) - rf*(np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))*np.cos(l)*np.cos(p)*np.sin(s)/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(1/2) + rf*(np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2*(-np.cos(s)*np.sin(l) - np.cos(l)*np.sin(p)*np.sin(s))*np.cos(l)*np.cos(p)*np.sin(s)/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(3/2)
    return F


@cython.boundscheck(False)
def df(np.ndarray[DTYPE_t, ndim=1] x, np.ndarray[DTYPE_t, ndim=1] params):
    """Evaluates holonomic constraint and its partial derivative
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
                 """    # Generated Thu Sep 10 13:14:27 2009
    cdef np.float64_t l = x[0]
    cdef np.float64_t p = x[1]
    cdef np.float64_t rr = params[0]
    cdef np.float64_t rrt = params[1]
    cdef np.float64_t rf = params[2]
    cdef np.float64_t rft = params[3]
    cdef np.float64_t lr = params[4]
    cdef np.float64_t lf = params[5]
    cdef np.float64_t ls = params[6]
    cdef np.float64_t s = params[7]
    cdef np.ndarray[DTYPE_t, ndim=2] dF = np.zeros([2,2], dtype=DTYPE)
    dF[0, 0] = lf*(np.cos(l)*np.sin(s) + np.cos(s)*np.sin(l)*np.sin(p)) + rr*np.sin(l) + lr*np.sin(l)*np.sin(p) - ls*np.cos(p)*np.sin(l) + rf*(np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))*(-np.cos(l)*np.cos(s) + np.sin(l)*np.sin(p)*np.sin(s))/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(1/2) + rf*(-np.cos(s)*np.sin(l) - np.cos(l)*np.sin(p)*np.sin(s))*(np.cos(l)*np.cos(s) - np.sin(l)*np.sin(p)*np.sin(s))/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(1/2) + rf*(np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))*(2*np.cos(l)*np.cos(s) - 2*np.sin(l)*np.sin(p)*np.sin(s))/(2*(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(3/2)) + rf*(np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2*(-np.cos(s)*np.sin(l) - np.cos(l)*np.sin(p)*np.sin(s))*(2*np.cos(l)*np.cos(s) - 2*np.sin(l)*np.sin(p)*np.sin(s))/(2*(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(3/2))
    dF[0, 1] = -lr*np.cos(l)*np.cos(p) - ls*np.cos(l)*np.sin(p) - lf*np.cos(l)*np.cos(p)*np.cos(s) + rf*(np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))*np.cos(l)*np.cos(p)*np.sin(s)/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(3/2) + rf*(-np.cos(s)*np.sin(l) - np.cos(l)*np.sin(p)*np.sin(s))*np.cos(l)*np.cos(p)*np.sin(s)/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(1/2) - rf*(np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))*np.cos(l)*np.cos(p)*np.sin(s)/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(1/2) + rf*(np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2*(-np.cos(s)*np.sin(l) - np.cos(l)*np.sin(p)*np.sin(s))*np.cos(l)*np.cos(p)*np.sin(s)/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(3/2)
    dF[1, 0] = lr*np.cos(p)*np.sin(l) + ls*np.sin(l)*np.sin(p) + lf*np.cos(p)*np.cos(s)*np.sin(l) + rf*(np.cos(l)*np.cos(s) - np.sin(l)*np.sin(p)*np.sin(s))*np.cos(l)*np.cos(p)*np.sin(s)/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(3/2) + rf*(np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))*np.cos(p)*np.sin(l)*np.sin(s)/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(1/2) + rf*(-np.cos(l)*np.cos(s) + np.sin(l)*np.sin(p)*np.sin(s))*np.cos(l)*np.cos(p)*np.sin(s)/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(1/2) - rf*(np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))*np.cos(p)*np.sin(l)*np.sin(s)/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(3/2) - rf*(-np.cos(s)*np.sin(l) - np.cos(l)*np.sin(p)*np.sin(s))*np.cos(p)*np.sin(l)*np.sin(s)/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(1/2) - rf*(np.cos(l)*np.cos(s) - np.sin(l)*np.sin(p)*np.sin(s))*np.cos(l)*np.cos(p)*np.sin(s)/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(1/2) + rf*(np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2*(-np.cos(l)*np.cos(s) + np.sin(l)*np.sin(p)*np.sin(s))*np.cos(l)*np.cos(p)*np.sin(s)/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(3/2) - rf*(np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2*(-np.cos(s)*np.sin(l) - np.cos(l)*np.sin(p)*np.sin(s))*np.cos(p)*np.sin(l)*np.sin(s)/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(3/2) - rf*(np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2*(2*np.cos(l)*np.cos(s) - 2*np.sin(l)*np.sin(p)*np.sin(s))*np.cos(l)*np.cos(p)*np.sin(s)/(2*(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(3/2)) + 3*rf*(np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2*(2*np.cos(l)*np.cos(s) - 2*np.sin(l)*np.sin(p)*np.sin(s))*np.cos(l)*np.cos(p)*np.sin(s)/(2*(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(5/2)) + 3*rf*(np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**3*(-np.cos(s)*np.sin(l) - np.cos(l)*np.sin(p)*np.sin(s))*(2*np.cos(l)*np.cos(s) - 2*np.sin(l)*np.sin(p)*np.sin(s))*np.cos(l)*np.cos(p)*np.sin(s)/(2*(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(5/2)) + 3*rf*(np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))*(-np.cos(s)*np.sin(l) - np.cos(l)*np.sin(p)*np.sin(s))*(2*np.cos(l)*np.cos(s) - 2*np.sin(l)*np.sin(p)*np.sin(s))*np.cos(l)*np.cos(p)*np.sin(s)/(2*(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(3/2))
    dF[1, 1] = lr*np.cos(l)*np.sin(p) - ls*np.cos(l)*np.cos(p) + lf*np.cos(l)*np.cos(s)*np.sin(p) + rf*(np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))*np.cos(l)*np.sin(p)*np.sin(s)/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(1/2) - rf*(np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))*np.cos(l)*np.sin(p)*np.sin(s)/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(3/2) - rf*(-np.cos(s)*np.sin(l) - np.cos(l)*np.sin(p)*np.sin(s))*np.cos(l)*np.sin(p)*np.sin(s)/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(1/2) - rf*(np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2*(-np.cos(s)*np.sin(l) - np.cos(l)*np.sin(p)*np.sin(s))*np.cos(l)*np.sin(p)*np.sin(s)/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(3/2) + 3*rf*(np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**3*np.cos(l)**2*np.cos(p)**2*np.sin(s)**2*(-np.cos(s)*np.sin(l) - np.cos(l)*np.sin(p)*np.sin(s))/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(5/2) + 3*rf*np.cos(l)**2*np.cos(p)**2*np.sin(s)**2*(np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))*(-np.cos(s)*np.sin(l) - np.cos(l)*np.sin(p)*np.sin(s))/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(3/2) + rf*np.cos(l)**2*np.cos(p)**2*np.sin(s)**2/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(3/2) - 2*rf*np.cos(l)**2*np.cos(p)**2*np.sin(s)**2/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(1/2) - 2*rf*(np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2*np.cos(l)**2*np.cos(p)**2*np.sin(s)**2/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(3/2) + 3*rf*(np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2*np.cos(l)**2*np.cos(p)**2*np.sin(s)**2/(1 - (np.cos(s)*np.sin(l) + np.cos(l)*np.sin(p)*np.sin(s))**2)**(5/2)
    return dF
@cython.boundscheck(False)
def convert_geometric(np.ndarray[DTYPE_t, ndim=1] x):
    cdef np.ndarray[DTYPE_t, ndim=1] lengths = np.zeros([3], dtype=DTYPE)
    cdef np.float64_t w = x[0]
    cdef np.float64_t c = x[1]
    cdef np.float64_t lmbda = x[2]
    cdef np.float64_t rr = x[3]
    cdef np.float64_t rrt = x[4]
    cdef np.float64_t rf = x[5]
    cdef np.float64_t rft = x[6]
    lengths[0] = (w+c)*np.cos(lmbda)-(rr+rrt)*np.sin(lmbda)
    lengths[1]  = (rf+rft)*np.sin(lmbda)-c*np.cos(lmbda)
    lengths[2] = w*np.sin(lmbda) + (rr+rrt-rf-rft)*np.cos(lmbda)
    return lengths
