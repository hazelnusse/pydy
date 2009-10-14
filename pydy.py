from sympy import (Symbol, symbols, Basic, Function, Mul, Pow, Matrix, sin,
        cos, tan, cot, S, eye, Add, trigsimp, expand, pretty, Eq, collect, sqrt,
        sympify, factor, zeros, simplify, solve_linear_system, ratsimp,
        powsimp, block_diag, Derivative)
from sympy.printing.pretty.pretty import PrettyPrinter
from sympy.printing.str import StrPrinter
import time

e1 = Matrix([1, 0, 0])
e2 = Matrix([0, 1, 0])
e3 = Matrix([0, 0, 1])
zero = Matrix([0, 0, 0])
t = Symbol("t")

Basic.__str__ = lambda self: PyDyStrPrinter().doprint(self)
Basic.__repr__ = lambda self: PyDyStrPrinter().doprint(self)
Matrix.__str__ = lambda self: PyDyStrPrinter().doprint(self)

class UnitVector(Basic):
    """A standard unit vector  with a symbolic and a numeric representation.
    """

    def __init__(self, frame, i=0): #=-1,num=None):
        self.frame = frame    # Parent reference frame
        self.common_frames = set([frame])
        self.i = i
        self.v = {}
        s = frame.name
        if i == 1:
            self.v['sym'] = Symbol(s.lower()+str(i))
            self.v['num'] = e1
        elif i == 2:
            self.v['sym'] = Symbol(s.lower()+str(i))
            self.v['num'] = e2
        elif i == 3:
            self.v['sym'] = Symbol(s.lower()+str(i))
            self.v['num'] = e3
        elif i == 0:
            self.v['sym'] = Symbol(s.lower()+str(0))
            self.v['num'] = zero

    def __str__(self):
        return PyDyStrPrinter().doprint(self)

    def __repr__(self):
        return PyDyStrPrinter().doprint(self)

    def __cmp__(self, other):
        if isinstance(other, UnitVector):
            if self.frame == other.frame:
                return cmp(self.i, other.i)
            else:
                return cmp(len(self.frame.ref_frame_list),
                       len(other.frame.ref_frame_list))
        else:
            raise NotImplementedError()

    def __eq__(self, other):
        if isinstance(other, UnitVector):
            if other.frame == self.frame:
                return (self.v['num'] == other.v['num'])
            else:
                other_selfframe = other.express(self.frame)
                if isinstance(other_selfframe, UnitVector):
                    return (self.v['num'] == other_selfframe.v['num'])
                else:
                    return False
        elif isinstance(other, Vector):
            other_selfframe = other.express(self.frame)
            if isinstance(other_selfframe, UnitVector):
                return (self.v['num'] == other_selfframe.v['num'])
            else:
                return False
        elif isinstance(other, (Add, Mul)):
            other_as_Vector = Vector(other)
            return (self == other_as_Vector)
        else:
            return False

    def __mul__(self, other):
        if isinstance(other, Dyad):
            return NotImplemented
        else:
            return Basic.__mul__(self, other)

    def __neg__(self):
        return Vector({self: -1})

    def express(self, frame):
        """Express a UnitVector in a different reference frame.
        """

        if self.frame == frame:
            return self
        else:
            frame_list = self.frame.get_frames_list(frame)
            if len(self.common_frames) > 1:
                # Means that self is fixed in more than one frame
                for cf in self.common_frames:
                    if cf == self.frame:
                        continue
                    elif cf in frame_list:
                        fl = cf.get_frames_list(frame)
                        if len(fl) < len(frame_list):
                            frame_list = fl
                matrices = frame_list[0].get_rot_matrices(frame_list[-1])
            else:
                matrices = self.frame.get_rot_matrices(frame)
            if len(matrices) == 1 and matrices[0]*self.v['num'] == self.v['num']:
                return frame[self.i]
            else:
                u = self.v['num']
                for m in reversed(matrices):
                    u = m*u
                return Vector(u[0]*frame[1] + u[1]*frame[2] + u[2]*frame[3])

    def dot(self, other):
        """UnitVector dot product.
        """
        if isinstance(other, UnitVector):
            nrf = self.frame.NewtonianReferenceFrame
            if (self, other) in nrf.uv_dot_products:
                return nrf.uv_dot_products[(self, other)]
            elif (other, self) in nrf.uv_dot_products:
                return nrf.uv_dot_products[(other, self)]
            else:
                c = other.express(self.frame)
                if isinstance(c, UnitVector):
                    dp = (self.v['num'].T * c.v['num'])[0]
                    self.frame.NewtonianReferenceFrame.uv_dot_products[(self, \
                        other)] = dp
                    return dp
                elif isinstance(c, Vector):
                    s = c.dict.get(self, 0)
                    self.frame.NewtonianReferenceFrame.uv_dot_products[(self, \
                        other)] = s
                    return s
                else:
                    raise NotImplementedError()
        elif isinstance(other, Vector):
            s = S(0)
            for k, c in other.dict.items():
                s += c*self.dot(k)
            return s
        elif isinstance(other, Dyad):
            return other.ldot(self)
        else:
            raise NotImplementedError()

    def cross(self, other):
        """UnitVector cross product.
        """
        def cross_with_Vector(self, c):            # local function
            cp = {}
            for k, coef in c.dict.items():
                term = self.cross(k)
                if isinstance(term, UnitVector):
                    cp[term] = cp.get(term, 0) + coef
                elif isinstance(term, Vector):
                    for kt in term.dict:
                        cp[kt] = cp.get(kt, 0) + coef*term.dict[kt]
            return Vector(cp)

        if isinstance(other, UnitVector):
            nrf = self.frame.NewtonianReferenceFrame
            if (self, other) in nrf.uv_cross_products:
                return nrf.uv_cross_products[(self, other)]
            elif (other, self) in nrf.uv_cross_products:
                return -nrf.uv_cross_products[(other, self)]
            else:
                c = other.express(self.frame)
                #print self.frame
                if isinstance(c, UnitVector):
                    cp_list = [self.v['num'][1]*c.v['num'][2] - \
                            self.v['num'][2]*c.v['num'][1], \
                            -self.v['num'][0]*c.v['num'][2] + \
                            self.v['num'][2]*c.v['num'][0],  \
                            self.v['num'][0]*c.v['num'][1] - \
                            self.v['num'][1]*c.v['num'][0]]
                    cp = {}
                    for (c, i) in zip(cp_list, [1, 2, 3]):
                        if c != 0:
                            cp.update({self.frame[i] : c})
                    if len(cp) == 0:
                        return Vector(0)
                    elif len(cp) == 1:
                        if cp.values()[0] == 1:
                            return cp.keys()[0]  # Return a UnitVector object
                        else:
                            return Vector(cp)
                    else:
                        cp1 = Vector(cp)
                        cp2 = cp1.express(other.frame)
                        if isinstance(cp2, UnitVector):
                            return cp2
                        else:
                            for k in cp1.dict:
                                cp1.dict[k] = trigsimp(cp1.dict[k])
                            for k in cp2.dict:
                                cp2.dict[k] = trigsimp(cp2.dict[k])
                            if isinstance(cp1, UnitVector):
                                return cp1
                            elif isinstance(cp2, UnitVector):
                                return cp2
                            elif len(cp1.dict) <= len(cp2.dict):
                                return cp1
                            else:
                                return cp2
                elif isinstance(c, Vector):
                    cp1 = cross_with_Vector(self, c)
                    cp2 = cp1.express(other.frame)
                    if isinstance(cp1, UnitVector):
                        return cp1
                    elif isinstance(cp2, UnitVector):
                        return cp2
                    else:
                        for k in cp1.dict:
                            cp1.dict[k] = trigsimp(cp1.dict[k])
                        for k in cp2.dict:
                            cp2.dict[k] = trigsimp(cp2.dict[k])
                        if isinstance(cp1, UnitVector):
                            return cp1
                        elif isinstance(cp2, UnitVector):
                            return cp2
                        elif len(cp1.dict) <= len(cp2.dict):
                            return cp1
                        else:
                            return cp2
        elif isinstance(other, Vector):
            return cross_with_Vector(self, other)
        else:
            raise NotImplementedError()

    def dt(self, diff_frame):
        """UnitVector time derivative.
        """
        if isinstance(diff_frame, ReferenceFrame):
            if self.frame == diff_frame:
                return Vector(0)
            else:
                return cross(self.frame.ang_vel(diff_frame), self)
        else:
            raise TypeError("Must provide a ReferenceFrame to take the \
                derivative in")

class Dyad(object):
    """General dyad expression.
    """
    def __init__(self, v):
        """ v should be an additive expression of the form:
        (sympy expression)*UnitVector*UnitVector

        The sympy expression is optional, but typically would be an inertia
        scalar.
        """
        self.dict = {}
        if isinstance(v, dict):
            self.dict = v
        elif v == 0 or v == {}:
            pass
        else:
            raise NotImplementedError()

    def __add__(self, other):
        if isinstance(other, Dyad):
            nd = {}
            for k, v in self.dict.items():
                nd[k] = v
            for k, v in other.dict.items():
                nd[k] = nd.get(k, 0) + v
            for k, v in nd.items():
                if v == 0:
                    nd.pop(k)
        return Dyad(nd)

    def __neg__(self):
        nd = {}
        for k, v in self.dict.items():
            nd[k] = -v
        return Dyad(nd)

    def __sub__(self, other):
        if isinstance(other, Dyad):
            return self.__add__(other.__neg__())

    def subs(self, subs_dict):
        """Substitute into the scalar coeffiecients of each dyadic term.
        """
        nd = {}
        for k, v in self.dict.items():
            nd[k] = v.subs(subs_dict)
        return Dyad(nd)

    def expand(self):
        """Expand scalar coefficients of each dyadic term.
        """
        nd = {}
        for k, v in self.dict.items():
            nd[k] = v.expand()
        return Dyad(nd)

    def n(self):
        """Numerically evaluate each scalar coefficient.
        """
        nd = {}
        for k, v in self.dict.items():
            nd[k] = v.n()
        return Dyad(nd)


    def ldot(self, other):
        """Multplitcation by a UnitVector/Vector on the left.
        v * Dyad

        Returns a UnitVector / Vector
        """
        vec_dict = {}
        if isinstance(other, (UnitVector, Vector)):
            for d_term, coeff in self.dict.items():
                scalar_part = coeff*other.dot(d_term.args[0])
                if d_term.is_Mul:
                    vec_dict[d_term.args[1]] = vec_dict.get(d_term.args[1], 0) \
                            + scalar_part
                elif d_term.is_Pow:
                    vec_dict[d_term.args[0]] = (vec_dict.get(d_term.args[0], 0)
                            + scalar_part)
                else:
                    raise NotImplementedError()
        return Vector(vec_dict)

    def rdot(self, other):
        """Multplitcation by a UnitVector/Vector on the right.
        Dyad * v

        Returns a UnitVector / Vector
        """
        vec_dict = {}
        if isinstance(other, (UnitVector, Vector)):
            for d_term, coeff in self.dict.items():
                if d_term.is_Mul:
                    #print d_term
                    scalar_part = coeff*other.dot(d_term.args[1])
                    vec_dict[d_term.args[0]] = (vec_dict.get(d_term.args[0], 0)
                            + scalar_part)
                elif d_term.is_Pow:
                    #print d_term.args
                    scalar_part = coeff*other.dot(d_term.args[0])
                    vec_dict[d_term.args[0]] = (vec_dict.get(d_term.args[0], 0)
                            + scalar_part)
                else:
                    raise NotImplementedError()

        return Vector(vec_dict)

    def __str__(self):
        return PyDyStrPrinter().doprint(self)

    def express(self, frame):
        """Express a Dyad with Unit Vectors fixed in a specified frame."""
        dyad_dict = {}
        for d_term, coeff in self.dict.items():
            if d_term.is_Mul:   # Case of A[1]*A[2]
                t1 = d_term.args[0].express(frame)
                t2 = d_term.args[1].express(frame)
                if isinstance(t1, UnitVector) and isinstance(t2, UnitVector):
                    dyad_dict[t1*t2] = dyad_dict.get(t1*t2, 0) + coeff
                elif isinstance(t1, UnitVector) and isinstance(t2, Vector):
                    for k, v in t2.dict.items():
                        dyad_dict[t1*k] = dyad_dict.get(t1*k, 0) + coeff*v
                elif isinstance(t1, Vector) and isinstance(t2, UnitVector):
                    for k, v in t1.dict.items():
                        dyad_dict[k*t2] = dyad_dict.get(k*t2, 0) + coeff*v
                elif isinstance(t1, Vector) and isinstance(t2, Vector):
                    for k1, v1 in t1.dict.items():
                        for k2, v2 in t2.dict.items():
                            dyad_dict[k1*k2] = dyad_dict.get(k1*k2, 0) + coeff*v1*v2
            if d_term.is_Pow:
                t1 = d_term.args[0].express(frame)
                if isinstance(t1, UnitVector):
                    dyad_dict[t1**2] = dyad_dict.get(t1**2, 0) + coeff
                elif isinstance(t1, Vector):
                    t2 = Vector(t1.dict)  # make a copy
                    for k1, v1 in t1.dict.items():
                        for k2, v2 in t2.dict.items():
                            dyad_dict[k1*k2] = dyad_dict.get(k1*k2, 0) + coeff*v1*v2
        return Dyad(dyad_dict)

class Inertia(Dyad):
    """Inertia dyadic.
    """
    def __init__(self, F, scalars):
        """Specify frame, scalars as:
        frame - ReferenceFrame
        scalars - List or tuple of I11, I22, I33, I12, I23, I13 inertia scalars
        """
        I11, I22, I33, I12, I23, I13 = scalars
        self.dict = {}
        for i, s in enumerate(scalars):
            if s == 0:
                continue
            if i <= 2:
                self.dict[F[i+1]*F[i+1]] = s
            elif i == 3:
                self.dict[F[1]*F[2]] = s
                self.dict[F[2]*F[1]] = s
            elif i == 4:
                self.dict[F[2]*F[3]] = s
                self.dict[F[3]*F[2]] = s
            elif i == 5:
                self.dict[F[1]*F[3]] = s
                self.dict[F[3]*F[1]] = s

class Vector(Basic):
    """Symbolic vector expression.

    Internally represented as a dictionary whose keys are UnitVectors and whose
    values are the corresponding coefficient of that UnitVector.

    Example

    ::

        >>> N = ReferenceFrame("N")
        >>> x, y, z = symbols('x y z')
        >>> v = Vector(x*N[1] + y*N[2] + z*N[3])
        >>> v.dict == {N[1]: x, N[2]: y, N[3]: z}
        True

    """

    def __init__(self, v):
        """Initialize a Vector object.

        Example


        Method 1:

        ::

            >>> N = ReferenceFrame("N")
            >>> x, y, z = symbols('x y z')
            >>> v = Vector(x*N[1] + y*N[2] + z*N[3])
            >>> v
            x*n1> + y*n2> + z*n3>


        Method 2:

        ::

            >>> v = Vector({N[1]: x, N[2]: y, N[3]: z})
            >>> v
            x*n1> + y*n2> + z*n3>

        See also
        """

        if isinstance(v, dict):
            for k in v.keys():
                v[k] = sympify(v[k])
                if v[k] == 0:  v.pop(k)
            self.dict = v
        elif isinstance(v, Vector):
            self.dict = v.dict
        else:
            vdict = self.parse_terms(v)
            for k in vdict.keys():
                vdict[k] = sympify(vdict[k])
                if vdict[k] == 0:  vdict.pop(k)
            self.dict = vdict

    def __str__(self):
        return PyDyStrPrinter().doprint(self)

    def __repr__(self):
        return PyDyStrPrinter().doprint(self)

    def __add__(self, other):
        """Add two Vector objects.

        Example


        ::

            >>> N = ReferenceFrame('N')
            >>> v1 = Vector(2*N[1])
            >>> v2 = Vector(N[1] + 3*N[2])
            >>> v3 = v1 + v2
            >>> v3
            3*n1> + 3*n2>

        See Also

        L{__sub__}
        """
        if isinstance(other, Vector):
            sum = dict([(k, self.dict.get(k, 0) + other.dict.get(k, 0)) for k in
                    (self.dict.keys() + other.dict.keys())])

            if len(sum) == 1 and sum[sum.keys()[0]] == 1:
                return sum.keys()[0]
            else:
                return Vector(sum)
        elif isinstance(other, UnitVector):
            return self.__add__(Vector({other: S(1)}))
        elif isinstance(other, (Add, Mul)):
            return self.__add__(Vector(other))
        else:
            raise NotImplementedError()

    def __sub__(self, other):
        """Subtract two Vector objects.

        Example

        ::

            >>> N = ReferenceFrame('N')
            >>> v1 = Vector(2*N[1])
            >>> v2 = Vector(N[1] + 3*N[2])
            >>> v3 = v1 - v2
            >>> v3
            n1> - 3*n2>

        See Also

        L{__add__}
        """
        if isinstance(other, Vector):
            dif = dict([(k, self.dict.get(k, 0) - other.dict.get(k, 0)) for k in
                    (self.dict.keys() + other.dict.keys())])
            if len(dif) == 1 and dif[dif.keys()[0]] == 1:
                return dif.keys()[0]
            else:
                return Vector(dif)
        elif isinstance(other, UnitVector):
            return self.__sub__(Vector({other: S(1)}))
        elif isinstance(other, (Add, Mul)):
            return self.__sub__(Vector(other))
        else:
            raise NotImplementedError()

    def __eq__(self, other):
        """Compares two Vector objects for equality.
        """
        if isinstance(other, Vector):
            if self.dict == other.dict:         # Easy case
                return True
            else:                               # More difficult case
                l1 = len(self.dict)
                k1 = self.dict.keys()
                l2 = len(other.dict)
                k2 = other.dict.keys()
                if l1 == 0 or l2 == 0:
                    return False
                else:
                    self_in_first_key_frame = self.express(k1[0].frame)
                    other_in_first_key_frame = other.express(k1[0].frame)
                    # Now the two vectors are both expressed in the same
                    # coordinate frame
                    if self_in_first_key_frame.dict == \
                            other_in_first_key_frame.dict:
                        return True
                    else:
                        return False
        elif isinstance(other, UnitVector):
            v = self.express(other.frame)
            if isinstance(v, UnitVector):
                if v == other:  return True
            else:
                return False
        else:
            other_as_Vector = Vector(other)
            return self == other_as_Vector

    """
    def __rmul__(self, other):
        if isinstance(other, Dyad):
            return NotImplemented
        elif isinstance(other, Basic) and not isinstance(other, UnitVector):
            product = {}
            for k in self.dict:
                product[k] = other*self.dict[k]
            return Vector(product)
        else:
            return NotImplemented

    def __mul__(self, other):
        if isinstance(other, Dyad):
            return NotImplemented
        else:
            return Basic.__mul__(self, other)

    def __lmul__(self, other):
        if isinstance(other, Dyad):
            return NotImplemented
        elif isinstance(other, Symbol):
            prod = {}
            for k in self.dict:
                prod[k] = other*self.dict[k]
        else:
            return Basic.__mul__(self, other)
    """

    def __neg__(self):
        return Vector(dict([(k, -self.dict[k]) for k in self.dict]))

    def coeffv(self, scalar):
        """Vector coefficient of a scalar
        """

        return Vector(dict([(k, self.dict[k].coeff(scalar)) for k in
            self.dict if self.dict[k].coeff(scalar) != None]))


    def cross(self, other):
        if isinstance(other, Vector):
            vcp = {}
            for k in self.dict:
                for ko in other.dict:
                    kcrossko = k.cross(ko)
                    if isinstance(kcrossko, UnitVector):
                        vcp[kcrossko] = (vcp.get(kcrossko, 0) +
                                self.dict[k]*other.dict[ko])
                    else:
                        for uv_term in kcrossko.dict:
                            vcp[uv_term] = (vcp.get(uv_term, 0) +
                                    self.dict[k]*other.dict[ko]*
                                    kcrossko.dict[uv_term])
            return Vector(vcp)
        elif isinstance(other, UnitVector):
            vcp = {}
            for k in self.dict:
                k_cross_other = k.cross(other)
                if isinstance(k_cross_other, UnitVector):
                    vcp[k_cross_other] = (vcp.get(k_cross_other, 0) +
                        self.dict[k])
                else:
                    for uv_term in k_cross_other.dict:
                        vcp[uv_term] = (vcp.get(uv_term, 0) +
                            self.dict[k]*k_cross_other.dict[uv_term])
            return Vector(vcp)
        elif isinstance(other, (Add, Mul)):
            return self.cross(Vector(other))
        else:
            raise NotImplementedError()

    def dot(self, other):
        """Vector dot product.
        """
        if isinstance(other, Vector):
            s = S(0)
            for k in self.dict:
                s += sum([self.dict[k]*other.dict[ko]*k.dot(ko) for ko in
                    other.dict])
            return s
        elif isinstance(other, UnitVector):
            s = sum([self.dict[k]*k.dot(other) for k in self.dict])
            return s
        elif isinstance(other, (Add, Mul)):
            return self.dot(Vector(other))
        elif isinstance(other, Dyad):
            return other.ldot(self)
        else:
            raise NotImplementedError()

    def dt(self, frame):
        if isinstance(frame, ReferenceFrame):
            dt_self = {}
            for k in self.dict:
                # First term comes from time differentiating in the frame of
                # the UnitVector frame of k
                dt_self[k] = dt_self.get(k, 0) + (self.dict[k]).diff(t)
                # Second term comes from the omega cross term
                t2 = k.frame.ang_vel(frame).cross(k)
                if isinstance(t2, UnitVector):
                    dt_self[k] = dt_self.get(k, 0) + self.dict[k]
                else:       # Must be a Vector
                    for term in t2.dict:
                        dt_self[term] = (dt_self.get(term, 0) +
                            self.dict[k]*t2.dict[term])
            if len(dt_self) == 1:
                if dt_self.values()[0] == 1:
                    return dt_self.keys()[0]        # Return a UnitVector
                else:
                    return Vector(dt_self)
            else:
                return Vector(dt_self)

    def express(self, frame):
        """Expresses a Vector with UnitVectors fixed in the specified frame.
        """

        new = {}

        for uv in self.dict.keys():
            # Convert each unit vector term to the desired frame
            uv_in_frame = uv.express(frame)

            # Case for UnitVectors
            if isinstance(uv_in_frame, UnitVector):
                new[uv_in_frame] = new.get(uv_in_frame, 0) + self.dict[uv]

            # Case for Vectors
            elif isinstance(uv_in_frame, Vector):
                # Go through each term
                for uv_term, coef in uv_in_frame.dict.items():
                    new[uv_term] = (new.get(uv_term, 0) +
                            self.dict[uv]*coef)

        for uv in new.keys():
            new[uv] = new[uv].expand().subs(uv.frame.NewtonianReferenceFrame.csqrd_dict).expand()
            #new[uv] = expand(trigsimp(new[uv]))
            #new[uv] = trigsimp(expand(trigsimp(new[uv])))
            if new[uv] == 0: new.pop(uv)

        if len(new) == 1 and new.values()[0] == 1:
            return new.keys()[0]
        else:
            return Vector(new)

    @property
    def mag(self):
        """Magnitude of a Vector.
        """
        return sqrt(self.mag_sqr)

    @property
    def mag_sqr(self):
        """Magnitude squared of a Vector.
        """
        m = 0
        s = set([])
        for k1, k2 in ((x,y) for x in self.dict for y in self.dict):
            if (k2, k1) in s:
                continue
            else:
                s.add((k1, k2))
        for k1, k2 in s:
            if k1 == k2:
                    #m += expand(self.dict[k1]**2)
                    m += self.dict[k1]**2
            else:
                    #m += 2*expand(self.dict[k1]*self.dict[k2]*dot(k1, k2))
                    m += 2*self.dict[k1]*self.dict[k2]*dot(k1, k2)

        # Try to factor things if possible
        """
        if isinstance(m, Add):
            trigterms = m.atoms(sin, cos)
            replacements = []
            for i, t in enumerate(trigterms):
                replacements.append(Symbol('Trig%d' % i, dummy=True))
            subsdict = dict(zip(trigterms, replacements))
            rsubsdict = dict(zip(replacements, trigterms))
            trigadd = []
            otheradd = []
            for arg in m.args:
                if arg.atoms(sin, cos):
                    trigadd.append(arg)
                else:
                    otheradd.append(arg)
            trigexpr = S(0)
            otherexpr = S(0)
            for term in trigadd: trigexpr += term
            for term in otheradd: otherexpr += term
            trigexprs = trigexpr.subs(subsdict)
            trigexprsf = factor(trigexprs)
            trigexprf = trigexprsf.subs(rsubsdict)
            m = trigexprf + otherexpr
            #m = trigexpr + otherexpr
        """
        return m

    @property
    def normalized(self):
        v = Vector(0)
        m = self.mag
        for k in self.dict:
            v.dict[k] = self.dict[k] / m
        return v

    def parse_terms(self, v):
        """
        Given a Sympy expression with UnitVector terms, return a dictionary
        whose keys are the UnitVectors and whose values are the coeefficients
        of the UnitVectors
        """

        if v == 0:
            return {}
        elif isinstance(v, UnitVector):
            return {v: S(1)}
        elif isinstance(v, Vector):
            return v.dict
        elif isinstance(v, Mul):
            v = v.expand()  #  Could expand out to an Add instance
            if isinstance(v, Add):
                return self.parse_terms(v)  # If this happens, reparse
            elif isinstance(v, Mul):   # Otherwise it is still a Mul instance
                """if v.atoms(Vector):
                    i = v.args.index(list(v.atoms(Vector))[0])
                    scalarpart = S(1)
                    for j, term in enumerate(v.args):
                        if j == i: continue
                        else:
                            scalarpart *= term
                    newvdict = {}
                    for k in v.args[i].dict:
                        newvdict[k] = scalarpart*v.args[i].dict[k]
                    return newvdict
                """
                args = v.args
                term_types_list = [type(terms) for terms in args]
                if UnitVector in term_types_list:
                    i = term_types_list.index(UnitVector)
                    coefs = args[:i] + args[i+1:]
                    prod_coefs = S(1)
                    for p in coefs:
                        prod_coefs *= p
                    return {args[i]: prod_coefs}
                elif Vector in term_types_list:
                    i = term_types_list.index(Vector)
                    coefs = args[:i] + args[i+1:]
                    prod_coefs = S(1)
                    for p in coefs:
                        prod_coefs *= p
                    for k in args[i].dict.keys():
                        args[i].dict[k] *= prod_coefs
                    return args[i].dict
                else:
                    raise NotImplementedError()
            else:
                raise NotImplementedError()
        elif isinstance(v, Pow):
            v = v.expand()
            #  I don't think this will ever be entered into.
            #  You would have to have something like A[1]*A[2],
            #  which isn't a valid vector expression.
            #  Or q1*q2, which is a scalar expression.
            for b in v.args:
                if isinstance(b, UnitVector):
                    return {b: v.coeff(b)}
        elif isinstance(v, Add):
            v = v.expand()
            terms = {}
            for add_term in v.args:
                if isinstance(add_term, (Mul, Pow)):
                    add_term_dict = self.parse_terms(add_term)
                elif isinstance(add_term, UnitVector):
                    add_term_dict = {add_term: S(1)}
                elif isinstance(add_term, Vector):
                    add_term_dict = add_term.dict
                else:
                    raise NotImplementedError()
                for k in add_term_dict:
                    terms[k] = terms.get(k, 0) + add_term_dict[k]
            return terms
        else:
            return NotImplemented

    def partials(self, u_list):
        """Computes partial velocities.
        """
        return [self.coeffv(u) for u in u_list]

    def subs(self, subs_dict):
        return Vector(dict([(k, self.dict[k].subs(subs_dict)) for k in
            self.dict]))

    def expandv(self):
        """Expands each coefficient of a Vector's UnitVectors
        """
        ex = {}
        for uv, c in self.dict.items():
            ex[uv] = c.expand()
        return Vector(ex)

class Point(object):
    """
    A class for keeping track of the relative position, velocity, and
    acceleration of points.  Points can be created in three ways:

    Method 1:
    P = Point('P')

    Method 2:
    Q = P.locate('Q', r)

    Method 3:
    Q = P.locate('Q', r, frame)

    Methods 2 and 3 automatically form the velocity and acceleration of the new
    point Q relative to the parent point P.

    Method 1 is used to create the 'base' point from which all other points are
    derived.  Method 1 is automatically called when a 'base' reference frame is
    created, this point corresponds to what is commonly known as the inertial
    origin, and this syntax is typically never used explicitly.  Method 2
    creates a new point Q, located relative to point P by the Vector r.  Method
    3 creates a new point Q, located relative to point P by the Vector r, but
    it is assumed that this point is fixed in frame, so that the velocity of Q
    relative to P is the velocity of P plus the cross product of the angular
    velocity of frame relative to the Inertial Frame with r.
    """

    def __init__(self, name, relativeposition=None, parentpoint=None,
            fixedinframe=None, mass=None, force=None):
        # When instantiated by ReferenceFrame
        if not any([relativeposition, parentpoint]):
            self.name = name
            self.point_list = [self]
            self.pos = {self: Vector(0)}
            self._vrel = Vector(0)
            self.NewtonianFrame = fixedinframe
            self._fixedin = set([fixedinframe])
            self.parentpoint = None
            self.children = []
            self.mass = 0
            self.force = Vector(0)
        # When instantiated by locate method
        elif all([name, relativeposition, parentpoint]):
            relativeposition = Vector(relativeposition)
            self.name = name
            self.parentpoint = parentpoint
            self.children = []
            parentpoint.children.append(self)
            # Initialize the inertial velocity, relative to the parent point
            self._vrel = {}
            # Assign the vector pointing back to the parent to the new point's
            # position
            self.pos = {parentpoint: -relativeposition}
            # Add the position vector pointing from the parent to the new point
            # to the parent's pos dictionary
            parentpoint.pos[self] = relativeposition
            # Add self to the begin of the parent's point list
            self.point_list = [self] + parentpoint.point_list
            # Copy the Newtonian frame from the parent point.
            self.NewtonianFrame = parentpoint.NewtonianFrame
            # Set the mass and force
            self.mass = 0 if mass == None else mass
            self.force = Vector(0) if force == None else force
            # If an optional frame argument is specified, append it to the
            # _fixedin list of the parent and create the list for the new point
            if isinstance(fixedinframe, ReferenceFrame):
                parentpoint._fixedin.add(fixedinframe)
                self._fixedin = set([fixedinframe])
                self._vrel = cross(fixedinframe.ang_vel(self.NewtonianFrame),
                        relativeposition)
            elif fixedinframe == None:
                self._fixedin = set([])
                self._vrel = relativeposition.dt(self.NewtonianFrame)
            else:
                raise TypeError('fixedinframe must be a ReferenceFrame type')
        else:
            raise NotImplementedError()

    def apply_force(self, force, other=None, reset=False):
        """Apply force to a point or particle.

        Can be used to apply a force to a point or particle.

        Repeated calls to apply_force are additive.

        If optional argument other is specified, a torque or equal magnitude
        but opposite sign will applied to the other point.

        If you need to reset the applied force to zero, use optional parameter
        reset=True
        """
        if other==None:
            if reset==False:
                self.force += Vector(force)
            elif reset==True:
                self.force = Vector(force)
            else:
                raise TypeError('reset must be a boolean')
        elif isinstance(other, ReferenceFrame):
            if reset==False:
                force = Vector(force)
                self.force += force
                other.force -= force
            elif reset==True:
                force = Vector(force)
                self.force = force
                other.force = -force
            else:
                raise TypeError('reset must be a boolean')
        else:
            raise TypeError('other must be a ReferenceFrame')

    def locate(self, name, r, frame=None, mass=None, force=None):
        """Returns a new Point located relative to the parent point.

        Introduces the concept of a point fixed in a frame.

        Method 1:

            P1 = N.O.locate('P1', r)

        Method 2:

            P2 = N.O.locate('P1', r, frame)

        Both methods assign the relative position vector r in an identical way,
        they differ in how the velocity of the point is determined.  Method 1
        takes the time derivative of the supplied vector r in the
        NewtonianFrame.  Method 2 treats the new point as fixed in the supplied
        frame, so the velocity is calculated as the cross product of the
        angular velocity of frame relative to the NewtonianFrame with the
        supplied position vector r.
        """
        r = Vector(r)
        newpoint = Point(name, relativeposition=r, \
                parentpoint=self, fixedinframe=frame)
        newpoint.mass = mass if mass is not None else 0
        newpoint.force = Vector(force) if force is not None else Vector(0)
        return newpoint

    def rel(self, other):
        """
        Returns the position from Point other to Point self, i.e., the position
        of self relative to other.
        """
        if isinstance(other, Point):
            pl = self.get_point_list(other)
            pos = Vector(0)
            for i, p in enumerate(pl[:-1]):
                pos -= pl[i].pos[pl[i+1]]
            return pos
        elif other is None:
            return Vector(0)

    def vel(self, point=None, frame=None):
        """Calculate the velocity of a point.

        Used without arguments, vel() returns the velocity of self
        relative to the Newtonian Frame, taking into account whether points
        have been declared as fixed in a frame or not.

        Used with arguments, .vel(point, frame) returns the velocity of self
        relative to point, as view by an observer fixed in frame.
        """

        v = Vector(0)
        if point == frame == None:
            if hasattr(self, 'abs_vel'):
                return self.abs_vel
            else:
                for p in self.point_list:
                    v += p._vrel
        elif isinstance(point, Point) and isinstance(frame, ReferenceFrame):
            # Get the point list from point to self
            point_list = point.get_point_list(self)
            for i, pa in enumerate(point_list[:-1]):
                pb = point_list[i+1]
                set_intersect = pa._fixedin & pb._fixedin
                # Case when the two points are not fixed in the same frame
                if len(set_intersect) == 0:
                    v += dt(pb.rel(pa), frame)
                # Case when the two points are fixed in the same frame
                elif len(set_intersect) == 1:
                    v += cross(set_intersect.pop().ang_vel(frame),
                            pb.rel(pa))
                else:
                    raise NotImplementedError('Somehow these two points are \
                        both fixed in 2 or more of the same frames')
        return v

    def get_point_list(self, other=None):
        """
        Gets the list of Points between Point self and Point other, including both.
        """
        if other == None:
            return self.point_list
        elif self == other:
            return [self]
        else:
            r2t = list(reversed(other.point_list))

            if len(self.point_list) == 1:
                return r2t
            elif len(r2t) == 1:
                return self.point_list

            r1t = list(reversed(self.point_list))
            i = 1

            while r1t[i] == r2t[i]:
                del r1t[i-1]
                del r2t[i-1]
                if len(r1t)<2 or len(r2t)<2:
                    break

            r1t.reverse()
            return r1t[:-1] + r2t

    def __str__(self):
        return '<Point %s>' % self.name

    def __repr__(self):
        return '<Point %s>' % self.name

class ReferenceFrame(object):
    """
    A standard reference frame with 3 mutually perpendicular unit vectors.
    Reference frames can be created in two ways:

    Method 1:
    A = ReferenceFrame('A')

    Method 2:
    B = A.rotate('B', axis, angle)
    where:
    axis = 1, 2 or 3
    angle is the radian measure of the rotation.

    Method 1 typically is used to create the 'base' frame from which all other
    frames are derived.  Method 2 is used to create all subsequent frames.  In
    doing so, circular definitions of frames are avoided and a tree structure
    of reference frames is created.  The first argument is a string and
    determines how the basis UnitVectors are printed.
    """

    def __init__(self, s, matrix=None, frame=None, omega=None):
        """
        If instantiated without the optional arguments, the 'base'
        ReferenceFrame is created.  The optional arguments are automatically
        generated by the rotate() method for the purpose of creating a new
        ReferenceFrame object.  See rotate() method for details of the optional
        arguments.
        """
        self.children = []
        self.name = s
        self.triad = [UnitVector(self, i) for i in (1,2,3)]
        self.transforms = {}
        self.parentframe = frame
        self.torque = Vector(0)

        if not any([matrix, frame, omega]):
            self.ref_frame_list = [self]
            self.O = Point(s + 'O', fixedinframe=self)
            self.point_list = [self.O]
            self.NewtonianReferenceFrame = self
            self.uv_dot_products = {}
            self.uv_cross_products = {}
            self.inertia = Inertia(self, (0,0,0,0,0,0))
        else:
            self.ref_frame_list = [self] + frame.ref_frame_list[:]
            self.NewtonianReferenceFrame = frame.NewtonianReferenceFrame

        if isinstance(omega, Vector):
            self._wrel = omega
            frame._wrel_children[self] = -omega
        elif isinstance(omega, tuple):
            # Case of simple rotations about 1 axis
            self._wrel = Vector(omega[1] * self.triad[omega[0]-1])
            frame._wrel_children[self] = Vector(-omega[1] *
                    self.triad[omega[0]-1])
        else:
            self._wrel = Vector(0)

        self._wrel_children = {}


        if frame is not None:
            frame.children.append(self)
            self.append_transform(frame, matrix)
            frame.append_transform(self, matrix.T)

    def __getitem__(self, i):
        """
        Reference the UnitVectors that are attached to a reference frame.
        Example:
        A = ReferenceFrame('A')

        A[1], A[2], A[3] are the three basis vectors of the A frame.
        """
        if i == 1 or i == 2 or i ==3:
            return self.triad[i-1]
        else:
            raise NotImplementedError()

    def append_transform(self, frame, matrix):
        """
        Appends 'matrix' to the transforms dict which transform vectors
        expressed in self basis vectors to the basis vectors of the frame
        'frame'
        """
        # We just append it to our "transforms" dict.
        self.transforms[frame] = matrix

    def rotate(self, name, axis, angle, I=None, I_frame=None):
        """Returns a new rotated reference frame.

        Perform simple rotations, Euler rotations, space fixed rotations, axis
        angle rotations, Euler parameter rotations, or Rodrigues parameter
        rotations.

        Perform a simple rotation about the 1, 2, or 3 axis, by an amount
        specified by angle, to create a new ReferenceFrame object.
        Automatically generates the angular velocity of the new frame with
        respect to the parent frame.

        Currently the orientation is stored as the direction cosine matrix,
        further work should implement quaternions.  Extra functionality such as
        the ability to specify a set of Euler angles or an arbitrary axis and
        angle needs to be implemented.

        When rotate is used to generate a new ReferenceFrame, the orientation
        of that reference frame relative to its parent reference frame is
        stored in both frames in the form of the direction cosine matrix.
        Additionally, the angular velocity of the new reference frame relative
        to the parent frame (and vice versa) is stored with the new reference
        frame (and the parent reference frame).
        """

        if not isinstance(angle, (list, tuple)):
            if axis in set((1, 2, 3)):
                matrix = self._rot(axis, angle)
                omega = (axis, angle.diff(t))
            elif axis in set((-1, -2, -3)):
                matrix = self._rot(-axis, -angle)
                omega = (-axis, -angle.diff(t))
            elif isinstance(axis, (UnitVector, Vector)):
                raise NotImplementedError("Axis angle rotations not \
                    implemented.")
            else:
                raise ValueError("Invalid axis")
            newframe = ReferenceFrame(name, matrix, self, omega)
            self[abs(axis)].common_frames.add(newframe)
            newframe[abs(axis)].common_frames.add(self)
            if I == None and I_frame==None:
                newframe.inertia = Inertia(newframe, (0, 0, 0, 0, 0, 0))
            elif I != None and I_frame != None:
                newframe.inertia = Inertia(I_frame, I)
            elif I != None and I_frame == None:
                newframe.inertia = Inertia(newframe, I)
            else:
                raise NotImplementedError()
            return newframe
        else:
            if len(angle) == 3:
                csqrd_dict = {cos(angle[0])**2:1-sin(angle[0])**2,\
                              cos(angle[1])**2:1-sin(angle[1])**2,\
                              cos(angle[2])**2:1-sin(angle[2])**2}
                rot_type = str(axis)
                if rot_type in set(('BODY123', 'BODY132', 'BODY231', 'BODY213',
                    'BODY312', 'BODY321', 'BODY121', 'BODY131', 'BODY232',
                    'BODY212', 'BODY313', 'BODY323', 'SPACE123', 'SPACE132',
                    'SPACE231', 'SPACE213', 'SPACE312', 'SPACE321', 'SPACE121',
                    'SPACE131', 'SPACE232', 'SPACE212', 'SPACE313', 'SPACE323')):
                    if rot_type[0] == 'B':  # Body fixed (Euler) angles
                        a1 = int(rot_type[4])
                        a2 = int(rot_type[5])
                        a3 = int(rot_type[6])
                        C = self._rot(a1, angle[0]) * self._rot(a2, angle[1]) * \
                            self._rot(a3, angle[2])
                    else:                   # Space fixed angles
                        a1 = int(rot_type[5])
                        a2 = int(rot_type[6])
                        a3 = int(rot_type[7])
                        C = self._rot(a3, angle[2]) * self._rot(a2, angle[1]) * \
                            self._rot(a1, angle[0])
                    # From Spacecraft Dynamics, by Kane, Likins, Levinson
                    # Eqns 1.10.(5-7), pg. 47
                    # Angular velocity components in new frame's basis vectors

                    w1 = (C[0,2]*(C[0,1].diff(t)) + C[1,2]*(C[1,1].diff(t)) + \
                        C[2,2]*(C[2,1].diff(t))).expand().subs(csqrd_dict).expand()

                    w2 = (C[1,0]*(C[1,2].diff(t)) + C[2,0]*(C[2,2].diff(t)) + \
                        C[0,0]*(C[0,2].diff(t))).expand().subs(csqrd_dict).expand()

                    w3 = (C[2,1]*(C[2,0].diff(t)) + C[0,1]*(C[0,0].diff(t)) + \
                        C[1,1]*(C[1,0].diff(t))).expand().subs(csqrd_dict).expand()


                    # First initialize with zero angular velocity
                    newFrame = ReferenceFrame(name, C, self, Vector({}))
                    # Angular velocity vector of newFrame relative to self
                    omega  = Vector({newFrame[1]: w1, newFrame[2]: w2,
                        newFrame[3]: w3})
                    newFrame.set_omega(omega, self, force=True)
                    if I == None and I_frame==None:
                        newFrame.inertia = Inertia(newFrame, (0, 0, 0, 0, 0, 0))
                    elif I != None and I_frame != None:
                        newFrame.inertia = Inertia(I_frame, I)
                    elif I != None and I_frame == None:
                        newFrame.inertia = Inertia(newFrame, I)
                    else:
                        raise NotImplementedError()
                    return newFrame
            else:
                raise NotImplementedError("angle must be a list/tuple of \
                        length 3")

    def _rot(self, axis, angle):
        """Returns direction cosine matrix for simple 1,2,3 rotations

        """
        if axis == 1:
            return Matrix([[1, 0, 0],
                [0, cos(angle), -sin(angle)],
                [0, sin(angle), cos(angle)]])
        elif axis == 2:
            return Matrix([[cos(angle), 0, sin(angle)],
                [0, 1, 0],
                [-sin(angle), 0, cos(angle)]])
        elif axis == 3:
            return Matrix([[cos(angle), -sin(angle), 0],
                [sin(angle), cos(angle), 0],
                [0, 0, 1]])

    def __str__(self):
        return '<Frame %s>' % self.name

    def __repr__(self):
        return "<Frame %s>" % self.name

    def get_frames_list(self, frame):
        """
        Returns a list of frames from "self" to "frame", including both.

        Example::

            N - A - D - E - F
                 \
                  B - C

        Then:

        C.get_frames_list(F) == [C, B, A, D, E, F]
        F.get_frames_list(C) == [F, E, D, A, B, C]
        """
        if self == frame:
            return [self]
        else:
            r2t = list(reversed(frame.ref_frame_list))

            if len(self.ref_frame_list) == 1:
                return r2t
            elif len(r2t) == 1:
                return self.ref_frame_list

            r1t = list(reversed(self.ref_frame_list))
            i = 1
            while r1t[i] == r2t[i]:
                del r1t[i-1]
                del r2t[i-1]
                if len(r1t)<2 or len(r2t)<2:
                    break

            r1t.reverse()
            return r1t[:-1] + r2t

    def apply_torque(self, torque, other=None, reset=False):
        """Apply torque to a reference frame or rigid body.

        Can be used to apply a torque to a reference frame or rigid body.

        Repeated calls to apply_torqe are additive.

        If optional argument other is specified, a torque or equal magnitude
        but opposite sign will applied to the other reference frame.

        If you need to reset the applied torque to zero, use optional parameter
        reset=True
        """
        if other==None:
            if reset==False:
                self.torque += Vector(torque)
            elif reset==True:
                self.torque = Vector(torque)
            else:
                raise TypeError('reset must be a boolean')
        elif isinstance(other, ReferenceFrame):
            if reset==False:
                torque = Vector(torque)
                self.torque += torque
                other.torque -= torque
            elif reset==True:
                torque = Vector(torque)
                self.torque = torque
                other.torque = -torque
            else:
                raise TypeError('reset must be a boolean')
        else:
            raise TypeError('other must be a ReferenceFrame')

    def get_rot_matrices(self, frame):
        """
        Returns a list of matrices to get from self to frame.
        # local function
        def shrink(rm_list):
            new_list = []
            for i, rm in enumerate(rm_list[:-1]):
                rmn = rm_list[i+1]
                # All true if simple rotation about 1 axis
                cmp1 = [rm[0,0]==rmn[0,0], rm[0,1]==rmn[0,1], rm[0,2]==rmn[0,2],
                        rm[1,0]==rmn[1,0], rm[2,0]==rmn[2,0]]
                # All true if simple rotation about 2 axis
                cmp2 = [rm[0,1]==rmn[0,1], rm[1,0]==rmn[1,0], rm[1,1]==rmn[1,1],
                        rm[1,2]==rmn[1,2], rm[2,1]==rmn[2,1]]
                # All true if simple rotation about 3 axis
                cmp3 = [rm[0,2]==rmn[0,2], rm[1,2]==rmn[1,2], rm[2,0]==rmn[2,0],
                        rm[2,1]==rmn[2,1], rm[2,2]==rmn[2,2]]
                if all(cmp1):
                    # create the matrix
                    break
                elif all(cmp2):
                    # create the matrix
                    break
                elif all(cmp3):
                    # create the matrix
                    break
                new_list.append(rm)
        """

        frames = self.get_frames_list(frame)
        if frames == [self]:
            return [eye(3)]
        result = []
        for i, f in enumerate(frames[:-1]):
            result.append(f.transforms[frames[i+1]])
        result.reverse()
        return result

    def set_omega(self, omega, frame, force=False):
        """Sets the angular velocity relative to another frame.
        """
        if self._wrel == Vector(0) or force:
            self._wrel = omega
        #if self.W == {} or force:
        #    self.W[frame] = omega
        else:
            raise ValueError("set_omega has already been called.")

    def ang_vel(self, frame=None):
        """Angular velocity relative to another frame.
        """

        if frame == self:
            return Vector(0)
        else:
            if frame == None:
                if hasattr(self, 'abs_ang_vel'):
                    return self.abs_ang_vel
                else:
                    frame = self.NewtonianReferenceFrame
            elif frame == self.NewtonianReferenceFrame and hasattr(self,
                'abs_ang_vel'):
                return self.abs_ang_vel

            om = Vector(0)
            fl = frame.get_frames_list(self)
            n = len(fl)
            for i, f in enumerate(fl[:-1]):
                if f == fl[i+1].parentframe:
                    om += fl[i+1]._wrel
                else:
                    om -= fl[i]._wrel
            return om

    def ang_acc(self, frame=None):
        """Angular acceleration relative to another frame.
        """

        if frame == self:
            return Vector(0)
        else:
            if frame == None: frame = self.NewtonianReferenceFrame
            if hasattr(self, 'abs_ang_acc'):
                return self.abs_ang_acc
            else:
                alpha = Vector(0)
                fl = frame.get_frames_list(self)
                n = len(fl)
                for i, f in enumerate(fl[:-1]):
                    if f == fl[i+1].parentframe:
                        alpha += fl[i+1]._alpharel
                    else:
                        alpha -= fl[i]._alpharel
            return alpha

    def get_omega_list(self, frame):
        """
        Returns a list of simple angular velocities from self to frame.
        """
        frames = self.get_frames_list(frame)
        if frames == [self]:
            return [Vector({})]
        result = []
        for i, f in enumerate(frames[:-1]):
            #result.append(f.W[frames[i+1]])
            result.append(f._wrel)
        return result

class NewtonianReferenceFrame(ReferenceFrame):
    """A Newtonian Reference Frame class.

    Includes a dextral set of UnitVectors and an origin:

    >>> from pydy import *
    >>> N = NewtonianReferenceFrame('N')
    >>> N[1]
    n1>
    >>> N[2]
    n2>
    >>> N[3]
    n3>
    >>> N.O
    <Point NO>
    >>>
    """
    def __init__(self, s):
        ReferenceFrame.__init__(self, s)
        # Holonomic constraint equations
        self.hc_eqns = []
        # Differentiatated holonomic constraint equations
        self.dhc_eqns = []
        # Nonholonomic constraint equations
        self.nhc_eqns = []
        self.symbol_dict = {}
        self.symbol_dict_back = {}
        self.trig_func_set = set([])
        self.cos_func_set = set([])
        self.sin_func_set = set([])
        self.tan_func_set = set([])
        self.csqrd_dict = {}
        self.crossterms = set([])

    def setkindiffs(self, eqn_list):#, dependent_speeds=None, acc=True):
        """Set the kinematic differential equations of the system.

        Must be given as a list of Relational objects.

        """
        self.kindiffs = eqn_list


    def setdyndiffs(self, eqns):
        """
        Sets the dynamic equations of motion.
        """
        self.dyndiffs = eqns

    def recursive_subs(self, PorF, expr_dict):
        # Substitute into appropriate velocity/angular velocity
        if isinstance(PorF, Point):
            PorF._vrel = PorF._vrel.subs(expr_dict)
        elif isinstance(PorF, ReferenceFrame):
            PorF._wrel = PorF._wrel.subs(expr_dict)
        else:
            raise NotImplementedError()

        #  Initiate recursion
        if PorF.children == []:
            return
        else:
            for child in PorF.children:
                self.recursive_subs(child, expr_dict)

    def recursive_partials(self, PorF):
        """Recursively form the relative partial velocities of each point and
        partial angular velocity of each reference frame.
        """
        # Substitute into appropriate velocity/angular velocity
        if isinstance(PorF, (Point, ReferenceFrame)):
            if isinstance(PorF, Point): pv = PorF.vel().partials(self.u_list)
            if isinstance(PorF, ReferenceFrame): pv = PorF.ang_vel().partials(self.u_list)
            if hasattr(self, 'u_dependent') and len(self.u_dependent) != 0:
                pv_i = [pv[i] for i in self.independent_ci]
                pv_d = Matrix([pv[i] for i in self.dependent_ci]).T
                #  Maybe should use a dummy matrix instead of
                #  u_dependent_transform... then back substitute once final
                #  equations are derived.
                con = matrixv_multiply(pv_d, self.T_con).tolist()[0]
                PorF.partialv = [pv_i[i] + con[i] for i in range(len(con))]
            else:
                # Case when the system has no constraints
                PorF.partialv = pv
        else:
            raise NotImplementedError

        #  Initiate recursion
        if PorF.children == []:
            return
        else:
            for child in PorF.children:
                self.recursive_partials(child)

    def declare_coords(self, string, number, list=True):
        """Declare the generalized coordinates and their time derivatives.
        """
        q_list, qdot_list = gcs(string, number, list)
        self.q_list = q_list
        self.qdot_list = qdot_list
        # Generate lists of Symbol objects instead of Function objects
        self.csqrd_dict = {}
        self.tan_dict = {}
        self.cot_dict = {}
        for q in q_list:
            self.csqrd_dict[cos(q)**2] = 1 - sin(q)**2
            self.tan_dict[sin(q)/cos(q)] = tan(q)
            self.cot_dict[cos(q)/sin(q)] = cot(q)
        self.q_list_s = [Symbol(str(q.func)) for q in q_list]
        sin_q_list = [sin(qs) for qs in self.q_list]
        cos_q_list = [cos(qs) for qs in self.q_list]
        tan_q_list = [tan(qs) for qs in self.q_list]
        trig_q_list = sin_q_list + cos_q_list + tan_q_list
        sin_q_list_s = [Symbol('s'+str(qs)) for qs in self.q_list_s]
        cos_q_list_s = [Symbol('c'+str(qs)) for qs in self.q_list_s]
        tan_q_list_s = [Symbol('t'+str(qs)) for qs in self.q_list_s]
        trig_q_list_s = sin_q_list_s + cos_q_list_s + tan_q_list_s
        self.qdot_list_s = [Symbol(str(q.func)+'d') for q in q_list]
        self.q_list_dict = dict(zip(q_list, self.q_list_s))
        self.q_list_dict_back = dict(zip(self.q_list_s, q_list))
        trig_q_dict = dict(zip(trig_q_list, trig_q_list_s))
        trig_q_dict_back = dict(zip(trig_q_list_s, trig_q_list))

        self.qdot_list_dict = dict(zip(qdot_list, self.qdot_list_s))
        self.qdot_list_dict_back = dict(zip(self.qdot_list_s, qdot_list))
        # Update the comprehensive symbol dictionaries
        for d in (self.q_list_dict, self.qdot_list_dict):
            self.symbol_dict.update(d)
        for d in (self.q_list_dict_back, self.qdot_list_dict_back):
            self.symbol_dict_back.update(d)
        self.trig_subs_dict = trig_q_dict
        self.trig_subs_dict_back = trig_q_dict_back
        return q_list, qdot_list

    def declare_speeds(self, string, number, lst=True):
        """Declare the generalized speeds and their time derivatives.
        """
        u_list, udot_list = gcs(string, number, lst)
        self.u_list = u_list
        self.udot_list = udot_list
        self.u_independent = u_list
        self.u_dependent = []
        self.udot_independent = udot_list
        self.udot_dedependent = []
        self.independent_ci = list(range(len(u_list)))
        self.dependent_ci = []

        # Generate lists of Symbol objects instead of Function objects
        self.u_list_s = [Symbol(str(u.func)) for u in u_list]
        self.udot_list_s = [Symbol(str(u.func)+'p') for u in u_list]

        # Generate a set of cross terms
        for ui in u_list:
            for uj in u_list:
                self.crossterms.update(set([ui*uj]))
            for qd in self.qdot_list:
                self.crossterms.update(set([ui*qd]))
        for qd_i in self.qdot_list:
            for qd_j in self.qdot_list:
                self.crossterms.update(set([qd_i*qd_j]))
        self.crossterms = list(self.crossterms)
        self.crossterms.sort()
        # Generate substitution dictionaries between Symbol and Function
        # representation of the coordinates, generalized speeds, and their
        # respective time derivatives
        self.u_list_dict = dict(zip(u_list, self.u_list_s))
        self.u_list_dict_back = dict(zip(self.u_list_s, u_list))
        self.udot_list_dict = dict(zip(udot_list, self.udot_list_s))
        self.udot_list_dict_back = dict(zip(self.udot_list_s, udot_list))

        for d in (self.u_list_dict, self.udot_list_dict):
            self.symbol_dict.update(d)
        for d in (self.u_list_dict_back, self.udot_list_dict_back):
            self.symbol_dict_back.update(d)

        return u_list, udot_list

    def declare_parameters(self, string):
        """Declare the parameters (constants) of the system.
        """
        self.parameter_list = symbols(string)
        return self.parameter_list

    def set_motion_constraint_matrix(self, T_con, T_con_dt, u_dep, u_indep,
            d_ci, i_ci):
        """Set the motion constraint matrix and it's time derivative.

        Returns two substitution dictionaries, one for each matrix.
        """
        m, n = T_con.shape
        assert (m, n) == T_con_dt.shape
        T_con_dict = {}
        T_con_dt_dict = {}
        T_con_sym = zeros((m, n))
        T_con_dt_sym = zeros((m, n))
        for i in range(m):
            for j in range(n):
                tconij = Symbol('T_con%d%d'%(i,j))
                tcondtij = Symbol('T_con%d%dd'%(i,j))
                T_con_dict[tconij] = T_con[i, j]
                T_con_dt_dict[tcondtij] = T_con_dt[i, j]
                T_con_sym[i, j] = tconij
                T_con_dt_sym[i, j] = tcondtij
        self.T_con = T_con_sym
        self.T_con_dt = T_con_dt_sym

        # Set the dependent and independent speeds
        self.u_dependent = u_dep
        self.u_independent = u_indep
        self.udot_dependent = [ud.diff(t) for ud in u_dep]
        self.udot_independent = [ui.diff(t) for ui in u_indep]
        self.dependent_ci = d_ci
        self.independent_ci = i_ci

        return T_con_dict, T_con_dt_dict

    def frstar(self):
        """Computes the generalized inertia forces of the system.
        """
        n = len(self.u_list)
        # p == n in systems where no dependent speeds have been introduced
        p = len(self.u_independent)
        # m == 0 in systems where no dependent speeds have been introduced
        m = len(self.u_dependent)
        self.mass_matrix = zeros((p, p))
        # When dependent speeds are used, the mass matrix must be formed
        # carefully, taking into account the transformation matrix between the
        # independent and the independent speeds:
        # M_i * d(u_i)/dt + M_d * d(u_d)/dt + ... + Fr = 0
        # M_i * d(u_i)/dt + M_d * (d(T_ud)/dt * u_i + T_ud * d(u_i)/dt) + ... + Fr = 0
        # (M_i + M_d * T_ud) * d(u_i)/dt + M_d * d(T_ud)/dt * u_i + ... + Fr = 0
        if m != 0:
            self.mass_matrix_i = zeros((p, p))
            self.mass_matrix_d = zeros((p, m))

        self.recursive_frstar(self.O)
        self.recursive_frstar(self)

    def recursive_frstar(self, PorF):
        """Recursively computes generalized inertia forces for each particle
        and rigid body in the system.

        Generalized inertia forces will be linear in the time derivatives of
        the generalized speeds and the gyroscopic terms of the form u_j*u_k.
        As such, when computing the generalized inertia forces it makes sense
        to collect all like terms, so that simplifications on each of the
        coefficients of these linear terms can be simplified individually.

        We store these coefficients as a Matrix of length: n + (n**2+n)/2.
        Where n is the number of generalized speeds (both dependent and
        independent).  The (n**2 + n)/2 comes from the numer of unique possible
        gyroscopic terms.  The coefficients are ordered from u1,..., un, while
        the gyroscopic terms are ordered according to the ordering of the
        "crossterms" attribute.
        """
        # List of d(u_i)/dt, u_i * u_j, and u_i * qd_j terms.  The acceleration
        # of every point is linear in all these terms.
        udgyro_list = self.udot_list + self.crossterms
        n = len(self.udot_list)
        m = len(self.u_dependent)
        p = len(self.u_independent)
        assert isinstance(PorF, (ReferenceFrame, Point))
        if isinstance(PorF, Point) and PorF.mass == 0:
            PorF.gen_inertia_force = [(0, 0)] * p
        elif isinstance(PorF, ReferenceFrame) and PorF.inertia.dict == {}:
            PorF.gen_inertia_force = [(0, 0)] * p
        else:
            PorF.gen_inertia_force = []
            # Compute the generalized inertia forces
            if isinstance(PorF, Point):
                acc = PorF.abs_acc
                inertia_force = {}
                for k, v in acc.dict.items():
                    inertia_force[k] = -PorF.mass * v
                inertia_force = Vector(inertia_force).expandv()
            else:
                alph = PorF.ang_acc()
                I = PorF.inertia
                w = PorF.ang_vel()
                #inertia_force = (-alph.dot(I)-w.cross(I.rdot(w))).expandv()
                inertia_force = (-dot(I, alph) - cross(w, dot(I, w))).expandv()

            # List of coefficients of all linear terms
            coef_list = inertia_force.partials(udgyro_list)
            # Loop through all partial velocities / partial angular velocites
            for i, pv in enumerate(PorF.partialv):
                sum_ud = 0
                sum_gyro = 0
                coef_list_d_pv = []
                for c in coef_list:
                    coef_list_d_pv.append(c.dot(pv))
                if n == p:  # Case for no motion constraints
                    for j, udot in enumerate(udgyro_list[:p]):
                        self.mass_matrix[i, j] += coef_list_d_pv[j]
                        sum_ud += coef_list_d_pv[j] * udot
                    for j, gyro in enumerate(udgyro_list[p:]):
                        sum_gyro += coef_list_d_pv[j+p] * gyro
                    PorF.gen_inertia_force.append((sum_ud, sum_gyro))
                else:       # Case for systems with motion constraints
                    mm_row = zeros((1, p+m))
                    mm_i_row = zeros((1, p))
                    mm_d_row = zeros((1, m))
                    for j in range(p):
                        mm_row[j] += coef_list_d_pv[j]
                    for j, jt in enumerate(self.dependent_ci):
                        mm_d_row[j] = mm_row[jt]
                    for j, jt in enumerate(self.independent_ci):
                        mm_i_row[j] = mm_row[jt]

                    # Normal gyroscopic terms that appear in GIF's
                    for j, gyro in enumerate(udgyro_list[p:]):
                        sum_gyro += coef_list_d_pv[j+p] * gyro
                    # Extra gyroscopic terms that appear in GIF's due to
                    # constraints
                    sum_gyro += (mm_d_row * self.T_con_dt * Matrix(self.u_independent))[0]

                    # Mass matrix, constrained
                    mm_con = mm_i_row +  mm_d_row*self.T_con
                    sum_ud = (mm_con * Matrix(self.udot_independent))[0]
                    self.mass_matrix[i, :] += mm_con
                    PorF.gen_inertia_force.append((sum_ud, sum_gyro))

        #  Initiate recursion
        if PorF.children == []:
            return
        else:
            for child in PorF.children:
                self.recursive_frstar(child)

    def fr(self):
        """Computes the generalized active forces of the system.
        """
        self.recursive_fr(self.O)
        self.recursive_fr(self)

    def recursive_fr(self, PorF):
        """Recursively computes generalized active forces for each particle
        and rigid body in the system.
        """
        if isinstance(PorF, Point):
            if PorF.force == Vector(0):
                PorF.gen_active_force = [0] * len(self.u_independent)
            else:
                PorF.gen_active_force = [PorF.force.dot(pv) for pv in
                    PorF.partialv]
        elif isinstance(PorF, ReferenceFrame):
            if PorF.torque == Vector(0):
                PorF.gen_active_force = [0] * len(self.u_independent)
            else:
                PorF.gen_active_force = [PorF.torque.dot(pw) for pw in
                    PorF.partialv]
        else:
            raise NotImplementedError()

        #  Initiate recursion
        if PorF.children == []:
            return
        else:
            for child in PorF.children:
                self.recursive_fr(child)

    def gravity(self, v):
        """Applies a gravitational force to each particle and rigid body in the
        system.
        """
        v = Vector(v)
        self.recursive_gravity(self.O, v)

    def recursive_gravity(self, Point, v):
        """Recursively apply gravity to all points which have been assigned a
        nonzero mass."""

        gf = {}
        for k in v.dict:
            gf[k] = Point.mass * v.dict[k]
        Point.force += Vector(gf)
        #  Initiate recursion
        if Point.children == []:
            return
        else:
            for child in Point.children:
                self.recursive_gravity(child, v)

    def form_kanes_equations(self):
        """Forms Kanes equations in a slightly modified form.

        Rather than returning:
        Fr + Fr* = 0

        It returns a list of equations which have the udot's on the left hand
        side and everything else on the opposite side.
        """
        # Form Partial Velocities and Partial Angular velocities for every
        # Point and Reference Frame
        self.recursive_partials(self)
        self.recursive_partials(self.O)
        # Compute the generalized active forces
        self.fr()
        # Compute the generalized inertia forces
        self.frstar()
        p = len(self.u_independent)
        self.kanes_equations = []
        self.ke_lhs = [0] * p
        self.ke_rhs_if = [0] * p
        self.ke_rhs_af = [0] * p
        ke = []
        for i in range(p):
            self.kanes_equations.append([0,0])
        self.recursive_eoms(self.O)
        self.recursive_eoms(self)
        for i in range(p):
            s = S(0)
            for j, ud in enumerate(self.udot_independent):
                c_ud = self.mass_matrix[i, j]
                if c_ud != 0:
                    #if c_ud.could_extract_minus_sign():
                    #    s -= -c_ud * ud
                    #else:
                    s += c_ud * ud
            self.ke_lhs[i] = s
            s = S(0)
            for uu in self.crossterms:
                c_uu = self.ke_rhs_if[i].coeff(uu)
                if c_uu is not None:
                    #if c_uu.could_extract_minus_sign():
                    #    s -= -c_uu * uu
                    #else:
                    s += c_uu * uu
            self.ke_rhs_if[i] = s
            af = factor(self.ke_rhs_af[i].subs(self.trig_subs_dict).\
                    subs(self.symbol_dict)).subs(self.symbol_dict_back).\
                    subs(self.trig_subs_dict_back)
            ke.append(Eq(self.ke_lhs[i], self.ke_rhs_if[i] + af))
        #kes = []
        #for i in range(p):
        #    lhs = collect(self.kanes_equations[i][0], self.udot_independent)
        #    kes.append(Eq(lhs, self.kanes_equations[i][1]))
        self.kanes_equations = ke
        return ke

    def set_kanes_equations(self, eqns):
        self.kanes_equations = eqns

    def solve_kanes_equations(self, dummy_vars=None):
        """Solves Kane's equations for the time derivatives of the generalized
        speeds.

        Forms the adjugate matrix, factors out common terms along the rows,
        performs the matrix multiplication, and then multiplies by the common
        terms and divides by the determinant.

        If optional dummy_mass_matrix argument is eqaul to True, returns the
        result with dummy symbols, and a dictionary with the dummy symbols as
        keys and the expression they represent as the corresponding values.
        """
        m, n = self.mass_matrix.shape
        assert m == n
        mm, mm_dict = dummy_matrix(self.mass_matrix, 'M')
        ke_rhs, ke_dict = dummy_matrix(Matrix([ke.rhs for ke in
            self.kanes_equations]), 'rhs')
        mm_dict.update(ke_dict)
        assert ke_rhs.shape == (n, 1)

        # Form the adjugate and the determinant
        mm_adj = mm.adjugate().expand()
        for i in range(m):
            for j in range(n):
                if mm_adj[i, j] != 0:
                    mm_adj[i, j] = factor(mm_adj[i, j])
        mm_det = factor(mm.det(method="berkowitz").expand())
        assert mm_det != 0, "Mass matrix is singular!!!"

        soln = []
        # Try to factor each row
        for i in range(m):
            row = mm_adj[i,:]
            row_bool = [False] * n
            for j, entry in enumerate(row):
                if isinstance(entry, Mul) or entry == 0:
                    row_bool[j] = True
            if all(row_bool):       # All factored or are zero
                flag = 0
                for entry in row:
                    if entry == 0:
                        continue
                    elif flag == 0:
                        flag = 1
                        row_factor = set(entry.args)
                    else:
                        row_factor &= set(entry.args)
                if row_factor:
                    f = Mul(*list(row_factor))
                    row /= f
                else:
                    f = 1
                soln.append( (f / mm_det) * ((row * ke_rhs)[0]))
            else:
                soln.append( (1 / mm_det) * ((row * ke_rhs)[0]))

        dyndiffs = []
        for i, udot in enumerate(self.udot_independent):
            rhs = soln[i]
            if dummy_vars == None:
                rhs = rhs.subs(mm_dict).expand()
            dyndiffs.append(Eq(udot, rhs))

        if dummy_vars:
            return dyndiffs, mm_dict
        else:
            return dyndiffs

    def recursive_eoms(self, PorF):
        """Traverse Point and ReferenceFrame tree and sum the generalized
        inertia forces and generalized active forces.
        """
        for r in range(len(self.u_independent)):
            # Term on the left hand side of Fr* = -Fr
            self.kanes_equations[r][0] += PorF.gen_inertia_force[r][0]
            # Term on the right hand side of Fr* = -Fr
            self.kanes_equations[r][1] -=  PorF.gen_active_force[r] + PorF.gen_inertia_force[r][1]
            self.ke_lhs[r] += PorF.gen_inertia_force[r][0]
            self.ke_rhs_if[r] -= PorF.gen_inertia_force[r][1]
            self.ke_rhs_af[r] -= PorF.gen_active_force[r]
        if PorF.children == []:
            return
        else:
            for child in PorF.children:
                self.recursive_eoms(child)

    def output_eoms(self, filename, *args):
        """Output the equations of motion to a file as a function which can be
        integrated by scipy.odeint
        """
        ode_func_string = '# ' + time.asctime() + '\n'
        ode_func_string += "from numpy import sin, cos, tan, vectorize\n\n"
        ode_func_string += "def f(x, t, parameter_list):\n"
        ode_func_string += '    # Unpacking the parameters\n'
        s = ""
        for p in self.parameter_list:
            s += str(p) + ', '
        ode_func_string += '    ' + s[:-2] + ' = ' + 'parameter_list\n'
        ode_func_string += '    # Unpacking the states (q\'s and u\'s)\n'
        s = ""

        for q in self.q_list:
            s += str(q) + ', '
        for u in self.u_independent:
            s += str(u) + ', '
        ode_func_string += '    ' + s[:-2] + ' = ' + 'x\n'

        trig_func_string = ""

        for tf in self.trig_func_set:
            trig_func_string += '    ' + str(tf) + ' = '
            if str(tf)[0] == 's':
                trig_func_string += 'sin(' + str(tf.args[0]) + ')\n'
            if str(tf)[0] == 'c':
                trig_func_string += 'cos(' + str(tf.args[0]) + ')\n'
            if str(tf)[0] == 't':
                trig_func_string += 'tan(' + str(tf.args[0]) + ')\n'

        ode_func_string += trig_func_string

        dxdt_list = ""

        if hasattr(self, 'u_dependent'):
            ode_func_string += '    # Dependent generalized speeds\n'
            for u in self.u_dependent:
                ode_func_string += '    ' + str(u) + ' = ' +\
                    str(self.u_dependent_eqs[u].subs(self.qdot_list_dict)) + '\n'

        ode_func_string += '    # Kinematic differential equations\n'
        qdl = []
        for qd in self.qdot_list:
            if qd in self.kindiffs: qdl.append(qd)

        for qd in qdl:
            ode_func_string += '    ' + str(qd)[:-1] + 'p' + ' = ' + str(self.kindiffs[qd]) + '\n'
            dxdt_list += str(qd)[:-1] + 'p, '

        ode_func_string += '    # Dynamic differential equations\n'

        for ud in self.udot_independent:
            ode_func_string += '    ' + str(ud)[:-1] + 'p' +  ' = ' + str(self.dyndiffs[ud]) + '\n'
            dxdt_list += str(ud)[:-1] + 'p, '

        ode_func_string += '    ' + 'return [' + dxdt_list[:-2] + ']'

        qdot2u_func_string = ""
        qdot2u_func_string += "def qdot2u(q, qd, parameter_list):\n"
        qdot2u_func_string += '    # Unpacking the parameters\n'
        s = ""
        for p in self.parameter_list:
            s += str(p) + ', '
        qdot2u_func_string += '    ' + s[:-2] + ' = ' + 'parameter_list\n'
        qdot2u_func_string += '    # Unpacking the q\'s and qdots\n'
        s = ""
        for q in self.q_list:
            s += str(q) + ', '
        qdot2u_func_string += '    ' + s[:-2] + ' = ' + 'q\n'
        s = ""
        for qd in qdl:
            s += str(qd)[:-1] + 'p, '
        qdot2u_func_string += '    ' + s[:-2] + ' = ' + 'qd\n'

        trig_func_string = ""

        for tf in self.trig_func_set:
            trig_func_string += '    ' + str(tf) + ' = '
            if str(tf)[0] == 's':
                trig_func_string += 'sin(' + str(tf.args[0]) + ')\n'
            if str(tf)[0] == 'c':
                trig_func_string += 'cos(' + str(tf.args[0]) + ')\n'
            if str(tf)[0] == 't':
                trig_func_string += 'tan(' + str(tf.args[0]) + ')\n'

        qdot2u_func_string += trig_func_string

        dxdt_list = ""

        qdot2u_func_string += '    # Kinematic differential equations\n'
        if hasattr(self, 'dependent_rates'):
            qdot_i = [qd.subs(self.qdot_list_dict) for qd in qdl]
            for i, u in enumerate(self.u_list):
                qdot2u_func_string += '    ' + str(u) + ' = ' +\
                        str((self.transform_matrix[i,:]*Matrix(qdot_i))[0]) + '\n'
                dxdt_list += str(u) + ', '
            qdot2u_func_string += '    return [' + dxdt_list[:-2] + ']'
        else:
            for i, u in enumerate(self.u_list):
                qdot2u_func_string += '    ' + str(u) + ' = ' +\
                        str((self.transform_matrix[i,:]*Matrix(self.qdot_list_s))[0]) + '\n'
                dxdt_list += str(u) + ', '
            qdot2u_func_string += '    return [' + dxdt_list[:-2] + ']'

        f = open(filename, 'w')
        f.write(ode_func_string + '\n\n' + qdot2u_func_string)

        if args:
            n = len(args)
            a =  ""
            a += "def animate(q, parameter_list):\n"
            a += '    # Unpacking the parameters\n'
            s = ""
            for p in self.parameter_list:
                s += str(p) + ', '
            a += '    ' + s[:-2] + ' = parameter_list\n'
            a += '    # Unpacking the coordinates\n'
            s = ""
            for q in self.q_list:
                s += str(q) + ', '
            a += '    ' + s[:-2] + ' = q\n'

            trig_func_set = set([])
            a_temp = ""
            ret_string = ""
            for k, arg in enumerate(args):
                ret_string += "["
                if isinstance(arg, (UnitVector, Vector)):
                    pos_or_axis = [arg.dot(self[i]) for i in (1,2,3)]
                    for i, p in enumerate(pos_or_axis):
                        nv = "out_%d_%d"%(k,i)
                        a_temp += "    " + nv + " = " + str(p) + "\n"
                        ret_string += nv + ", "
                        trig_terms = p.atoms(sin, cos, tan)
                        if trig_terms:
                            trig_func_set.update(trig_terms)
                    ret_string = ret_string[:-2] + "], "
                else:
                    raise TypeError('Optional parameters must be Vector/UniVectors')

            a += "    # Trigonometric functions needed\n"
            trig_func_string = ""
            for tf in trig_func_set:
                trig_func_string += '    ' + str(tf) + ' = '
                if str(tf)[0] == 's':
                    trig_func_string += 'sin(' + str(tf.args[0]) + ')\n'
                if str(tf)[0] == 'c':
                    trig_func_string += 'cos(' + str(tf.args[0]) + ')\n'
                if str(tf)[0] == 't':
                    trig_func_string += 'tan(' + str(tf.args[0]) + ')\n'
            a += trig_func_string

            a += "    # Position of Points and Axis/Angle Calculations\n"
            if ret_string != "":
                a_temp += "    return " + ret_string[:-2]
            a += a_temp
            f.write('\n\n' + a)
        f.close()

    def define_speeds(self, eqns):
        """Defines the generalized speeds equations, conditioning them so that
        they are more easily inverted to determined the kinematic differential
        equations.
        """
        eqns_cond = []
        for e in eqns:
            rhs = collect(e.rhs.expand().subs(self.csqrd_dict).expand(), self.qdot_list)
            eqns_cond.append(Eq(e.lhs, rhs))
        return eqns_cond


def matrixv_multiply(A, B):
    """For multplying a matrix of PyDy Vector/UnitVectors with matrices of
    Sympy expressions.

    Normal matrix_multiply doesn't work because PyDy vectors are not derived
    from Basic."""
    ma, na = A.shape
    mb, nb = B.shape
    if na != mb:
        raise ShapeError()
    product = Matrix(ma, nb, lambda i,j: 0)
    for i in xrange(ma):
            for j in xrange(nb):
                s = Vector(0)
                for k in range(na):
                    aik = A[i, k]
                    bkj = B[k, j]
                    if isinstance(aik, Vector):
                        assert not isinstance(bkj, (UnitVector, Vector))
                        p = {}
                        for uv, val in aik.dict.items():
                            p[uv] = bkj*val
                    elif isinstance(aik, UnitVector):
                        assert not isinstance(bkj, (UnitVector, Vector))
                        p = bkj*aik
                    elif isinstance(bkj, Vector):
                        assert not isinstance(aik, (UnitVector, Vector))
                        p = {}
                        for uv, val in bkj.dict.items():
                            p[uv] = aik*val
                    elif isinstance(bkj, UnitVector):
                        assert not isinstance(aik, (UnitVector, Vector))
                        p = aik*bkj
                    else:
                        raise NotImplementedError()
                    s += Vector(p)
                product[i, j] = s
    return product

def most_frequent_frame(vector):
    """Determines the most frequent frame of all unitvector terms in a vector.
    """
    frame_counter = {}
    for uv in vector.dict:
        frame_counter[uv.frame] = frame_counter.get(uv.frame, 0) + 1
    return max([(frame_counter[x], x) for x in frame_counter])[1]

def express(v, frame):
    """Expresses a vector in terms of UnitVectors fixed in a specified frame.
    """
    v = Vector(v)
    if (isinstance(v, UnitVector) or isinstance(v, Vector)) and \
            (isinstance(frame, ReferenceFrame)):
        return v.express(frame)
    else:
        raise TypeError('v must be UnitVector or Vector object and frame must \
            be a ReferenceFrame object')

def dot(v1, v2):
    """Vector dot product.

    between UnitVector, Vector, and Dyad classes

    Returns a scalar sympy expression in the case of the dot product between
    two UnitVectors/Vectors.  Returns a UnitVector/Vector in the case of
    the dot product between a Dyad and a UnitVector/Vector.

    In the scalar dot product, the operation commutes, i.e. dot(v1, v2) dot(v2,
    v1).  In the vector/dyad dot product, the operation is noncommutative,
    i.e., dot(v1, v2) != dot(v2, v1)

    """

    if isinstance(v1, (UnitVector, Vector)) and isinstance(v2, (UnitVector,
        Vector)):
        return v1.dot(v2)
    elif isinstance(v1, Dyad) and isinstance(v2, (UnitVector, Vector)):
        return v1.rdot(v2)
    elif isinstance(v2, Dyad) and isinstance(v1, (UnitVector, Vector)):
        return v2.ldot(v1)
    else:
        if not isinstance(v1, (UnitVector, Vector)):
            v1 = Vector(v1)
        if not isinstance(v2, (UnitVector, Vector)):
            v2 = Vector(v2)
        return v1.dot(v2)

def cross(v1, v2):
    """Vector cross product.

    Parameters
    v1, v2: PyDy UnitVector or Vector objects.

    Returns
    A UnitVector or Vector object.

    See Also
    L{dot}, L{express}
    """

    if (isinstance(v1, UnitVector) or isinstance(v1, Vector)) and \
            (isinstance(v2, UnitVector) or isinstance(v2, Vector)):
                return v1.cross(v2)
    else:
        if not (isinstance(v1, UnitVector) or isinstance(v1, Vector)):
            v1 = Vector(v1)
        if not (isinstance(v2, UnitVector) or isinstance(v2, Vector)):
            v2 = Vector(v2)
        return v1.cross(v2)

def coeffv(v, scalar):
    if isinstance(v, Vector):
        return v.coeffv(scalar)
    elif isinstance(v, UnitVector):
        return S(1)
    else:
        raise NotImplementedError()

def dt(v, frame):
    """Time derivative of a vector as viewed by an observer fixed in a frame.
    """
    v = Vector(v)
    if isinstance(frame, ReferenceFrame):
        if isinstance(v, (UnitVector, Vector)):
            res = v.dt(frame)
            return res
        else:
            raise TypeError('First argument must be a Vector or \
            UnitVector, instead a %s object was given' % str(type(v)))
    else:
        raise TypeError('Second argument must be a ReferenceFrame, \
                instead a %s object was given' % str(type(v)))

def GeneralizedCoordinate(s, constant=False):
    gc = Symbol(s)(Symbol('t'))
    gc.is_gc = True
    if constant==True:
        gc.fdiff = lambda argindex: 0
    gc.__repr__ = lambda self: PyDyStrPrinter().doprint(self)
    gc.__str__ = lambda self: PyDyStrPrinter().doprint(self)
    return gc

def gcs(s, number=1, list=False):
    gc_list = [GeneralizedCoordinate(s[0]+str(i)) for i in range(1, number+1)]
    if list == False:
        if number == 1:
            return gc_list[0]
        else:
            return gc_list
    elif list == True:
        gcd_list = [gc.diff(t) for gc in gc_list]
        return (gc_list, gcd_list)

class PyDyStrPrinter(StrPrinter):
    #printmethod = '_sympystr_'
    def _print_UnitVector(self, e):
        s = str(e.v['sym'])
        name = s[:-1]
        index = s[-1]
        r = "%s%s>" % (name, index)
        return r

    def _print_Vector(self, e):
        s = ''
        i = 0
        small_dot = "*"
        if e.dict.keys() != []:
            uv_list = e.dict.keys()
            uv_list.sort(sort_UnitVector)
            for k in uv_list:
                # Case when the scalar coefficient is 1 or -1
                if (e.dict[k] == 1) or (e.dict[k] == -1):
                    # First term don't print a leading + if positive
                    if i == 0:
                        if e.dict[k] == 1: sign = ''
                        if e.dict[k] == -1: sign = '-'
                        s += sign + self.doprint(k)
                        i += 1
                    # All other terms put the sign and pad with spaces
                    else:
                        if e.dict[k] == 1: sign = '+'
                        if e.dict[k] == -1: sign = '-'
                        s += ' ' + sign + ' ' + self.doprint(k)
                # Case when the scalar coefficient is a Sympy expression
                else:
                    # First term
                    if i == 0:
                        # Put parenthesis around Add terms
                        if isinstance(e.dict[k], Add):
                            s += ('(' + self.doprint(e.dict[k]) +
                                    ')' + small_dot + self.doprint(k))
                        else:
                            s += (self.doprint(e.dict[k]) +
                                small_dot + self.doprint(k))
                        i += 1
                    # All subsequent terms pad with spaces and add parenthesis
                    else:
                        if isinstance(e.dict[k], Add):
                            if e.dict[k].expand().could_extract_minus_sign():
                                s += (' - (' + self.doprint(-e.dict[k].expand())
                                    + ')' + small_dot + self.doprint(k))
                            else:
                                s += (' + (' + self.doprint(e.dict[k])
                                    + ')' + small_dot + self.doprint(k))
                        elif isinstance(e.dict[k], Mul):
                            mulcoef = S(1)
                            sign_counter = 0
                            for arg in e.dict[k].args:
                                if arg.expand().could_extract_minus_sign():
                                    mulcoef *= -arg.expand()
                                    sign_counter += 1
                                else:
                                    mulcoef *= arg
                            if sign_counter % 2 == 0:
                                s += ' + ' + (self.doprint(mulcoef) + small_dot +
                                        self.doprint(k))
                            else:
                                s += ' - ' + (self.doprint(mulcoef) + small_dot +
                                        self.doprint(k))
                        elif e.dict[k].is_negative:
                            s += ' - ' + (self.doprint(-e.dict[k]) + small_dot
                                    + self.doprint(k))
                        else:
                            s += ' + ' + (self.doprint(e.dict[k]) + small_dot +
                                self.doprint(k))
            return s
        else:
            return "0>"

    def _print_Function(self, e):
        """
        Print ui(t) as ui, where is i is an index number of the generalized
        speed.
        """
        if hasattr(e, 'is_gc'):
            return str(e.func)
        else:
            return StrPrinter().doprint(e)

    def _print_Symbol(self, e):
        """
        Print ui(t) as ui, where is i is an index number of the generalized
        speed.
        """
        if hasattr(e, 'is_gc'):
            return str(e.func)
        else:
            return StrPrinter().doprint(e)

    def _print_sin(self, e):
        """
        Print sin(qi(t)) as si, where i is any number.
        """
        if str(e.args[0].func)[0] == 'q':
            return 's' + str(e.args[0].func)[1:]
        else:
            return StrPrinter().doprint(e)

    def _print_cos(self, e):
        """
        Print cos(qi(t)) as si, where i is any number.
        """
        if str(e.args[0].func)[0] == 'q':
            return 'c' + str(e.args[0].func)[1:]
        else:
            return StrPrinter().doprint(e)

    def _print_tan(self, e):
        """
        Print tan(qi(t)) as si, where i is any number.
        """
        if str(e.args[0].func)[0] == 'q':
            return 't' + str(e.args[0].func)[1:]
        else:
            return StrPrinter().doprint(e)

    def _print_Derivative(self, expr):
        if len(expr.args) == 2:
            return str(expr.args[0].func) + "d"*len(expr.args[1:])
        elif len(expr.args) == 3:
            return str(expr.args[0].func) + "d"*len(expr.args[1:])
        else:
            return StrPrinter().doprint(expr)

    def _print_Dyad(self, expr):
        s = ""
        if expr.dict == {}:
            return "0>>"
        for k, v in expr.dict.items():
            if isinstance(v, Add):
                v_str = "(" + str(v) + ")"
            else:
                v_str = str(v)
            if k.is_Mul:
                s += v_str + "*" + str(k.args[0]) + "*" + str(k.args[1]) + " + "
            elif k.is_Pow:
                s += v_str + "*" + str(k.args[0]) + "*" + str(k.args[0]) + " + "
        return s[:-3]

    #def _print_Matrix(self, expr):
    #    return expr._format_str(lambda elem: elem.__str__())

class PyDyPrettyPrinter(PrettyPrinter):
    def _print_UnitVector(self, e):
        class Fake(object):
            baseline = 0
            def render(self, *args, **kwargs):
                one = "\xe2\x82\x81"
                two = "\xe2\x82\x82"
                three = "\xe2\x82\x83"
                bold = "\033[1m"
                reset = "\033[0;0m"
                s = str(e.v['sym'])
                name = s[:-1]
                index = s[-1]
                r = "%s%s" % (bold, name)
                if index == "1":
                    r += one
                elif index == "2":
                    r += two
                elif index == "3":
                    r += three
                r += reset
                return r
        return Fake()

    def _print_Vector(self, e):
        class Fake(object):
            def render(self, *args, **kwargs):
                s = ''
                i = 0
                small_dot = "\xC2\xB7"
                if e.dict.keys() != []:
                    uv_list = e.dict.keys()
                    uv_list.sort(sort_UnitVector)
                    for k in uv_list:
                        # Case when the scalar coefficient is 1 or -1
                        if (e.dict[k] == 1) or (e.dict[k] == -1):
                            # First term don't print a leading + if positive
                            if i == 0:
                                if e.dict[k] == 1: sign = ''
                                if e.dict[k] == -1: sign = '-'
                                s += sign + ppuv(k)
                                i += 1
                            # All other terms put the sign and pad with spaces
                            else:
                                if e.dict[k] == 1: sign = '+'
                                if e.dict[k] == -1: sign = '-'
                                s += ' ' + sign + ' ' + ppuv(k)
                        else:
                            # First term
                            if i == 0:
                                # Put parenthesis around Add terms
                                if isinstance(e.dict[k], Add):
                                    s += ('(' + pretty(e.dict[k]) +
                                            ')' + small_dot + ppuv(k))
                                else:
                                    s += (pretty(e.dict[k]) +
                                        small_dot + ppuv(k))
                                i += 1
                            # All other terms pad with spaces and add parenthesis
                            else:
                                if isinstance(e.dict[k], Add):
                                    s += (' + (' + pretty(e.dict[k])
                                        + ')' + small_dot + ppuv(k))
                                elif isinstance(e.dict[k], (Mul, Pow)):
                                    coef = (pretty(e.dict[k]) + small_dot +
                                            ppuv(k))
                                    if coef[0] == '-':
                                        s += ' - ' + coef[1:]
                                    else:
                                        s += ' + ' + coef
                                else:
                                    s += (' + ' + pretty(e.dict[k]) + small_dot +
                                            ppuv(k))

                    return s
                else:
                    return "\033[1m" + "0" + "\033[0;0m"
        return Fake()

    def _print_Derivative(self, expr):
        return str(expr.args[0].func) + "'"*len(expr.args[1:])

    def _print_Mul(self, e):
        s = ''
        i = 0
        e_ts = trigsimp(expand(trigsimp(e)))
        if e_ts.is_Mul:
            N = len(e_ts.args)
            for a in e_ts.args:
                if i == 0:
                    if a == -1:
                        s += '-'
                    else:
                        s += self.doprint(a) + '\xC2\xB7'
                    i += 1
                elif i < N-1:
                    #if a.is_Pow and a.args[1]<0:
                    #    s += '\b/' + self.doprint(a.args[0]) + '\xC2\xB7'
                    #else:
                    s += self.doprint(a) + '\xC2\xB7'
                    i += 1
                else:
                    s += self.doprint(a)
            return s
        else:
            return self.doprint(e_ts)

    def _print_sin(self, e):
        """
        Print sin(qi(t)) as si, where i is any number.
        """
        class Fake(object):
            def render(self, *args, **kwargs):
                if str(e.args[0].func)[0] == 'q':
                    return u's' + unicode_subscript(str(e.args[0].func)[1:])
                else:
                    return PrettyPrinter().doprint(e)
        return Fake()

    def _print_cos(self, e):
        """
        Print cos(qi(t)) as si, where i is any number.
        """
        class Fake(object):
            def render(self, *args, **kwargs):
                if str(e.args[0].func)[0] == 'q':
                    return u'c' + unicode_subscript(str(e.args[0].func)[1:])
                else:
                    return PrettyPrinter().doprint(e)
        return Fake()

    def _print_tan(self, e):
        """
        Print tan(qi(t)) as si, where i is any number.
        """
        class Fake(object):
            def render(self, *args, **kwargs):
                if str(e.args[0].func)[0] == 'q':
                    return u't' + unicode_subscript(str(e.args[0].func)[1:])
                else:
                    return PrettyPrinter().doprint(e)
        return Fake()

def unicode_subscript(num):
    """Converts an integer to the unicode subscript representation of that
    integer.

    Reference

    """
    n = str(num)
    subscript_dict = {
            '0': u'\u2080',
            '1': u'\u2081',
            '2': u'\u2082',
            '3': u'\u2083',
            '4': u'\u2084',
            '5': u'\u2085',
            '6': u'\u2086',
            '7': u'\u2087',
            '8': u'\u2088',
            '9': u'\u2089'}
    uni = u""
    for u in n:
        uni += subscript_dict[u]
    return uni

def pprint(e):
    print PyDyPrettyPrinter().doprint(e)

def ppuv(e):
    """Pretty print a UnitVector.
    """
    one = "\xe2\x82\x81"
    two = "\xe2\x82\x82"
    three = "\xe2\x82\x83"
    bold = "\033[1m"
    reset = "\033[0;0m"
    s = str(e.v['sym'])
    name = s[:-1]
    index = s[-1]
    r = "%s%s" % (bold, name)
    if index == "1":
        r += one
    elif index == "2":
        r += two
    elif index == "3":
        r += three
    r += reset
    return r

def sort_UnitVector(a, b):
    """Sort UnitVector objects by how many rotations their reference frame is away from
    the Newtonian frame.
    """
    if a.frame == b.frame:
        return (a.i > b.i) - (a.i < b.i)
    else:
        return (len(a.frame.ref_frame_list) > len(b.frame.ref_frame_list)) -\
            (len(a.frame.ref_frame_list) < len(b.frame.ref_frame_list))

def coefficient_matrix(eqns, linear_terms):
    """Given a list of equations linear in some specified terms, form the
    matrix of coefficients of those linear terms.

    """
    m = len(eqns)
    n = len(linear_terms)
    B = zeros((m, n))
    d = {}
    for i in range(m):
        for j in range(n):
            B_ij = eqns[i].expand().coeff(linear_terms[j])
            if B_ij is not None:
                B[i, j] = B_ij
    return B

def generate_function(name, Eq_list, func_args, params=None, nested_terms=None,
        docstring=None, triples=None, time=None):
    """Generate a Python function in string form.

    Input:
        name:       A string name for the function to be generated
        Eq_list:    A list of Sympy relational objects (lhs == rhs)
        func_args:  Quantities neccessary to compute all the right hand sides
            of the equations in Eq_list.
        params:     A list of Sympy Symbols/Functions which are to be treated
                    as function parameters
        nested_terms:   Quantities which appear in the right hand sides of the
                        equations in Eq_list and should be evaluated before
                        those expressions.
        docstring:  A string to be used as the functions docstring.
        triples:    Boolean value which will cause the return values to be
                    returned as a list of length 3 lists.
        time:       Boolean which will cause time to be an arugment in the
                    function signature.

    Output:
        A Python string which is exec()-able and defines a function that
        calculates the right hand sides of the expressions in Eq_listand returns
        them as a list.

    """

    fs = ""
    if time:
        time_string = ", t"
    else:
        time_string = ""
    if params:
        fs += "def " + name + "(_x" + time_string + ", _params):\n"
    else:
        fs += "def " + name + "(_x" + time_string + "):\n"

    if docstring:
        fs += '    """' + docstring + '\n    """\n'
    arg_string = "    "
    for a in func_args:
        if isinstance(a, (Symbol, Derivative)):
            arg_string += str(a) + ", "
        elif isinstance(a, Function):
            arg_string += str(a.func) + ", "
    arg_string = arg_string[:-2] + " = _x\n"
    fs += "    # Unpack function arguments\n"
    fs += arg_string
    if params:
        param_string = "    "
        for p in params:
            if isinstance(p, Symbol):
                param_string += str(p) + ", "
            elif isinstance(p, Function):
                param_string += str(p.func) + ", "
        param_string = param_string[:-2] + " = _params\n"
        fs += "\n    # Unpack function parameters\n"
        fs += param_string

    m = len(Eq_list)
    # Trig terms
    trig_set = set([])
    for eqn in Eq_list:
        for i in range(m):
            trig_set.update(eqn.rhs.atoms(sin, cos, tan))
        if nested_terms:
            for nest in nested_terms:
                for v in nest.values():
                    trig_set.update(v.atoms(sin, cos, tan))
    trig_set = list(trig_set)
    trig_set.sort()
    if trig_set:
        trig_string = "\n    # Trigonometric functions\n"
        for tt in trig_set:
            trig_string += "    " + str(tt) + " = " + str(type(tt)) + "(" + str(tt.args[0]) + ")\n"
        fs += trig_string

    # Nested terms
    if nested_terms:
        nested_string = "\n    # Nested terms\n"
        for nest in nested_terms:
            ntk = nest.keys()
            ntk.sort()
            for nt in ntk:
                nested_string += "    " + str(nt) + " = " + str(nest[nt]) + "\n"
        fs += nested_string
    ret_string = "    return ["
    fs += "\n    # Calculate return values\n"
    if triples:
        ret_string_d = ""
        i = 1
        for eqn in Eq_list:
            fs += "    " + str(eqn.lhs) + " = " + str(eqn.rhs) + "\n"
            ret_string_d += str(eqn.lhs) + ", "
            if i % 3 == 0:
                ret_string += "[" + ret_string_d[:-2] + "], "
                ret_string_d = ""
                i = 1
                continue
            i += 1
    else:
        for eqn in Eq_list:
            fs += "    " + str(eqn.lhs) + " = " + str(eqn.rhs) + "\n"
            ret_string += str(eqn.lhs) + ", "
    fs += "\n    # Return calculated values\n"
    fs += ret_string[:-2] + "]\n\n"

    return fs



def linear_transform(B, params, name, det=None, nested_terms=None, x=None,\
        y=None, docstring=None):
    """Given a m x n matrix of Sympy expressions, return an exec-able string
    which would define the Python function mapping x \in R^n to y \in R^m.

    Required arguments:
        B:  A Sympy matrix

        params: A list of Symbol or Function objects upon which the entries of
        B depend.  The order of the list will govern the order of the function
        signature.

        name:  The desired name of the automatically generated function.

    Optional arguments:
        det:  When matrix inverses are constructed by forming the adjugate
        matrix and the determinant, perform the matrix multiplication first,
        then divide each element by the determinant.

        nested_terms:  When the entries of B have been defined in terms of
        quantities which are not in the parameter list, but depend upon
        quantities in the parameter list.

    Returns:
        A string with the function signature:
            def name(x, params):
                ...
                return B*x

    """

    fs = ""
    fs += "def " + name + "(_x, _params):\n"
    if docstring:
        fs += '    """' + docstring + '\n    """\n'
    m,n = B.shape
    param_string = ""
    for p in params:
        if isinstance(p, Symbol):
            param_string += str(p) + ", "
        elif isinstance(p, Function):
            param_string += str(p.func) + ", "
        elif isinstance(p, Derivative) and len(p.args) == 2:
            param_string += str(p.args[0]) + "p, "
    param_string = param_string[:-2] + " = _params\n"

    x_string = ""
    if x:
        x_var = []
        for j in range(n):
            if str(x[j])[-1] == "'":
                x_var.append(Symbol(str(x[j])[:-1] + "p"))
                x_string += str(x_var[-1]) + ", "
            else:
                x_var.append(x[j])
                x_string += str(x[j]) + ", "
        x_var = Matrix(x_var)
    else:
        x_var = Matrix([n, 1], lambda i,j: Symbol("_x%d"%j))
        for j in range(n):
            x_string += "_x%d"%j + ", "
    x_string = x_string[:-2] + " = _x\n"

    fs += "    " + x_string
    fs += "    " + param_string

    # Trig terms
    trig_set = set([])
    for i in range(m):
        for j in range(n):
            trig_set.update(B[i, j].atoms(sin, cos, tan))
    if nested_terms:
        for nest in nested_terms:
            for v in nest.values():
                trig_set.update(v.atoms(sin, cos, tan))

    if trig_set:
        trig_string = ""
        for tt in trig_set:
            trig_string += "    " + str(tt) + " = " + str(type(tt)) + "(" + str(tt.args[0]) + ")\n"
        fs += trig_string

    # Nested terms
    if nested_terms:
        nested_string = ""
        for nest in nested_terms:
            ntk = nest.keys()
            ntk.sort()
            for nt in ntk:
                nested_string += "    " + str(nt) + " = " + str(nest[nt]) + "\n"
        fs += nested_string

    if det:
        fs += "    det = " + str(det) + "\n"

    # Perform the matrix multiplication
    ret_string = "    return ["
    for i in range(m):
        if y:
            if str(y[i])[-1] == "'":
                fs += "    " + str(y[i])[:-1] + "p = "
                ret_string += str(y[i])[:-1] + "p, "
            else:
                fs += "    " + str(y[i]) + " = "
                ret_string += str(y[i]) + ", "
        else:
            fs += "    _y%d = "%i
            ret_string += "_y%d, "%i
        if det:
            fs += "("
        #row = B[i, :]
        #for j in range(n):
        #    Bij = simplify(row[j])
        #    n, d = Bij.as_numer_denom()
        from sympy import together
        prod = together((B[i, :]*x_var)[0])
        #num, den = prod.as_numer_denom()
        #num = factor(num)
        #den = factor(den)
        #prod = num / den
        fs += str(prod)
        if det:
            fs += ")/det"
        fs += "\n"
    ret_string = ret_string[:-2] + "]\n\n"
    fs += ret_string

    return fs

def transform_matrix(B, x, x_dependent, subs_dict=None):
    """Given an m x n coefficent matrix B, n linear terms x, and m linear terms
    xd taken to be dependent, return the transform matrix between the
    independent linear terms and the dependent ones.

    Given:

        B*x = 0

    Where:
        B \in R^{m x n}
        x \in R^{n x 1}
        0 \in R^{m x 1}

    we can partition x into dependent and indepent terms and rewrite as:

        Bd*xd + Bi*xi = 0

    Where:
        Bd \in R^{m x m}
        Bi \in R^{m x (n-m)}
        xd \in R^{m x 1}
        xi \in R^{(n-m) x 1}

    so:

    xd = -inv(Bd)*Bi*xi
       = T*xi

    Returns: -inv(Bd), Bi, dependent_index, independent_index

    """
    m, n = B.shape
    md = len(x_dependent)
    if m != md:
        raise ValueError('Number of equations must equal number of ' +
                         'dependent terms.')
    independent_ci = []
    dependent_ci = []
    try:
        for xd in x_dependent:
            dependent_ci.append(x.index(xd))
        dependent_ci.sort()
        independent_ci = list(set(range(n)) - set(dependent_ci))
        independent_ci.sort()
    except ValueError:
        print('Each of the dependent speeds must be in the speed list used' +
              'in the declare_speeds.')

    # Create a matrix with dummy symbols representing non-zero entries
    B_dummy, d = dummy_matrix(B, 'b')

    # Generate the independent and dependent matrices
    Bd = zeros((m, m))
    Bi = zeros((m, n-m))
    for j, jd in enumerate(dependent_ci):
        Bd[:, j] = B_dummy[:, jd]
    for j, ji in enumerate(independent_ci):
        Bi[:, j] = B_dummy[:, ji]

    # Invert the Bd matrix and determine
    # xd = -inv(Bd) * Bi * xi = T * xi
    # Form the adjugate and matrix multiply by Bi
    Bd_adj = Bd.adjugate().expand()

    # Form the negative of the determinant
    Bd_det = -factor(Bd.det().expand())
    assert Bd_det != 0, "Equations are singular."

    # Form inv(Bd)
    # inv(Bd) = adjugate(Bd) / det(Bd)
    Bd_inv = zeros((m,m))
    for i in range(m):
        for j in range(m):
            if Bd_adj[i,j] != 0:
                Bd_inv[i,j] = factor(Bd_adj[i,j]) / Bd_det
    if subs_dict==None:
        Bd_inv = Bd_inv.subs(d)
        Bi = Bi.subs(d)
        return Bd_inv, Bi
    else:
        return Bd_inv, Bi, d

def eqn_list_to_dict(eqn_list, reverse=None):
    """Convert a list of Sympy Relational objects to a dictionary.

    Most commonly, this will be for converting things like:

    [y == a*x + b, z == c*x + d]

    to a Python dictionary object with keys corresponding to the left hand side
    of the relational objects and values corresponding to the right hand side
    of the relational objects:

    {y: a*x+b, z: c*x + d}

    Optional argument 'reverse' swaps the keys for the values, i.e:

    {a*x + b: y, c*x +d: z}

    Remember that dictionaries are *NOT* ordered, so if the list you pass
    requires a special ordering, it will *NOT* be presever by converting it to
    a dictionary.

    """
    eqn_dict = {}
    for eqn in eqn_list:
        if reverse:
            eqn_dict[eqn.rhs] = eqn.lhs
        else:
            eqn_dict[eqn.lhs] = eqn.rhs

    return eqn_dict

def dict_to_eqn_list(eqn_dict, reverse=None):
    """Convert a Python dictionary to a list of Sympy Relational objects.
    """
    eqn_list = []
    keys = eqn_dict.keys()
    keys.sort()
    for k in keys:
        if reverse:
            eqn_list.append(Eq(eqn_dict[k], k))
        else:
            eqn_list.append(Eq(k, eqn_dict[k]))

    return eqn_list

def animate(Frame, *args):
    """Generate a list of equations useful for creating animations.

    """
    n1 = Frame[1]
    n2 = Frame[2]
    n3 = Frame[3]
    eqn_list = []
    for a in args:
        for i, ni in enumerate((n1, n2, n3)):
            eqn_list.append(Eq(Symbol(a[0]+"_%d"%(i+1)), dot(a[1], ni)))

    return eqn_list

def mass_center(O, points):
    """Calculate the mass center of a list of points relative to the point O.

    The position of each point in the list, relative to the point O, must be
    defined.

    The points list can either be of the form:
    [P1, P2, ..., Pn]
    or
    [(P1, m1), (P2, m2), ..., (Pn, mn)]

    The second form is useful when you want to form the center of mass of the
    system and not assign mass to individual points.
    """
    assert isinstance(O, Point), "First argument must be a Point object"
    mt = S(0)
    cm = {}
    for p in points:
        if isinstance(p, Point):
            pi = p.rel(O)       # Position vector from O to P_i
            mi = p.mass         # Mass of point P_i
        elif isinstance(p, tuple):
            pi = p[0].rel(O)
            mi = p[1]
        # Compute the total mass of all the points/particles in the given list.
        mt += mi
        for k, v in pi.dict.items():
            cm[k] = cm.get(k, S(0)) + mi*v
        # Divide UnitVector coefficient by the total mass
    for k in cm:
        cm[k] /= mt
    return Vector(cm)

def inertia_of_point_mass(m, p, F):
    """Determine the Inertia dyadic of a particle.

    Input:
        m:  Mass
        p:  Position from point O to mass m
        F:  Reference frame to express the dyad entries with respect to.

    Output:
        Inertia Dyadic relative to O of a particle of mass m, located relative
        to the point O by the position vector p.

    """
    I11 = m*dot(cross(p, F[1]), cross(p, F[1]))
    I22 = m*dot(cross(p, F[2]), cross(p, F[2]))
    I33 = m*dot(cross(p, F[3]), cross(p, F[3]))
    I12 = m*dot(cross(p, F[1]), cross(p, F[2]))
    I23 = m*dot(cross(p, F[2]), cross(p, F[3]))
    I13 = m*dot(cross(p, F[1]), cross(p, F[3]))
    return Inertia(F, [I11, I22, I33, I12, I23, I13])

def dummy_matrix(mat, char):
    """Returns a matrix of dummy symbols for non-zero and non-unity entries.

    char specifies the string to use in the beginning of the names of the dummy
    symbols.

    Also returns a substitution dictionary with the dummy symbols as the
    keys and the symbolic expression they represent as the values.
    """
    m, n = mat.shape
    new_mat = zeros((m, n))
    d = {}
    dr = {}
    for i in range(m):
        for j in range(n):
            mij = mat[i, j]
            if mij != 0:
                if mij == 1 or mij == -1:
                    new_mat[i, j] = mij
                    continue
                if mij in dr:
                    new_mat[i, j] = dr[mij]
                elif -mij in dr:
                    new_mat[i, j] = -dr[-mij]
                else:
                    ds = Symbol(char + '%d%d'%(i,j), dummy=True)
                    d[ds] = mij
                    dr[mij] = ds
                    new_mat[i, j] = ds
    return new_mat, d

if __name__ == "__main__":
        import doctest
        doctest.testmod()
