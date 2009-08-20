from sympy import (Symbol, symbols, Basic, Function, Mul, Pow, Matrix, sin,
        cos, tan, cot, S, eye, Add, trigsimp, expand, pretty, Eq, collect, sqrt,
        sympify, factor, zeros, simplify, solve_linear_system, ratsimp,
        powsimp)
from sympy.printing.pretty.pretty import PrettyPrinter
from sympy.printing.str import StrPrinter
import time

e1 = Matrix([1, 0, 0])
e2 = Matrix([0, 1, 0])
e3 = Matrix([0, 0, 1])
e1n = Matrix([-1, 0, 0])
e2n = Matrix([0, -1, 0])
e3n = Matrix([0, 0, -1])
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

class Dyad(Basic):
    """General dyad expression.
    """
    def __init__(self, v):
        """ v should be an additive expression of the form:
        (sympy expression)*UnitVector*UnitVector

        The sympy expression is optional, but typically would be an inertia
        scalar.
        """
        self.dict = {}
        if v.is_Add:
            for term in v.args:
                if term.args[-1].is_Pow:
                    self.dict.update({term.args[-1]: term.coeff(term.args[-1])})
                else:
                    self.dict.update({term.coeff(term.args[0]): term.args[0]})
        elif v.is_Mul:
            self.dict.update({v.args[-2:]: v.args[:-2]})
        elif v.is_Pow:
            self.dict.update({term.args[-1]: term.coeff(term.args[-1])})
        elif v == 0:
            pass
        else:
            raise NotImplementedError()

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
                    scalar_part = coeff*other.dot(d_term.args[1])
                    vec_dict[d_term.args[1]] = (vec_dict.get(d_term.args[1], 0)
                            + scalar_part)
                elif d_term.is_Pow:
                    scalar_part = coeff*other.dot(d_term.args[0])
                    vec_dict[d_term.args[0]] = (vec_dict.get(d_term.args[0], 0)
                            + scalar_part)
                else:
                    raise NotImplementedError()

        return Vector(vec_dict)

    def __str__(self):
        return PyDyStrPrinter().doprint(self)
        #return pydy_str(self)

class Inertia(Dyad):
    """Inertia dyadic.
    """
    def __new__(cls, frame, scalars):
        """Specify frame, scalars as:
        frame - ReferenceFrame
        scalars - List or tuple of I11, I22, I33, I12, I23, I13 inertia scalars
        """
        I11, I22, I33, I12, I23, I13 = scalars
        return Dyad(I11*frame[1]**2 + I22*frame[2]**2 + I33*frame[3]**2 +
                I12*frame[1]*frame[2] + I12*frame[2]*frame[1] +
                I23*frame[2]*frame[3] + I23*frame[3]*frame[2] +
                I13*frame[1]*frame[3] + I13*frame[3]*frame[1])

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
            new[uv] = expand(trigsimp(new[uv]))
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
        frame = most_frequent_frame(self)
        m = 0
        s = set([])
        for k1, k2 in ((x,y) for x in self.dict for y in self.dict):
            if (k2, k1) in s:
                continue
            else:
                s.add((k1, k2))
        for k1, k2 in s:
            if k1 == k2:
                    m += expand(self.dict[k1]**2)
            else:
                    m += 2*expand(self.dict[k1]*self.dict[k2]*dot(k1, k2))

        # Try to factor things if possible
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

    def acc(self, point=None, frame=None):
        """Calculate the acceleration of a point.

        Used without arguments, vel() returns the velocity of self
        relative to the Newtonian Frame, taking into account whether points
        have been declared as fixed in a frame or not.

        Used with arguments, .vel(point, frame) returns the velocity of self
        relative to point, as view by an observer fixed in frame.
        """

        a = Vector(0)
        if point == frame == None:
            for p in self.point_list:
                a += p._arel
        elif isinstance(point, Point) and isinstance(frame, ReferenceFrame):
            # Get the point list from point to self
            point_list = point.get_point_list(self)
            for i, pa in enumerate(point_list[:-1]):
                pb = point_list[i+1]
                set_intersect = pa._fixedin & pb._fixedin
                # Case when the two points are not fixed in the same frame
                if len(set_intersect) == 0:
                    a += dt(pb.vel(pa), frame)
                # Case when the two points are fixed in the same frame
                elif len(set_intersect) == 1:
                    a += cross(set_intersect.pop().ang_vel(frame),
                            pb._vrel)
                else:
                    raise NotImplementedError('Somehow these two points are \
                        both fixed in 2 or more of the same frames')
        return a

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

                    w1 = C[0,2]*(C[0,1].diff(t)) + C[1,2]*(C[1,1].diff(t)) + \
                        C[2,2]*(C[2,1].diff(t))

                    w2 = C[1,0]*(C[1,2].diff(t)) + C[2,0]*(C[2,2].diff(t)) + \
                        C[0,0]*(C[0,2].diff(t))

                    w3 = C[2,1]*(C[2,0].diff(t)) + C[0,1]*(C[0,0].diff(t)) + \
                        C[1,1]*(C[1,0].diff(t))

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
            if frame == None: frame = self.NewtonianReferenceFrame
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



#class Particle:
#    def __init__(self, p, m):
#        self.mass = Mass
#        self.point = p

#class RigidBody(ReferenceFrame, Particle):
#    def __init__(self, s, matrix=None, frame=None, omega=None, m=0, I=0):
#        ReferenceFrame.__init__(self, s, Matrix=None, frame=None, omega=None)
#        Particle.__init__(self,s+'O', m)
#        self.inertia = I


#class InertialSystem():
#    def __init__(self, s):
#        self.name = s
#
#    def __getitem__(self, i):
#        return self.triad[i-1]


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

    def setkindiffs(self, expr_dict, dependent_speeds=None):
        """Recursivly apply kinematic differential equations to velocity and
        angular velocity expressions of every Point and ReferenceFrame,
        respectively.

        Also forms the acceleration and angular acceleration of each point.

        """
        self.kindiffs = expr_dict
        for eqn in self.kindiffs.values():
            sins = eqn.atoms(sin)
            coss = eqn.atoms(cos)
            if sins is not None:
                self.trig_func_set.update(sins)
                self.sin_func_set.update(sins)
            if coss is not None:
                self.trig_func_set.update(coss)
                self.cos_func_set.update(coss)

        if dependent_speeds is not None:
            #self.dependent_speeds = dependent_speeds
            #ind_speeds = list(set(self.u_list) - set(dependent_speeds))
            #self.independent_speeds = []
            #for u in self.u_list:
            #    if u in ind_speeds:
            #        self.independent_speeds.append(u)
            for eqn in dependent_speeds.values():
                sins = eqn.atoms(sin)
                coss = eqn.atoms(cos)
                if sins is not None:
                    self.trig_func_set.update(sins)
                    self.sin_func_set.update(sins)
                if coss is not None:
                    self.trig_func_set.update(coss)
                    self.cos_func_set.update(coss)
        else:
            self.independent_speeds = self.u_list

        for c in self.cos_func_set:
            if c**2 not in self.csqrd_dict:
                self.csqrd_dict[c**2] = 1 - sin(c.args[0])**2

        # Form partial velocity expressions
        self.recursive_partials(self)
        self.recursive_partials(self.O)

        # Form angular accelerations of ReferenceFrames and accelerations of
        # Points
        self.recursive_acc(self)
        self.recursive_acc(self.O)

    def setdyndiffs(self, eqns):
        """
        Sets the dynamic equations of motion.
        """
        self.dyndiffs = eqns

        for ud, rhs in eqns.items():
            sins = rhs.atoms(sin)
            coss = rhs.atoms(cos)
            if sins is not None:
                self.trig_func_set.update(sins)
                self.sin_func_set.update(sins)
            if coss is not None:
                self.trig_func_set.update(coss)
                self.cos_func_set.update(coss)

        for c in self.cos_func_set:
            if c**2 not in self.csqrd_dict:
                self.csqrd_dict[c**2] = 1 - sin(c.args[0])**2

    def recursive_acc(self, PorF):
        """Recursively form acceleration of Points and angular acceleration of
        ReferenceFrames.
        """
        if isinstance(PorF, Point):
            if PorF._fixedin == set([]):
                PorF._arel = PorF._vrel.dt(self)#.subs(self.kindiffs)
            elif len(PorF._fixedin) == 1:
                frame = list(PorF._fixedin)[0]
                PorF._arel = frame.ang_vel(self).cross(PorF._vrel) + \
                    frame.ang_acc(self).cross(PorF.rel(PorF.parentpoint))#.subs(self.kindiffs)
            else:
                frame_counter = {}
                for frame in PorF._fixedin:
                    frame_counter[frame] = len(frame.get_frames_list(self))
                closest = min([(frame_counter[x], x) for x in frame_counter])[1]
                PorF._arel = (closest.ang_vel(self).cross(PorF._vrel) +
                        closest.ang_acc(self).cross(PorF.rel(
                        PorF.parentpoint)))
        elif isinstance(PorF, ReferenceFrame):
            PorF._alpharel = PorF._wrel.subs(self.kindiffs).dt(PorF.NewtonianReferenceFrame)

        #  Initiate recursion
        if PorF.children == []:
            return
        else:
            for child in PorF.children:
                self.recursive_acc(child)

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
        """Recursively form the relative partial velocities of each point.
        """
        # Substitute into appropriate velocity/angular velocity
        if isinstance(PorF, Point):
            PorF._partialvrel = PorF._vrel.partials(self.independent_speeds)
        elif isinstance(PorF, ReferenceFrame):
            PorF._partialwrel = PorF._wrel.partials(self.u_list)
        else:
            raise NotImplementedError()

        #  Initiate recursion
        if PorF.children == []:
            return
        else:
            for child in PorF.children:
                self.recursive_partials(child)

    def declare_coords(self, string, number, list=True):
        """Declare the generalized coordinates and their time derivatives.
        """
        q_list, q_list, qdot_list = gcs(string, number, list)
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
        self.qdot_list_s = [Symbol(str(q.func)+'p') for q in q_list]
        self.q_list_dict = dict(zip(q_list, self.q_list_s))
        self.q_list_dict_back = dict(zip(self.q_list_s, q_list))
        self.qdot_list_dict = dict(zip(qdot_list, self.qdot_list_s))
        self.qdot_list_dict_back = dict(zip(self.qdot_list_s, qdot_list))
        # Update the comprehensive symbol dictionaries
        for d in (self.q_list_dict, self.qdot_list_dict):
            self.symbol_dict.update(d)
        for d in (self.q_list_dict_back, self.qdot_list_dict_back):
            self.symbol_dict_back.update(d)
        return q_list, q_list, qdot_list

    def declare_speeds(self, string, number, list=True):
        """Declare the generalized speeds and their time derivatives.
        """
        u_list, u_list, udot_list = gcs(string, number, list)
        self.u_list = u_list
        self.udot_list = udot_list

        # Generate lists of Symbol objects instead of Function objects
        self.u_list_s = [Symbol(str(u.func)) for u in u_list]
        self.udot_list_s = [Symbol(str(u.func)+'p') for u in u_list]

        # Generate a set of cross terms
        for ui in u_list:
            for uj in u_list:
                self.crossterms.update(set([ui*uj]))

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

        return u_list, u_list, udot_list

    def declare_parameters(self, string):
        """Declare the parameters (constants) of the system.
        """
        self.parameter_list = symbols(string)
        return self.parameter_list

    def impose_constraints(self, eqns, dependent, method='ADJ'):
        """Form the constraint matrix associated with linear velocity
        constraints.

        Must be of the form:
        B*u == 0
        Where B is (m, n)

        Includes both differentiated holonomic constraints and nonholonomic
        constraints.  Both types of constraints must be linear in the time
        derivatives of the generalized coordinates.

        Requires that both kinematic_chain and/or self.set_nhc_eqns has been
        called.
        """
        #  Check to make sure the user passes compatible equations and
        #  determine the columns associated with the independent and dependent
        #  speeds
        m = len(eqns)
        md = len(dependent)
        n = len(self.u_list)
        if m != md:
            raise ValueError('Number of equations must equal number of\
                dependent speeds.')
        independent_ci = []
        dependent_ci = []
        try:
            for ud in dependent:
                dependent_ci.append(self.u_list.index(ud))
            dependent_ci.sort()
            independent_ci = list(set(range(n)) - set(dependent_ci))
            independent_ci.sort()
        except ValueError:
            print('Each of the dependent speeds must be in the speed list used\
                in the declare_speeds.')
        # Assign the dependent and independent speeds
        self.dependent_speeds = [self.u_list[jd] for jd in dependent_ci]
        self.independent_speeds = [self.u_list[ji] for ji in independent_ci]

        # Put the equations in Matrix form and create a matrix with dummy
        # symbols representing non-zero entries
        B = zeros((m, n))
        B_dummy = zeros((m, n))
        d = {}
        for i in range(m):
            for j in range(n):
                bij = eqns[i].lhs.expand().coeff(self.u_list[j])
                if bij is not None:
                    dummy_sym = Symbol('a', dummy=True)
                    d[dummy_sym] = bij
                    B[i, j] = bij
                    B_dummy[i, j] = dummy_sym

        # Generate the independent and dependent matrices
        # i.e. Bd * ud + Bi * ui = 0
        # Where Bd is m x m and Bi is m x (n-m)
        Bd = zeros((m, m))
        Bi = zeros((m, n-m))
        for j, jd in enumerate(dependent_ci):
            Bd[:, j] = B_dummy[:, jd]
        for j, ji in enumerate(independent_ci):
            Bi[:, j] = B_dummy[:, ji]

        # Invert the Bd matrix and determine the dependent speeds
        # ud = -inv(Bd) * Bi * ui
        Bdinv = Bd.inv(method=method)
        T = - (Bdinv * Bi).subs(d)

        # Create dictionary for dependent speeds
        dep = {}
        for u in dependent:
            dep[u] = (T[i, :] * Matrix(self.independent_speeds))[0]
        self.dependent_speeds_eqs = dep
        return B, T, dep

    """
    def solve_constraint_matrix(self, dependent_qdots, method='GE'):
        Solve the constraint matrix for the dependent qdots in terms of the
        independent qdots.

        When solving for the dependent speeds in terms of a the independent
        ones, a symbolic matrix inversion is necessary.  The method parameter
        can be either 'GE', 'ADJ', or 'LU'.  See inv() for more information on
        the differences.

        d_column_index = []
        for qd in dependent_qdots:
            d_column_index.append(self.qdot_list.index(qd))
        i_column_index = list(set(range(0, len(self.qdot_list))) -
                set(d_column_index))

        rows, columns = self.constraint_matrix.shape
        if rows != len(d_column_index):
            raise ValueError('Number of dependent qdots should equal number of constraint equations')
        if (columns-rows) != len(i_column_index):
            raise ValueError('Number of independent qdots should equal number of qdots minus number of dependent qdots')
        ds_mat = zeros([rows, rows])
        is_mat = zeros([rows, columns-rows])
        for i, j in enumerate(d_column_index):
            ds_mat[:, i] = self.constraint_matrix[:,j]
        for i, j in enumerate(i_column_index):
            is_mat[:, i] = self.constraint_matrix[:,j]
        # Need to now invert ds_mat, and form -inv(ds_mat)*is_mat

        ds_mat_dummy = zeros([rows, rows])
        is_mat_dummy = zeros([rows, columns-rows])
        d = {}
        for i in range(rows):
            for j in range(rows):
                if ds_mat[i, j] != 0:
                    s = Symbol('a', dummy=True)
                    d[s] = ds_mat[i, j]
                    ds_mat_dummy[i, j] = s
            for j in range(columns-rows):
                if is_mat[i, j] != 0:
                    s = Symbol('a', dummy=True)
                    d[s] = is_mat[i, j]
                    is_mat_dummy[i, j] = s
        matr_dummy = -ds_mat_dummy.inv(method=method)*is_mat_dummy

        independent_speeds =  Matrix([[self.qdot_list[i]] for i in i_column_index])

        dependent_rates = {}
        for r, qd in enumerate(dependent_qdots):
            dependent_rates[qd] =\
                (matr_dummy[r,:]*independent_speeds)[0].subs(d)

        return dependent_rates
    """

    def frstar(self):
        """Computes the generalized inertia forces of the system.
        """
        self.recursive_frstar(self.O)
        self.recursive_frstar(self)

    def recursive_frstar(self, PorF):
        """Recursively computes generalized inertia forces for each particle
        and rigid body in the system.
        """
        if isinstance(PorF, Point):
            if not hasattr(PorF, 'partialv'): self._partialv(PorF)
            if PorF.mass == 0:
                PorF.gen_inertia_force = [0] * len(self.u_list)
            else:
                PorF.gen_inertia_force = [-PorF.mass * PorF.acc().dot(pv)
                        for pv in PorF.partialv]
        elif isinstance(PorF, ReferenceFrame):
            if not hasattr(PorF, 'partialw'): self._partialw(PorF)
            if PorF.inertia == Inertia(self, (0,0,0,0,0,0)):
                PorF.gen_inertia_force = [0] * len(self.independent_speeds)
            else:
                alph = PorF.ang_acc()
                I = PorF.inertia
                w = PorF.ang_vel()
                PorF.gen_inertia_force = [pw.dot(-alph.dot(I)
                    - w.cross(I.rdot(w))) for pw in PorF.partialw]
        else:
            raise NotImplementedError()

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
            if not hasattr(PorF, 'partialv'): self._partialv(PorF)
            if PorF.force == Vector(0):
                PorF.gen_active_force = [0] * len(self.independent_speeds)
            else:
                PorF.gen_active_force = [PorF.force.dot(pv) for pv in
                    PorF.partialv]
        elif isinstance(PorF, ReferenceFrame):
            if not hasattr(PorF, 'partialw'): self._partialw(PorF)
            if PorF.torque == Vector(0):
                PorF.gen_active_force = [0] * len(self.independent_speeds)
            else:
                PorF.gen_active_force = [PorF.torque.dot(pw) for pw in
                    PorF.partialw]
        else:
            raise NotImplementedError()

        #  Initiate recursion
        if PorF.children == []:
            return
        else:
            for child in PorF.children:
                self.recursive_fr(child)

    def _partialv(self, point):
        """Computes the r absolute partial velocities of a point
        """
        if point.parentpoint == None:
            point.partialv = [Vector(0)] * len(self.independent_speeds)
        else:
            point.partialv = [Vector(point.parentpoint.partialv[r] +\
                    point._partialvrel[r]) for r in \
                    range(len(self.independent_speeds))]

    def _partialw(self, frame):
        """Computes the r absolute partial velocities of a ReferenceFrame
        """
        if frame.parentframe == None:
            frame.partialw = [Vector(0)] * len(self.independent_speeds)
        else:
            frame.partialw = [Vector(frame.parentframe.partialw[r] +
                    frame._partialwrel[r]) for r in \
                    range(len(self.independent_speeds))]

    def set_nhc_eqns(self, *args):
        """Assigns nonholonomic constraint equations, forms constraint matrix.
        """
        for arg in args:
            for eqn in arg:
                self.nhc_eqns.append(eqn)
        self.form_constraint_matrix()

    def gravity(self, v):
        """Applies a gravitational force to each particle and rigid body in the
        system.
        """
        v = Vector(v)
        self.recursive_gravity(self.O, v)

    def recursive_gravity(self, Point, v):
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
        self.fr()
        self.frstar()

        self.kanes_equations = [0] * len(self.independent_speeds)
        self.recursive_eoms(self.O)
        self.recursive_eoms(self)
        return self.kanes_equations

    def recursive_eoms(self, PorF):
        for r in range(len(self.independent_speeds)):
            self.kanes_equations[r] += PorF.gen_active_force[r] +\
                PorF.gen_inertia_force[r]

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
        udsort = []
        usort = []
        for i, ud in enumerate(self.udot_list):
            if ud in self.dyndiffs:
                udsort.append(ud)
                usort.append(self.u_list[i])

        for q in self.q_list:
            s += str(q) + ', '
        for u in usort:
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

        if hasattr(self, 'dependent_speeds'):
            ode_func_string += '    # Dependent generalized speeds\n'
            for u in self.dependent_speeds:
                ode_func_string += '    ' + str(u) + ' = ' +\
                    str(self.dependent_speeds_eqs[u].subs(self.qdot_list_dict)) + '\n'

        ode_func_string += '    # Kinematic differential equations\n'
        qdl = []
        for qd in self.qdot_list:
            if qd in self.kindiffs: qdl.append(qd)

        for qd in qdl:
            ode_func_string += '    ' + str(qd)[:-1] + 'p' + ' = ' + str(self.kindiffs[qd]) + '\n'
            dxdt_list += str(qd)[:-1] + 'p, '

        ode_func_string += '    # Dynamic differential equations\n'

        for ud in udsort:
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
            a += '    ' + s[:-2] + ' = ' + 'parameter_list\n'
            a += '    # Unpacking the coordinates\n'
            s = ""
            for q in self.q_list:
                s += str(q) + ', '
            a += '    ' + s[:-2] + ' = ' + 'q\n'

            trig_func_set = set([])
            a_temp = ""
            ret_string = ""
            for arg in args:
                ret_string += "["
                if isinstance(arg, tuple) and len(arg) == 2:
                    if isinstance(arg[0], Point) and isinstance(arg[1], Point):
                        pos = [arg[0].rel(arg[1]).dot(self[i]) for i in (1,2,3)]
                        for i, p in enumerate(pos):
                            nv = "p_" + arg[1].name + "_" + arg[0].name +\
                                "_%d"%(i+1)
                            a_temp += "    " + nv + " = " + str(p) + "\n"
                            ret_string += nv + ", "
                            trig_terms = p.atoms(sin, cos, tan)
                            if trig_terms:
                                trig_func_set.update(trig_terms)
                        ret_string = ret_string[:-2] + "], "
                    elif isinstance(arg[0], (UnitVector, Vector)):
                        if isinstance(arg[0], (UnitVector, Vector)):
                            axis = [arg[0].dot(self[i]) for i in (1,2,3)]
                        else:
                            raise ValueError('Axis must be a Vector or UnitVector')
                        # XXX TODO allow for a general sympy expression for the
                        # angle
                        if arg[1] in self.q_list or arg[1] == 0:
                            angle = arg[1] or S.Zero 
                            trig_terms = angle.atoms(sin, cos, tan)
                            if trig_terms:
                                trig_func_set.update(trig_terms)
                        else:
                            raise ValueError('Angle must be in the coordinate\
                                    list')
                        for j, p in enumerate(axis):
                            nv = arg[0].frame.name + str(arg[0].i) + "_%d"%(j+1)
                            a_temp += "    " + nv + " = " + str(p) + "\n"
                            ret_string += nv + ", "
                            trig_terms = p.atoms(sin, cos, tan)
                            if trig_terms:
                                trig_func_set.update(trig_terms)
                        ret_string += str(angle) +"], "
                else:
                    raise TypeError('Optional parameters must be 2-tuples with\
                        either two Point objects or an Axis and an Angle')

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

    def form_kindiffs(self, eqns, qdot_list, method='ADJ'):
        """Given a list of equations and linear terms, form the tranformation
        matrix.
        """
        m = len(eqns)
        n = len(qdot_list)
        if m != n:
            raise ValueError('Number of equations must equal number of unkowns.')
        M = zeros((m,n))
        for i in range(m):
            for j in range(n):
                mij = eqns[i].rhs.coeff(qdot_list[j])
                M[i,j] = mij if mij is not None else S(0)
        self.transform_matrix = M
        d = {}
        dummy = zeros((m,n))
        for i in range(m):
            for j in range(n):
                if M[i,j] != 0:
                    s = Symbol('a%d%d'%(i,j), dummy=True)
                    d[s] = M[i,j]
                    dummy[i,j] = s
        Minv = dummy.inv(method=method)
        for i in range(n):
            for j in range(n):
                if Minv[i,j] != 0:
                    num,den = Minv[i,j].as_numer_denom()
                    Mij = simplify(num.expand() / den.expand()).subs(d)
                    num,den = Mij.as_numer_denom()
                    Minv[i, j] = \
                        (num.expand().subs(self.csqrd_dict).expand() / \
                        den.expand().subs(self.csqrd_dict).expand())#.subs(self.tan_dict)
                    sins = Minv[i,j].atoms(sin)
                    coss = Minv[i,j].atoms(cos)
                    tans = Minv[i,j].atoms(tan)
                    if sins:
                        self.trig_func_set.update(sins)
                        self.sin_func_set.update(sins)
                    if coss:
                        self.trig_func_set.update(coss)
                        self.cos_func_set.update(coss)
                    if tans:
                        self.trig_func_set.update(tans)
                        self.tan_func_set.update(tans)

        self.transform_matrix_inv = Minv
        kindiffs = {}
        rhs = Matrix([eqn.lhs for eqn in eqns])
        for i, qd in enumerate(qdot_list):
            kindiffs[qd] = (Minv[i,:]*rhs)[0]
        return M, Minv, kindiffs

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

#    """
#    def form_kindiffs(self, T, qdot_list, u_list):
#        Invert the transformation u = T*q' and determine the kinematic
#        differential equations.
#
#        qdot_list and u_list must be of the same length.  For holonomic systems
#       with no unreduced coordinates, this will always be the case.  For 
#       systems with unreduced coordinates and/or nonholonomic constraints,
#       there will generally be fewer u's than qdots, so one must choose the u's
#       in a way such that they only depend on the same number of qdots.
#
#        m, n = T.shape
#        system = T.col_insert(n, Matrix(u_list))
#        d = {}
#        system2 = zeros(system.shape)
#        for i in range(system.shape[0]):
#            for j in range(system.shape[1]):
#                if system[i, j] != 0:
#                    s = Symbol("a", dummy=True)
#                    d[s] = system[i, j]
#                    system2[i, j] = s
#        kindiffs = solve_linear_system(system2, *qdot_list)
#
#        for qd in qdot_list:
#            kindiffs[qd] = trigsimp(kindiffs[qd].subs(d).expand())
#
#        return kindiffs
#    """

def most_frequent_frame(vector):
    """Determines the most frequent frame of all unitvector terms in a vector.
    """
    frame_counter = {}
    for uv in vector.dict:
        frame_counter[uv.frame] = frame_counter.get(uv.frame, 0) + 1
    return max([(frame_counter[x], x) for x in frame_counter])[1]

def kinematic_chain(point1, point2, r=None, vec_list=None):
    """Close a kinematic loop and form the associated constraint equations.

    point1: begin point
    point2: end point
    r: vector from point2 to point1 which closes the loop kinematic loop
    vec_list:  List of UnitVectors to dot the constraint with

    Returns the kinematic constraint equations, along with the time derivatives
    of those equations.
    """
    if r is None:
        loop = point2.rel(point1)
    else:
        loop = point2.rel(point1) + r

    if vec_list is None:
        frame = most_frequent_frame(loop)
        vec_list = [frame[i] for i in (1, 2, 3)]
    hc_eqs = []
    dhc_eqs = []

    q_list = point1.NewtonianFrame.q_list
    q_list_s = point1.NewtonianFrame.q_list_s
    q_list_dict = point1.NewtonianFrame.q_list_dict
    q_list_dict_back = point1.NewtonianFrame.q_list_dict_back

    qdot_list = point1.NewtonianFrame.qdot_list
    qdot_list_s = point1.NewtonianFrame.qdot_list_s
    qdot_list_dict = point1.NewtonianFrame.qdot_list_dict
    qdot_list_dict_back = point1.NewtonianFrame.qdot_list_dict_back

    cm_row_list = []
    # Generate the scalar holonomic constraint equations
    for uv in vec_list:
        hc = dot(uv, loop)
        if hc != 0:
            hc_eqs.append(hc)
            hc_s = hc.subs(q_list_dict)
            # Get a list of the coordinates involved in the constraint equation
            coords_s = list(hc_s.atoms(Symbol) & set(q_list_s))
            print coords_s
            coords_gc = [cd_s.subs(q_list_dict_back) for
                    cd_s in coords_s]
            coords_gc_dot = [c_gc.diff(t) for c_gc in coords_gc]

            dhc = S(0)
            cm_row_i = [0]*len(qdot_list)
            for j, q in enumerate(coords_s):
                col_index =\
                    point1.NewtonianFrame.qdot_list.index(coords_gc_dot[j])
                cm_entry = hc_s.diff(q).subs(q_list_dict_back)
                cm_row_i[col_index] = cm_entry
                dhc += cm_entry*coords_gc_dot[j]
            if any(cm_row_i):
                cm_row_list.insert(-1, cm_row_i)

    point1.NewtonianFrame.constraint_matrix = Matrix(cm_row_list)

    subs_dhc_eqs = [collect(eq.subs(qdot_list_dict), qdot_list_s) for eq in dhc_eqs]
    if len(dhc_eqs) != 0:
        dhc_eqs = [eq.subs(qdot_list_dict_back) for eq in subs_dhc_eqs]
    # Append the constraints to the Newtonian Reference constraint list
    for hc_eqn in hc_eqs:
        point1.NewtonianFrame.hc_eqns.append(hc_eqn)

    for dhc_eqn in dhc_eqs:
        point1.NewtonianFrame.dhc_eqns.append(dhc_eqn)

    # Return the constraint equations so the user can look at them if they
    # want.
    return hc_eqs, subs_dhc_eqs

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

def InertiaForce(p):
    """Computes Inertia force of a point, defined as:

    R^* = -m*a

    Equation 4.11.3 or 4.11.6 from Dynamics: Theory and Application
    """
    if isinstance(p, Point):
        if p.mass == 0:
            return Vector(0)
        else:
            IF_dict = dict([(uv, -p.mass*coef) for uv, coef in \
                    p.acc().dict.items()])
            return Vector(IF_dict)
    else:
        raise TypeError('InertiForce requires a Point as its argument')

def InertiaTorque(frame):
    """Computes Inertia torque of a ReferenceFrame/Rigid Body, defined as:

    T^* = - cross(alpha, I) - cross(omega, dot(I, omega))

    Where I is the central inertia dyadic of the rigid body, omega is the
    angular velocity of the body, and alpha is the angular acceleration of the
    body.

    Equation 4.11.8 from Dynamics: Theory and Application
    """
    if isinstance(frame, ReferenceFrame):
            return dot(-frame.ang_acc(), frame.inertia)\
                    - cross(frame.ang_vel(),\
                    dot(frame.inertia, \
                        frame.ang_vel()))
    else:
        raise NotImplementedError("frame must be a Reference Frame")


def GeneralizedCoordinate(s, constant=False):
    gc = Symbol(s)(Symbol('t'))
    gc.is_gc = True
    if constant==True:
        gc.fdiff = lambda argindex: 0
    gc.__repr__ = lambda self: PyDyStrPrinter().doprint(self)
    gc.__str__ = lambda self: PyDyStrPrinter().doprint(self)
    return gc

def gcs(s, number=1, list=False):
    gc_list = [GeneralizedCoordinate(s[0]+str(i)) for i in range(1, number + 1)]
    if list == False:
        if number == 1:
            return gc_list[0]
        else:
            return gc_list
    elif list == True:
        gcd_list = [gc.diff(t) for gc in gc_list]
        return (gc_list, gc_list, gcd_list)

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
        return str(expr.args[0].func) + "'"*len(expr.args[1:])

    def _print_Matrix(self, expr):
        return expr._format_str(lambda elem: elem.__str__())

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

if __name__ == "__main__":
        import doctest
        doctest.testmod()
