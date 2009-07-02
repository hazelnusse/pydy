from sympy import Symbol, Basic, Function, Mul, Pow, Matrix, sin, \
        cos, S, eye, Add, trigsimp, expand
from sympy.printing.pretty.pretty import PrettyPrinter
from sympy.printing.str import StrPrinter

e1 = Matrix([1, 0, 0])
e2 = Matrix([0, 1, 0])
e3 = Matrix([0, 0, 1])
e1n = Matrix([-1, 0, 0])
e2n = Matrix([0, -1, 0])
e3n = Matrix([0, 0, -1])
zero = Matrix([0, 0, 0])
t = Symbol("t")



class UnitVector(Basic):
    """A standard unit vector  with a symbolic and a numeric representation"""

    # XXX: UnitVector should be noncommutative in Mul, but currently it is
    # probably commutative. However, we haven't found any case, where this
    # actually fails, so until we find some, let's leave it as is and keep this
    # in mind.
    def __init__(self, frame, i=0): #=-1,num=None):
        self.frame = frame    # Parent reference frame
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
        #elif i == -1:
        #    self.v['sym'] = Symbol('-'+s.lower()+str(abs(i)))
        #    self.v['num'] = e1n
        #elif i == -2:
        #    self.v['sym'] = Symbol('-'+s.lower()+str(abs(i)))
        #    self.v['num'] = e2n
        #elif i == -3:
        #    self.v['sym'] = Symbol('-'+s.lower()+str(abs(i)))
        #    self.v['num'] = e3n
        elif i == 0:
            self.v['sym'] = Symbol(s.lower()+str(0))
            self.v['num'] = zero

    def __str__(self):
        return pydy_str(self)

    def __repr__(self):
        return pydy_pretty(self)

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
        elif isinstance(other, Mul) or isinstance(other, Add):
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
        """Expresses a UnitVector with UnitVectors fixed in the specified frame.
        """

        if self.frame == frame:
            return self
        else:
            matrices = self.frame.get_rot_matrices(frame)
            if len(matrices) == 1 and matrices[0]*self.v['num'] == self.v['num']:
                return frame[self.i]
            else:
                u = self.v['num']
                for m in reversed(matrices):
                    u = m*u
                #u[0] = trigsimp(u[0])
                #u[1] = trigsimp(u[1])
                #u[2] = trigsimp(u[2])
                return Vector(u[0]*frame[1] + u[1]*frame[2] + u[2]*frame[3])

    def dot(self, other):
        if isinstance(other, UnitVector):
            c = other.express(self.frame)
            if isinstance(c, UnitVector):
                return (self.v['num'].T * c.v['num'])[0]
            elif isinstance(c, Vector):
                s = S(0)
                for k, coef in c.dict.items():
                    if self == k:
                        s += coef
                return trigsimp(s)
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
            c = other.express(self.frame)
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
                if len(cp) == 1:
                    if cp.values()[0] == 1:
                        return cp.keys()[0]  # Return a UnitVector object
                    else:
                        return Vector(cp)
                else:
                    cp1 = Vector(cp)
                    cp2 = cp1.express(other.frame)
                    if isinstance(cp2, UnitVector):
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
                elif len(cp1.dict) <= len(cp2.dict):
                    return cp1
                else:
                    return cp2
        elif isinstance(other, Vector):
            return cross_with_Vector(self, other)
        else:
            raise NotImplementedError()

    def dt(self, diff_frame):
        if isinstance(diff_frame, ReferenceFrame):
            return cross(self.frame.get_omega(diff_frame), self)
        else:
            raise NotImplementedError()

class Dyad(Basic):
    """
    General dyad expression
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
                    #if d_term.args[1] in vec_dict:
                    #    vec_dict[d_term.args[1]] += scalar_part
                    #else:
                    #    vec_dict.update({d_term.args[1]: scalar_part})
                elif d_term.is_Pow:
                    scalar_part = coeff*other.dot(d_term.args[0])
                    vec_dict[d_term.args[0]] = (vec_dict.get(d_term.args[0], 0)
                            + scalar_part)
                    #if d_term.args[0] in vec_dict:
                    #    vec_dict[d_term.args[0]] += scalar_part
                    #else:
                    #    vec_dict.update({d_term.args[0]: scalar_part})
                else:
                    raise NotImplementedError()

        return Vector(vec_dict)

    def __str__(self):
        return pydy_str(self)

class Inertia(Dyad):
    """Inertia dyadic.

    """
    def __new__(cls, frame, scalars):
        """Specify frame, scale as:
        frame - ReferenceFrame
        scalars - List or tuple of I11, I22, I33, I12, I23, I13 inertia scalars
        """
        I11, I22, I33, I12, I23, I13 = scalars
        return Dyad(I11*frame[1]**2 + I22*frame[2]**2 + I33*frame[3]**2 +
                I12*frame[1]*frame[2] + I12*frame[2]*frame[1] +
                I23*frame[2]*frame[3] + I23*frame[3]*frame[2] +
                I13*frame[1]*frame[3] + I13*frame[3]*frame[1])

class Vector(Basic):
    """
    General vector expression.  Internally represented as a dictionary whose
    keys are UnitVectors and whose key values are the corresponding coefficient
    of that UnitVector.  For example:
    N = ReferenceFrame("N")
    x, y, z = symbols('x y z')
    v = Vector(x*N[1]+y*N[2] + z*N[3])
    then v would be represented internally as:
    {N[1]: x, N[2]: y, N[3]: z}
    """

    def __init__(self, v):
        """
        Initialize a Vector in two ways:
        Method 1 (dictionary way):
        v = Vector({UnitVectors : Coefficients})
        for example:
        v = Vector({A[1] : sin(q1)})
        or
        Method 2 (sympy expression way):
        v = Vector(sin(q1)*A[1])

        Vector objects are internally represented as dictionaries whose keys
        are the UnitVectors and whose values are the coefficients of those
        UnitVectors.  In the above example, this would imply:
        v.dict == {A[1] : sin(q1)}
        """

        if isinstance(v, dict):
            for k in v.keys():
                if v[k] == 0:  v.pop(k)
            self.dict = v
        elif isinstance(v, Vector):
            self.dict = v.dict
        else:
            vdict = self.parse_terms(v)
            for k in vdict.keys():
                if vdict[k] == 0:  vdict.pop(k)
            self.dict = vdict

    def __str__(self):
        return pydy_str(self)

    def __repr__(self):
        return pydy_str(self)

    def __add__(self, other):
        """Adds two Vector objects and return a new Vector object.
        v1 + v2   <----> v1.__add__(v2)
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
        """Subtracts two Vector objects and return a new Vector object.
        v1 - v2   <----> v1.__sub__(v2)
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
                #print 'self.dict = ',self.dict
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

    def __mul__(self, other):
        if isinstance(other, Dyad):
            return NotImplemented
        else:
            return Basic.__mul__(self, other)

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
        if isinstance(other, Vector):
            s = S(0)
            for k in self.dict:
                s += sum([self.dict[k]*other.dict[ko]*k.dot(ko) for ko in
                    other.dict])
            return s
        elif isinstance(other, UnitVector):
            return sum([self.dict[k]*k.dot(other) for k in self.dict])
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
                t2 = k.frame.get_omega(frame).cross(k)
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

    def mag(self):
        m = 0
        for k in self.dict:
            m += self.dict[k]**S(2)
        return m**(S(1)/S(2))

    def mag_sqr(self):
        m = 0
        for k in self.dict:
            m += self.dict[k]**S(2.0)
        return m

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
                for k in add_term_dict.keys():
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

class Point:
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

    def __init__(self, s, r=None, point=None, mass=None, force=None, frame=None):
        # When instantiated by ReferenceFrame
        if not any([r, point, mass, force]):
            self.name = s
            self.point_list = [self]
            self.parent = None
            self.pos = {self: Vector(0)}
            self.vel = {self: Vector(0)}
            self.NewtonianFrame = frame
        # When instantiated by locate method
        elif all([s, r, point]) and not any([mass, force]):
            #print 'locate method instantiation, point.point_list = ', \
            #    point.point_list
            r = Vector(r)
            self.name = s
            self.pos = {point: -r}
            point.pos[self] = r
            self.point_list = [self] + point.point_list
            self.parent = point
            self.NewtonianFrame = point.NewtonianFrame
        else:
            raise NotImplementedError()

    def locate(self, s, r, frame=None):
        r = Vector(r)
        if frame == None:
            newpoint = Point(s, r, self)
            newpoint.vel =  {newpoint.NewtonianFrame:
                    r.dt(newpoint.NewtonianFrame)}
        elif isinstance(frame, ReferenceFrame):
            newpoint = Point(s, r, self)
            newpoint.vel = {newpoint.NewtonianFrame:
                    cross(frame.get_omega(newpoint.NewtonianFrame), r)}
        else:
            raise NotImplementedError()
        return newpoint

    def get_point_list(self, other=None):
        """
        Gets the list of Points between Point self and Point other, including both.
        """
        if other == None:
            return self.point_list
        elif self == other:
            return [self]
        else:
            r2t = [f for f in reversed(other.point_list)]

            if len(self.point_list) == 1:
                return r2t
            elif len(r2t) == 1:
                return self.point_list

            r1t = [f for f in reversed(self.point_list)]
            i = 1

            while r1t[i] == r2t[i]:
                del r1t[i-1]
                del r2t[i-1]
                if len(r1t)<2 or len(r2t)<2:
                    break

            r1t.reverse()
            return r1t[:-1] + r2t


            """
            self_list = [self]
            other_list = [other]
            #Build the list from self to top of the tree.
            while self_list[-1].parent is not None:
                self_list.append(self_list[-1].parent)
            #Build the list from top of the tree to other.
            while other_list[0].parent is not None:
                other_list.insert(0, other_list[0].parent)
            print 'self_list: ',self_list
            print 'other_list: ', other_list
            if len(self_list) == 1:
                return other_list
            elif len(other_list) == 1:
                return self_list
            else:
                i = 0
                while self_list[-(i+2)] == other_list[i+1]:
                    del self_list[-(i+1)]
                    del other_list[i]
                    if len(self_list)<=2 or len(other_list)<=2:
                        return self_list[:-1] + other_list[1:]
                    i += 1
                return self_list[:-1] + other_list
            """




    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

class ReferenceFrame:
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
        if not any([frame, matrix, omega]):
            self.ref_frame_list = [self]
            self.O = Point(s + 'O', frame=self)
            self.point_list = [self.O]
        else:
            self.ref_frame_list = [self] + frame.ref_frame_list[:]

        self.name = s
        self.triad = [UnitVector(self, i) for i in (1,2,3)]
        self.transforms = {}
        self.parent = frame
        self.W = {}
        if omega != None:
            self.set_omega(omega, self.parent)
            frame.set_omega(-omega, self, force=True)

        if frame is not None:
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

    def rotate(self, name, axis, angle):
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
                omega = Vector(angle.diff(t)*self[axis])
            elif axis in set((-1, -2, -3)):
                matrix = self._rot(-axis, -angle)
                omega = Vector(-angle.diff(t)*self[-axis])
            elif isinstance(axis, (UnitVector, Vector)):
                raise NotImplementedError("Axis angle rotations not \
                    implemented.")
            else:
                raise ValueError("Invalid axis")
            return ReferenceFrame(name, matrix, self, omega)
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

        Example:

        N - A - D - E - F
            |
            B
            |
            C

        Then:

        C.get_frames_list(F) == [C, B, A, D, E, F]
        F.get_frames_list(C) == [F, E, D, A, B, C]
        """
        if self == frame:
            return [self]
        else:
            r2t = [f for f in reversed(frame.ref_frame_list)]

            if len(self.ref_frame_list) == 1:
                return r2t
            elif len(r2t) == 1:
                return self.ref_frame_list

            r1t = [f for f in reversed(self.ref_frame_list)]
            i = 1
            while r1t[i] == r2t[i]:
                del r1t[i-1]
                del r2t[i-1]
                if len(r1t)<2 or len(r2t)<2:
                    break

            r1t.reverse()
            return r1t[:-1] + r2t

    def get_rot_matrices(self, frame):
        """
        Returns a list of matrices to get from self to frame.
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
        """Sets the angular velocity of self relative to frame.
        """
        if self.W == {} or force:
            self.W[frame] = omega
        else:
            raise ValueError("set_omega has already been called.")

    def get_omega(self, frame):
        """
        Returns the angular velocity of self relative to "frame".
        """

        if frame in self.W:
            return self.W[frame]
        else:
            om = {}
            for term in self.get_omega_list(frame):
                for k in term.dict:
                    if k in om:
                        om[k] += term.dict[k]
                    else:
                        om.update({k: term.dict[k]})
            self.W.update({frame: Vector(om)})
            return self.W[frame]
            #return sum(self.get_omega_list(frame))

    def get_omega_list(self, frame):
        """
        Returns a list of simple angular velocities from self to frame.
        """
        frames = self.get_frames_list(frame)
        if frames == [self]:
            #return [S(0)]
            return [Vector({})]
        result = []
        for i, f in enumerate(frames[:-1]):
            result.append(f.W[frames[i+1]])
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

def express(v, frame):
    if (isinstance(v, UnitVector) or isinstance(v, Vector)) and \
            (isinstance(frame, ReferenceFrame)):
        return v.express(frame)
    else:
        raise NotImplementedError()

def dot(v1, v2):
    """Dot product between UnitVector, Vector, and Dyad classes

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
    if (isinstance(v1, UnitVector) or isinstance(v1, Vector)) and \
            (isinstance(v2, UnitVector) or isinstance(v2, Vector)):
                return v1.cross(v2)
    else:
        if not (isinstance(v1, UnitVector) or isinstance(v1, Vector)):
            v1 = Vector(v1)
        if not (isinstance(v2, UnitVector) or isinstance(v2, Vector)):
            v2 = Vector(v2)
        return v1.cross(v2)

def identify(a):
    """
    Takes a Mul instance and parses it as

    a = c * UnitVector() * UnitVector()

    and returns c, v1, v2, where v1 and v2 are the UnitVectors.
    """
    if isinstance(a, Mul):
        unit_vectors = []
        for b in a.args:
            if isinstance(b, UnitVector):
                unit_vectors.append(b)
            if isinstance(b, Pow):
                if isinstance(b.args[0], UnitVector):
                    unit_vectors.append(b.args[0])
                    unit_vectors.append(b.args[0])
        if len(unit_vectors) == 2:
            v1 = unit_vectors[0]
            v2 = unit_vectors[1]
            c = coeff(a, v1*v2)
            #XXX here is a bug, if a=B[1]*Derivative()*B[1] and we do coeff for
            #B[1]**2
            #print "identify, coeff", a, v1*v2, c
            return c, v1, v2

    return a, None, None

def identify_v1(a):
    """
    Takes a Mul instance and parses it as

    a = c * UnitVector()

    and returns c, v1 where v1 is the UnitVector.
    """
    if isinstance(a, UnitVector):
        return S(1), a
    elif isinstance(a, Mul):
        unit_vectors = []
        for b in a.args:
            if isinstance(b, UnitVector):
                unit_vectors.append(b)
            if isinstance(b, Pow):
                if isinstance(b.args[0], UnitVector):
                    unit_vectors.append(b.args[0])
                    unit_vectors.append(b.args[0])
        if len(unit_vectors) == 1:
            v1 = unit_vectors[0]
            c = a.coeff(v1)
            return c, v1

    return a, None

def cross_vectors(u, v):
    c1 = u[1]*v[2] - u[2]*v[1]
    c2 = -(u[0]*v[2] - u[2]*v[0])
    c3 = u[0]*v[1] - u[1]*v[0]
    return c1, c2, c3

def dot_vectors(u, v):
    return u[0]*v[0]+u[1]*v[1]+u[2]*v[2]

def coeff(e, x):
    """
    Workaround the bug in sympy.
    """
    if isinstance(x, list):
        r = []
        for xi in x:
            ri = e.coeff(xi)
            if ri is None:
                r.append(S(0))
            else:
                r.append(ri)
        return r
    else:
        r = e.coeff(x)
        if r is None:
            return S(0)
        else:
            return r

def coeffv(v, scalar):
    if isinstance(v, Vector):
        return v.coeffv(scalar)
    else:
        raise NotImplementedError()


def expression2vector(e, frame):
    """
    Converts a sympy expression "e" to a coefficients vector in the frame "frame".
    """
    u1 = coeff(e, frame[1])
    u2 = coeff(e, frame[2])
    u3 = coeff(e, frame[3])
    return u1, u2, u3

def vector2expression(u, frame):
    """
    Converts a coefficients vector to a sympy expression in the frame "frame".
    """
    return u[0]*frame[1] + u[1]*frame[2] + u[2]*frame[3]

def dt(v, frame):
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

def InertiaForce(m, a):
    """Computes Inertia force, defined as:

    R^* = -m*a

    Equation 4.11.3 or 4.11.6 from Dynamics: Theory and Application
    """
    if isinstance(a, Vector):
        IF_dict = dict([(uv, -m*coef) for uv, coef in a.dict.items()])
    elif isinstance(a, UnitVector):
        IF_dict = {a: -m}
    else:
        raise TypeError('The acceleration a must be a Vector or UnitVector, \
            instead a %s object was given.' % str(type(a)))
    return Vector(IF_dict)

def InertiaTorque(I, omega, alpha):
    """Computes Inertia torque, defined as:

    T^* = - cross(alpha, I) - cross(omega, dot(I, omega))

    Where I is the central inertia dyadic of the rigid body, omega is the
    angular velocity of the body, and alpha is the angular acceleration of the
    body.

    Equation 4.11.8 from Dynamics: Theory and Application
    """
    if isinstance(I, Dyad):
        if isinstance(alpha, (UnitVector, Vector)) and isinstance(omega,
                (UnitVector, Vector)):
            return dot(-alpha, I) - cross(omega, dot(I, omega))
        else:
            raise NotImplementedError("Alpha and Omega must be UnitVector or \
                Vector Instances")
    else:
        raise NotImplementedError("I must be a Dyad")

def RigidBody(point, frame, mass, inertia):
    """
    Assign mass and inertia to a Point and a Reference frame, effectively
    creating a Rigid Body.

    """
    point.mass = mass
    frame.inertia = inertia

def Particle(point, mass):
    """
    Assign mass to a Point, effectively creating a particle
    """
#class GeneralizedCoordinate(Symbol):
    #def __init__(self, name, depends_on=Symbol('t'), *args):
    #    self.dv = depends_on

    #def dt(self):
    #    return self.diff(Symbol('t'))

    #def diff(self, var):
    #    if var == Symbol('t'):
    #        return GeneralizedCoordinate(self.name + "'")
    #    else:
    #        return Derivative(self, var)

    #def _eval_derivative(self, s):
    #    if s == self:
    #        return S.One
    #    elif s == self.dv:
    #        return GeneralizedCoordinate(self.name + "'")
    #    else:
    #        return S.Zero

def gcs(s, number=1):
    gc_list = [Function(s + str(i))(Symbol('t')) for i in
        range(number + 1)[1:]]

    for i, gc in enumerate(gc_list):
        gc_list[i].__str__ = lambda x: str(x.func)
    return gc_list

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
                    # All other terms pad with spaces and add parenthesis
                    else:
                        if isinstance(e.dict[k], Add):
                            s += (' + (' + self.doprint(e.dict[k])
                                + ')' + small_dot + self.doprint(k))
                        elif isinstance(e.dict[k], (Mul, Pow)):
                            coef = (self.doprint(e.dict[k]) + small_dot +
                                    self.doprint(k))
                            if coef[0] == '-':
                                s += ' - ' + coef[1:]
                            else:
                                s += ' + ' + coef
                        else:
                            s += (' + ' + self.doprint(e.dict[k]) + small_dot +
                                    self.doprint(k))

            return s
        else:
            return "0>"

    def _print_Function(self, e):
        return str(e.func)

    #def _print_Derivative(self, e):
    #    return "%s'" % str(e.args[0].func)

    def _print_Derivative(self, expr):
        return "%s'" % str(expr.args[0].func)


class PyDyPrettyPrinter(PrettyPrinter):
    def _print_UnitVector(self, e):
        class Fake(object):
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
                        s += sign + self.doprint(k)
                        i += 1
                    # All other terms put the sign and pad with spaces
                    else:
                        if e.dict[k] == 1: sign = '+'
                        if e.dict[k] == -1: sign = '-'
                        s += ' ' + sign + ' ' + self.doprint(k)
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
                    # All other terms pad with spaces and add parenthesis
                    else:
                        if isinstance(e.dict[k], Add):
                            s += (' + (' + self.doprint(e.dict[k])
                                + ')' + small_dot + self.doprint(k))
                        elif isinstance(e.dict[k], (Mul, Pow)):
                            coef = (self.doprint(e.dict[k]) + small_dot +
                                    self.doprint(k))
                            if coef[0] == '-':
                                s += ' - ' + coef[1:]
                            else:
                                s += ' + ' + coef
                        else:
                            s += (' + ' + self.doprint(e.dict[k]) + small_dot +
                                    self.doprint(k))

            return s
        else:
            return "\033[1m" + "0" + "\033[0;0m"


    def _print_Function(self, e):
        return "%s" % str(e.func)

    def _print_Derivative(self, e):
        return "%s'" % str(e.args[0].func)

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
        name = str(e.args[0])
        if name[0] == "q":
            index = name[1]
            return "s%s" % index
        else:
            return e

    def _print_cos(self, e):
        name = str(e.args[0])
        if name[0] == "q":
            index = name[1]
            return "c%s" % index
        else:
            return str(e)

    def _print_tan(self, e):
        name = str(e.args[0])
        if name[0] == "q":
            index = name[1]
            return "t%s" % index
        else:
            return str(e)

def pydy_str(e):
    p = PyDyStrPrinter()
    return p.doprint(e)

def pprint(e):
    p = PyDyPrettyPrinter()
    print p.doprint(e)

def sort_UnitVector(a, b):
    if a.frame == b.frame:
        return cmp(a.i, b.i)
    else:
        return cmp(len(a.frame.ref_frame_list),
               len(b.frame.ref_frame_list))
