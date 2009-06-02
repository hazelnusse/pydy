from sympy import Symbol, symbols, Basic, Derivative, Mul, Pow, Matrix, sin, \
        cos, S, eye, Add, trigsimp, expand
from sympy.printing.pretty.pretty import PrettyPrinter, xsym, pprint, pretty
# here is how to get a nice symbol for multiplication:
# print xsym("*")
from sympy.printing.str import StrPrinter
import sys
sys.displayhook = pprint

global _compact_trig
_compact_trig = False
e1 = Matrix([1, 0, 0])
e2 = Matrix([0, 1, 0])
e3 = Matrix([0, 0, 1])
e1n = Matrix([-1, 0, 0])
e2n = Matrix([0, -1, 0])
e3n = Matrix([0, 0, -1])
zero = Matrix([0, 0, 0])
#eye = Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
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
        elif i == -1:
            self.v['sym'] = Symbol('-'+s.lower()+str(abs(i)))
            self.v['num'] = e1n
        elif i == -2:
            self.v['sym'] = Symbol('-'+s.lower()+str(abs(i)))
            self.v['num'] = e2n
        elif i == -3:
            self.v['sym'] = Symbol('-'+s.lower()+str(abs(i)))
            self.v['num'] = e3n
        elif i == 0:
            self.v['sym'] = Symbol(s.lower()+str(0))
            self.v['num'] = zero

    def __str__(self):
        return pydy_str(self)

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
                u[0] = trigsimp(u[0])
                u[1] = trigsimp(u[1])
                u[2] = trigsimp(u[2])
                return Vector(u[0]*frame[1] + u[1]*frame[2] + u[2]*frame[3])

    def dot(self, other):
        if isinstance(other, UnitVector):
            c = other.express(self.frame)
            if isinstance(c, UnitVector):
                return (self.v['num'].T * c.v['num'])[0]
            elif isinstance(c, Vector):
                s = S(0)
                for k in c.dict.keys():
                    if self == k:
                        s += c.dict[k]
                return trigsimp(s)
            else:
                raise NotImplementedError()
        elif isinstance(other, Vector):
            s = S(0)
            for k in other.dict.keys():
                s += other.dict[k]*self.dot(k)
            return trigsimp(s)
        else:
            raise NotImplementedError()

    def cross(self, other):
        def cross_with_Vector(self, c):            # local function
            cp = {}
            for k in c.dict.keys():
                term = self.cross(k)
                coef = c.dict[k]
                if isinstance(term, UnitVector):
                    if cp.has_key(term):
                        cp[term] += coef
                    else:
                        cp.update({term : coef})
                elif isinstance(term, Vector):
                    for kt in term.dict.keys():
                        if cp.has_key(kt):
                            cp[kt] += coef*term.dict[kt]
                        else:
                            cp.update({kt : coef*term.dict[kt]})
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
            #print 'should be q3.diff(t)*A[3]',self.frame.get_omega(diff_frame)
            return cross(self.frame.get_omega(diff_frame), self)
        else:
            raise NotImplementedError()


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
        else:
            vdict = self.parse_terms(v)
            for k in vdict.keys():
                if vdict[k] == 0:  vdict.pop(k)
            self.dict = vdict

    def __str__(self):
        return pydy_str(self)

    def __add__(self, other):
        """Adds two Vector objects and return a new Vector object.
        v1 + v2   <----> v1.__add__(v2)
        """
        if isinstance(other, Vector):
            #print '__add__ was called, with', self, 'and', other
            s1 = set(self.dict.keys())
            s2 = set(other.dict.keys())
            sum = {}
            #print 's1', s1
            #print 's2', s2
            for k in s1.intersection(s2):
                sum.update({k: trigsimp(self.dict[k] + other.dict[k])})
            for k in s1.difference(s2):
                sum.update({k: trigsimp(self.dict[k])})
            for k in s2.difference(s1):
                sum.update({k: trigsimp(other.dict[k])})

            for k in sum.keys():
                sum[k] = trigsimp(sum[k])
                if sum[k] == 0: sum.pop(k)

            if len(sum) == 1 and sum[sum.keys()[0]] == 1:
                return sum.keys()[0]
            else:
                return Vector(sum)
        elif isinstance(other, UnitVector):
            return self + Vector({other: S(1)})
        elif isinstance(other, Add):
            return self + Vector(other)
        elif isinstance(other, Mul):
            return self + Vector(other)
        else:
            raise NotImplementedError()

    def __sub__(self, other):
        """Subtracts two Vector objects and return a new Vector object.
        v1 - v2   <----> v1.__sub__(v2)
        """
        if isinstance(other, Vector):
            #print 'v2 is a Vector'
            s1 = set(self.dict.keys())
            s2 = set(other.dict.keys())
            difference = {}

            for k in s1.intersection(s2):
                difference.update({k: trigsimp(self.dict[k] - other.dict[k])})
            for k in s1.difference(s2):
                difference.update({k: trigsimp(self.dict[k])})
            for k in s2.difference(s1):
                difference.update({k: -trigsimp(other.dict[k])})
            if not all(difference.values()):
                return Vector({})
            else:
                return Vector(difference)
        elif isinstance(other, UnitVector):
            return self - Vector({other: S(1)})
        elif isinstance(other, Add):
            return self - Vector(other)
        elif isinstance(other, Mul):
            return self - Vector(other)
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

    def __neg__(self):
        n = self.dict.copy()
        for k in n:
            n[k] *= -S(1)
        return Vector(n)

    def cross(self, other):
        if isinstance(other, Vector):
            vcp = {}
            for k in self.dict.keys():
                for ko in other.dict.keys():
                    kcrossko = k.cross(ko)
                    if isinstance(kcrossko, UnitVector):
                        if vcp.has_key(kcrossko):
                            vcp[kcrossko] += self.dict[k]*other.dict[ko]
                        else:
                            vcp.update({kcrossko: self.dict[k]*other.dict[ko]})
                    else:
                        for uv_term in kcrossko.dict.keys():
                            if vcp.has_key(uv_term):
                                vcp[uv_term] +=\
                                        self.dict[k]*other.dict[ko]*kcrossko.dict[uv_term]
                            else:
                                vcp.update({uv_term:
                                    self.dict[k]*other.dict[ko]*kcrossko.dict[uv_term]})
            return Vector(vcp)
        elif isinstance(other, UnitVector):
            vcp = {}
            for k in self.dict.keys():
                k_cross_other = k.cross(other)
                if isinstance(k_cross_other, UnitVector):
                    if vcp.has_key(k_cross_other):
                        vcp[k_cross_other] += self.dict[k]
                    else:
                        vcp.update({k_cross_other: self.dict[k]})
                else:
                    for uv_term in k_cross_other.dict.keys():
                        if vcp.has_key(uv_term):
                            vcp[uv_term] += self.dict[k]*k_cross_other.dict[uv_term]
                        else:
                            vcp.update({uv_term:
                                self.dict[k]*k_cross_other.dict[uv_term]})
            return Vector(vcp)
        elif isinstance(other, Mul) or isinstance(other, Add):
            return self.cross(Vector(other))
        else:
            raise NotImplementedError()

    def dot(self, other):
        if isinstance(other, Vector):
            s = S(0)
            for k in self.dict.keys():
                for ko in other.dict.keys():
                    s += self.dict[k]*other.dict[ko]*k.dot(ko)
            return s
        elif isinstance(other, UnitVector):
            s = S(0)
            for k in self.dict.keys():
                s += self.dict[k]*k.dot(other)
            return s
        elif isinstance(other, Mul) or isinstance(other, Add):
            return self.dot(Vector(other))
        else:
            raise NotImplementedError()

    def dt(self, frame):
        if isinstance(frame, ReferenceFrame):
            dt_self = {}
            for k in self.dict.keys():
                # First term comes from time differentiating in the frame of
                # the UnitVector frame of k
                if dt_self.has_key(k):
                    dt_self[k] += (self.dict[k]).diff(t)
                else:
                    dt_self.update({k: (self.dict[k]).diff(t)})
                # Second term comes from the omega cross term
                t2 = k.frame.get_omega(frame).cross(k)
                if isinstance(t2, UnitVector):
                    if dt_self.has_key(t2):
                        dt_self[k] += self.dict[k]
                    else:
                        dt_self.update({t2: self.dict[k]})
                else:       # Must be a Vector
                    for term in t2.dict.keys():
                        if dt_self.has_key(term):
                            dt_self[term] += self.dict[k]*t2.dict[term]
                        else:
                            dt_self.update({term: self.dict[k]*t2.dict[term]})
            if len(dt_self) == 1:
                if dt_self.values()[0] == 1:
                    return dt_self.keys()[0]        # Return a UnitVector
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
                if new.has_key(uv_in_frame):
                    new[uv_in_frame] += self.dict[uv]
                else:
                    new.update({uv_in_frame : self.dict[uv]})

            # Case for Vectors
            elif isinstance(uv_in_frame, Vector):
                # Go through each term
                for uv_term in uv_in_frame.dict.keys():
                    if new.has_key(uv_term):
                        new[uv_term] += self.dict[uv]*uv_in_frame.dict[uv_term]
                    else:
                        new.update({uv_term :
                            self.dict[uv]*uv_in_frame.dict[uv_term]})

        for uv in new.keys():
            new[uv] = trigsimp(new[uv])
            if new[uv] == 0: new.pop(uv)

        if len(new) == 1 and new.values()[0] == 1:
            return new.keys()[0]
        else:
            return Vector(new)

    def mag(self):
        m = 0
        for k in self.dict.keys():
            m += self.dict[k]**S(2)
        return m**(S(1)/S(2))

    def mag_sqr(self):
        m = 0
        for k in self.dict.keys():
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
            v.expand()  #  Could expand out to an Add instance
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
            v.expand()
            #  I don't think this will ever be entered into.
            #  You would have to have something like A[1]*A[2],
            #  which isn't a valid vector expression.
            #  Or q1*q2, which is a scalar expression.
            for b in v.args:
                if isinstance(b, UnitVector):
                    return {b: v.coeff(b)}
        elif isinstance(v, Add):
            #v.expand()
            terms = {}
            for add_term in v.args:
                if isinstance(add_term, Mul):
                    add_term_dict = self.parse_terms(add_term)
                elif isinstance(add_term, Pow):
                    add_term_dict = self.parse_terms(add_term)
                elif isinstance(add_term, UnitVector):
                    add_term_dict = {add_term: S(1)}
                elif isinstance(add_term, Vector):
                    add_term_dict = add_term.dict
                else:
                    raise NotImplementedError()
                for k in add_term_dict.keys():
                    if terms.has_key(k):
                        terms[k] += add_term_dict[k]
                    else:
                        terms.update(add_term_dict)
            return terms
        else:
            return NotImplemented

    def subs(self, subs_dict):
        new_dict = {}
        for k in self.dict.keys():
            new_dict.update({k: trigsimp(expand(self.dict[k].subs(subs_dict)))})
        return Vector(new_dict)

class Point:
    def __init__(self, s, r=None, frame=None):
        self.name = s
        # When instantiated by ReferenceFrame
        if r == S(0) and isinstance(frame, ReferenceFrame):
            self.pos = S(0)                # Should only happen for the base
            self.vel = S(0)                # Newtonian/Inertial Frame
            self.acc = S(0)
            self.NewtonianFrame = frame
        elif r != None and frame == None:  # When instantiated by locate method
            self.pos = r
        else:
            raise NotImplementedError()


    def locate(self, s, r, frame=None):
        r = Vector(r)
        if isinstance(r, UnitVector) or isinstance(r, Vector):
            if frame == None:
                newpoint = Point(s, r)
                newpoint.NewtonianFrame = self.NewtonianFrame
                #print 'r', r, 'type(r)', type(r)
                newpoint.vel = r.dt(newpoint.NewtonianFrame)
            elif isinstance(frame, ReferenceFrame):
                newpoint = Point(s, r)
                newpoint.NewtonianFrame = self.NewtonianFrame
                wn = frame.get_omega(newpoint.NewtonianFrame)
                #print 'computing vel'
                #print 'omega', wn
                #print 'self.vel', self.vel
                #print 'w x r', cross(wn, r)
                newpoint.vel = self.vel + cross(wn, r)
                #print 'self.vel + w x r = ', newpoint.vel
            else:
                raise NotImplementedError()
        return newpoint

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

    Method 1 typically is used to create the 'base'frame from which all other
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
        if frame == None:
            self.ref_frame_list = [self]
            self.O = Point(s, S(0), self)
        else:
            self.ref_frame_list = frame.ref_frame_list[:]
            self.ref_frame_list.insert(0,self)

        self.name = s
        self.triad = [UnitVector(self, i) for i in (1,2,3)]
        self.transforms = {}
        self.parent = frame
        self.W = {}
        if omega != None:
            #print 'omega in init', omega
            #print 'self.W before set_omega', self.W
            self.set_omega(omega, self.parent)
            #print 'omegas after set_omega',self.W
            #print 'omegas after set_omega, using get_omega', \
            #    self.get_omega(frame)
            #print 'self', self, 'parent frame', frame
            #print 'frame',frame.set_omega
            #print 'parent frame omega', frame.W
            #print 'arguments to frame.set_omega()', -omega, self
            frame.set_omega(-omega, self, force=True)
            #print 'after setting parent omega'
            #print 'omegas after set_omega, self.W',self.W
            #print 'omegas after set_omega, using get_omega', \
            #    self.get_omega(frame)
            #print 'parent frame omega', frame.W

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
                #print 'omega = ', omega
            elif isinstance(axis, (UnitVector, Vector)):
                raise NotImplementedError("Axis angle rotations not \
                    implemented")
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
                if len(r1t)<=2 or len(r2t)<=2:
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
        Returns the angular velocity of self relative to the frame "frame".
        """

        if self.W.has_key(frame):
            return self.W[frame]
        else:
            om = {}
            for term in self.get_omega_list(frame):
                for k in term.dict.keys():
                    if om.has_key(k):
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
    if (isinstance(v1, UnitVector) or isinstance(v1, Vector)) and \
            (isinstance(v2, UnitVector) or isinstance(v2, Vector)):
                return v1.dot(v2)
    else:
        if not (isinstance(v1, UnitVector) or isinstance(v1, Vector)):
            v1 = Vector(v1)
        if not (isinstance(v2, UnitVector) or isinstance(v2, Vector)):
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
        c = {}
        for k in v.dict.keys():
            s = v.dict[k].coeff(scalar)
            if s == None:
                continue
            else:
                c.update({k: s})
        return Vector(c)
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
    if isinstance(v, UnitVector) or isinstance(v, Vector):
        res = v.dt(frame)
        #print res
        return res
    else:
        raise NotImplementedError()


class PyDyPrinter(StrPrinter):
    #printmethod = "_pydystr_"

    #def _print_Symbol(self, e):
    #    return "|%s|" % str(e)

    def _print_UnitVector(self, e):
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
        #print type(r)
        r += reset
        return r

    def _print_Vector(self, e):
        s = ''
        i = 0
        small_dot = "\xC2\xB7"
        if e.dict.keys() != []:
            uv_list = e.dict.keys()
            uv_list.sort(sort_UnitVector)
            for k in uv_list:
                if (e.dict[k] == 1) or (e.dict[k] == -1):
                    if i == 0:
                        if e.dict[k] == 1: sign = ''
                        if e.dict[k] == -1: sign = '-'
                        s += sign + k.__str__()
                        i += 1
                    else:
                        if int(abs(e.dict[k])/e.dict[k]) == 1: sign = ''
                        if int(abs(e.dict[k])/e.dict[k]) == -1: sign = '-'
                        s += ' ' + sign + ' ' + k.__str__()
                else:
                    if i == 0:
                        if isinstance(e.dict[k], Add):
                            s += ('(' + self.doprint(e.dict[k]) +
                                    ')' + small_dot + k.__str__())
                        else:
                            s += (self.doprint(e.dict[k]) +
                                small_dot + k.__str__())
                        i += 1
                    else:
                        if isinstance(e.dict[k], Add):
                            s += (' + (' + self.doprint(e.dict[k])
                                + ')' + small_dot + k.__str__())
                        else:
                            if str(e.dict[k])[0] == '-':
                                sign = '-'
                                s += (' ' + sign + ' ' +
                                        self.doprint(e.dict[k])[1:] + small_dot + k.__str__())
                            else:
                                sign = '+'
                                s += (' ' + sign + ' ' + self.doprint(e.dict[k])
                                    + small_dot + k.__str__())
            return s
        else:
            return "\033[1m" + "0" + "\033[0;0m"


    def _print_Function(self, e):
        return "%s" % str(e.func)

    def _print_Derivative(self, e):
        return "%s'" % str(e.args[0].func)

    def _print_Mul(self, e):
        s = ''
        i = 1
        N = len(e.args)
        e_ts = trigsimp(expand(e))
        for a in e_ts.args:
            if i == 0:
                if a == -1:
                    s += '-'
                else:
                    s += self.doprint(a) + '\xC2\xB7'
                i += 1
            elif i < N:
                if a == -1:
                    s += '-'
                else:
                    s += self.doprint(a) + '\xC2\xB7'
                i += 1
            else:
                s += self.doprint(a)
        return s

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
    pp = PyDyPrinter()
    return pp.doprint(e)

def sort_UnitVector(a, b):
    if a.frame == b.frame:
        return cmp(a.i, b.i)
    else:
        return cmp(len(a.frame.ref_frame_list),
               len(b.frame.ref_frame_list))
