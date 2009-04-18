from sympy import Symbol, Basic, Derivative, Mul, Pow, Matrix, sin, cos, S, eye, Add, \
        trigsimp
from numpy import sign
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

    def _sympystr_(self):
        return str(self.v['sym']) + ">"

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

    def express(self, F):
        '''
        Express a UnitVector in the reference frame "F".
        '''
        if self.frame == F:
            return self
        else:
            matrices = self.frame.get_rot_matrices(F)
            if len(matrices) == 1 and matrices[0]*self.v['num'] == \
                    self.v['num']:
                return F[self.i]
            else:
                u = self.v['num']
                for m in reversed(matrices):
                    u = m*u
                return Vector(trigsimp(u[0])*F[1] + trigsimp(u[1])*F[2] + \
                        trigsimp(u[2])*F[3])

    def dot(self, other):
        if isinstance(other, UnitVector):
            c = other.express(self.frame)
            if isinstance(c, UnitVector):
                return (self.v['num'].T * c.v['num'])[0]
            elif isinstance(c, Vector):
                s = S(0)
                for k in c.dict.keys():
                    if self == k:
                        s += trigsimp(c.dict[k])
                return trigsimp(s)
            else:
                raise NotImplementedError()
        elif isinstance(other, Vector):
            s = S(0)
            for k in other.dict.keys():
                s += trigsimp(trigsimp(other.dict[k])*trigsimp(self.dot(k)))
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
            self.dict = v
        else:
            #print "Not a dict"
            self.dict = self.parse_terms(v)

    def express(self, Frame):
        """
        Expresses self in the basis vectors of 'Frame'.
        """

        self_Frame_dict = {}

        for unit_vector_term in self.dict.keys():
            uv_in_Frame = unit_vector_term.express(Frame)
            if isinstance(uv_in_Frame, UnitVector):
                if self_Frame_dict.has_key(uv_in_Frame):
                    self_Frame_dict[uv_in_Frame] += S(1)
                else:
                    self_Frame_dict.update({uv_in_Frame : \
                            trigsimp(self.dict[unit_vector_term])})
            elif isinstance(uv_in_Frame, Vector):
                for uv_term in uv_in_Frame.dict.keys():
                    if self_Frame_dict.has_key(uv_term):
                        self_Frame_dict[uv_term] +=\
                                trigsimp(self.dict[unit_vector_term] * \
                                uv_in_Frame.dict[uv_term])
                    else:
                        self_Frame_dict.update({uv_term : \
                                trigsimp(self.dict[unit_vector_term] \
                                * uv_in_Frame.dict[uv_term])})

        for key in self_Frame_dict.keys():
            self_Frame_dict[key] = trigsimp(self_Frame_dict[key])
            if self_Frame_dict[key] == 0:
                self_Frame_dict.pop(key)
        if len(self_Frame_dict.keys()) == 1:
            if self_Frame_dict[self_Frame_dict.keys()[0]] == 1:
                return self_Frame_dict.keys()[0]
            else:
                return Vector(self_Frame_dict)
        else:
            return Vector(self_Frame_dict)

    def _sympystr_(self):
        """
        Sympy printing of Vector objects.
        """

        s = ''
        i = 0
        if self.dict != {}:
            for k in self.dict.keys():
                if (self.dict[k] == 1) or (self.dict[k] == -1):
                    if i == 0:
                        #print 'hello'
                        if self.dict[k] == 1: sign = ''
                        if self.dict[k] == -1: sign = '-'
                        #s += k._sympystr_()
                        s += sign + k._sympystr_()
                        i += 1
                    else:
                        if int(abs(self.dict[k])/self.dict[k]) == 1: sign = ''
                        if int(abs(self.dict[k])/self.dict[k]) == -1: sign = '-'
                        s += ' ' + sign + ' ' + k._sympystr_()
                else:
                    if i == 0:
                        s += str(self.dict[k]) + '*' + k._sympystr_()
                        i += 1
                    else:
                        if str(self.dict[k])[0] == '-':
                            sign = '-'
                            s += ' ' + sign + ' ' + str(self.dict[k])[1:] + '*' + k._sympystr_()
                        elif str(self.dict[k])[0] != '-':
                            sign = '+'
                            s += ' ' + sign + ' ' + str(self.dict[k]) + '*' + k._sympystr_()
            return s
        else:
            return '0>'


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
                #print "b = ", b, "type(b): ", type(b)
                #print "b parsed = ", parse_terms(b)
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

    def __add__(self, other):
        """
        Adds two Vector objects and return a new Vector object.
        v1 + v2   <----> v1.__add__(v2)
        """
        if isinstance(other, Vector):
            s1 = set(self.dict.keys())
            s2 = set(other.dict.keys())
            sum = {}

            for k in s1.intersection(s2):
                sum.update({k: trigsimp(self.dict[k] + other.dict[k])})
            for k in s1.difference(s2):
                sum.update({k: trigsimp(self.dict[k])})
            for k in s2.difference(s1):
                sum.update({k: trigsimp(other.dict[k])})
            if not all(sum.values()):
                return Vector({})
            else:
                return Vector(sum)
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
        """
        Subtracts two Vector objects and return a new Vector object.
        v1 - v2   <----> v1.__sub__(v2)
        """
        if isinstance(other, Vector):
            print 'v2 is a Vector'
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
        '''
        Compares two Vector objects for equality.
        '''
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
                        #print self_in_first_key_frame.dict
                        #print other_in_first_key_frame.dict
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
            #print "hello"
            other_as_Vector = Vector(other)
            return self == other_as_Vector

    def __neg__(self):
        n = self.dict
        for k in n:
            n[k] *= -1
        return Vector(n)

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
                        for uv_term in kcrossko.keys():
                            if vcp.has_key(uv_term):
                                vcp[uv_term] +=\
                                        self.dict[k]*other.dict[ko]*kcrossko[uv_term]
                            else:
                                vcp.update({uv_term: 
                                    self.dict[k]*other.dict[ko]*kcrossko[uv_term]})
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

    def mag(self):
        m = 0
        for k in self.dict.keys():
            m += self.dict[k]**2
        return m**(1/2)

class Point:
    def __init__(self, s, r=None, frame=None):
        self.name = s
        if r == S(0) and isinstance(frame, ReferenceFrame):    # When instantiated by ReferenceFrame
            self.pos = S(0)                # Should only happen for the base
            self.vel = S(0)                # Newtonian Frame
            self.acc = S(0)
            self.NewtonianFrame = frame
        elif r != None and frame == None:  # When instantiated by locate method
            self.pos = r
        else:
            raise NotImplementedError()


    def locate(self, s, r, frame=None):
        if isinstance(r, UnitVector) or isinstance(r, Vector):
            if Frame == None:
                newpoint = Point(s, r)
                newpoint.NewtonianFrame = self.NewtonianFrame
                newpoint.vel = r.dt(newpoint.NewtonianFrame)
            elif isinstance(frame, ReferenceFrame):
                newpoint = Point(s, r)
                newpoint.NewtonianFrame = self.NewtonianFrame
                newpoint.vel = self.vel +\
                        frame.get_omega(newpoint.NewtonianFrame).cross(r)
            else:
                raise NotImplementedError()



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
            self.origin = Point(s, S(0), self)
        else:
            #print "inside __init__: ", frame.ref_frame_list
            self.ref_frame_list = frame.ref_frame_list[:]
            self.ref_frame_list.insert(0,self)

        self.name = s
        self.triad = [UnitVector(self, i) for i in (1,2,3)]
        #print "self.triad = ", self.triad, "type(self.triad) = ", \
                #        type(self.triad)
        self.transforms = {}
        self.parent = frame
        self.W = {}
        if omega != None:
            self.set_omega(omega, self.parent)
            self.parent.set_omega(-omega, self, force=True)

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
        """
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
        if axis == 1:
            matrix = Matrix([
                [1, 0, 0],
                [0, cos(angle), -sin(angle)],
                [0, sin(angle), cos(angle)],
                ])
            omega = angle.diff(t)*self[axis]
        elif axis == 2:
            matrix = Matrix([
                [cos(angle), 0, sin(angle)],
                [0, 1, 0],
                [-sin(angle), 0, cos(angle)],
                ])
            omega = angle.diff(t)*self[axis]
        elif axis == 3:
            matrix = Matrix([
                [cos(angle), -sin(angle), 0],
                [sin(angle), cos(angle), 0],
                [0, 0, 1],
                ])
            omega = angle.diff(t)*self[axis]
        else:
            raise ValueError("wrong axis")
        return ReferenceFrame(name, matrix, self, omega)

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

    def set_omega(self, W, frame, force=False):
        """
        Sets the angular velocity Omega with respect to the frame "frame".
        """
        if self.W == {} or force:
            self.W[frame] = W
        else:
            raise ValueError("set_omaga has already been called.")

    def get_omega(self, frame):
        """
        Returns the angular velocity of self relative to the frame "frame".
        """

        if self.W.has_key(frame):
            return self.W[frame]
        else:
            return sum(self.get_omega_list(frame))

    def get_omega_list(self, frame):
        """
        Returns a list of simple angular velocities from self to frame.
        """
        frames = self.get_frames_list(frame)
        if frames == [self]:
            return [S(0)]
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


'''
def dot(v1,v2):
    if isinstance(v1, UnitVector) and isinstance(v2, UnitVector):
        B = v2.frame
        u = express(v1, B)
        u1, u2, u3 = expression2vector(u, B)
        # second vector:
        v = Matrix(v2.v['num'])
        return dot_vectors((u1, u2, u3), v)
    else:
        v1v2 = (v1*v2).expand()
        #print "DOT:", v1v2
        if isinstance(v1v2, Add):
            #print v1v2.args
            e = 0
            for a in v1v2.args:
                c, v1, v2 = identify(a)
                #print "IDENTIFY", c, v1, v2, a
                if v1 is None or v2 is None:
                    raise Exception("!")
                else:
                    #print "c", c
                    #print "dot", dot(v1, v2)
                    e += c*dot(v1, v2)
            return e
        elif isinstance(v1v2, Mul):
            c, v1, v2 = identify(v1v2)
            if v1 is None or v2 is None:
                raise NotImplementedError()
            else:
                e = c*dot(v1, v2)
            return e
        elif v1v2 == 0:
            return v1v2
        else:
            raise NotImplementedError()
'''

'''
def parse_terms(v):
    #"""
    #Given a Sympy expression with UnitVector terms, return a dictionary whose
    #keys are the UnitVectors and whose values are the coeefficients of the
    #UnitVectors
    #"""
    v.expand()
    if v == 0:
        return {}
    elif isinstance(v, UnitVector):
        return {v : S(1)}
    elif isinstance(v, Mul):
        for b in v.args:
            if isinstance(b, UnitVector):
                return {b: v.coeff(b)}
        return NotImplemented
    elif isinstance(v, Pow):
        #  I don't think this will ever be entered into.
        #  You would have to have something like A[1]*A[2],
        #  which isn't a valid vector expression.
        #  Or q1*q2, which is a scalar expression.
        for b in v.args:
            if isinstance(b, UnitVector):
                return {b: v.coeff(b)}
    elif isinstance(v, Add):
        terms = {}
        for b in v.args:
            #print "b = ", b, "type(b): ", type(b)
            #print "b parsed = ", parse_terms(b)
            bp = parse_terms(b)
            #print "bp.keys(): ", bp.keys()[0]
            if terms.has_key(bp.keys()[0]):
                terms[bp.keys()[0]] += bp[bp.keys()[0]]
            else:
                terms.update(parse_terms(b))
        return terms
    else:
        return NotImplemented
'''


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
'''
def express(v, frame):
    """
    Express "v" in the reference frame "frame".
    """
    if isinstance(v, UnitVector):
        matrices = v.frame.get_rot_matrices(frame)
        #print matrices
        u = Matrix(v.v['num'])
        #XXX: Not sure why reversed is necessary
        for m in reversed(matrices):
            u = m*u
        return vector2expression(u, frame)
    elif isinstance(v, Add):
        e = 0
        for a in v.args:
            c, v1 = identify_v1(a)
            #print c, v1
            if v1 is None:
                pass
            else:
                e += c*express(v1, frame)
        e = e.expand()
        u = expression2vector(e, frame)
        u = [trigsimp(x) for x in u]
        return vector2expression(u, frame)
    elif isinstance(v, Mul):
        c, v1 = identify_v1(v)
        return (c*express(v1, frame)).expand()
    elif v == 0:
        return v
    else:
        #print "XXX", v
        raise NotImplementedError()
'''

def express(v, frame):
    if (isinstance(v, UnitVector) or isinstance(v, Vector)) and \
            (isinstance(frame, ReferenceFrame)):
        return v.express(frame)
    else:
        raise NotImplementedError()


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

def dt(u, frame, t):
    if isinstance(u, Add):
        r = 0
        for a in u.args:
            c, v = identify_v1(a)
            #print "c = ", c, "type(c)", type(c), "v = ",v, "type(v)", type(v)
            dc_dt = c.diff(t)
            W = v.frame.get_omega(frame)
            #print "W = ", W, "type(W):", type(W)
            if W == 0:
                print "W = 0"
                print "r = ", r, "dc_dt = ", dc_dt, "v = ", v
                r += dc_dt * v
                continue
            #print "dc_dt*v + W x v", dc_dt, v, W, v, cross(W, v)
            r += dc_dt * v + c*cross(W, v)
        r = r.expand()
        return r

    elif isinstance(u, Mul):
        c, v = identify_v1(u)
        #print "c = ", c, "v = ", v
        dc_dt = c.diff(t)
        W = v.frame.get_omega(frame)
        #print "W = ", W
        #print "W x v = ", cross(W,v)
        r = dc_dt * v + c*cross(W, v)
        return r
    else:
        raise NotImplementedError()

from sympy.printing.pretty.pretty import PrettyPrinter, xsym
# here is how to get a nice symbol for multiplication:
# print xsym("*")
from sympy.printing.str import StrPrinter

class PyDyPrinter(StrPrinter):
    printmethod = "_pydystr_"

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
        r = "%s%s%s" % (bold, name, reset)
        if index == "1":
            r += one
        elif index == "2":
            r += two
        elif index == "3":
            r += three
        return r

    def _print_sin(self, e):
        name = str(e.args[0])
        if name[0] == "q":
            index = name[1]
            return "s_%s" % index
        else:
            return str(e)


def print_pydy(e):
    pp = PyDyPrinter()
    print pp.doprint(e)
