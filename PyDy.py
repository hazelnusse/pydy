import numpy as N
from sympy import Symbol, Basic, Mul, Pow, Matrix, sin, cos, S, eye
e1 = Matrix([1, 0, 0])
e2 = Matrix([0, 1, 0])
e3 = Matrix([0, 0, 1])
e1n = Matrix([-1, 0, 0])
e2n = Matrix([0, -1, 0])
e3n = Matrix([0, 0, -1])
zero = Matrix([0, 0, 0])


class UnitVector(Basic):
    """A standard unit vector  with a symbolic and a numeric representation"""

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


#    def __neg__(self):
#        return UnitVector(self.name,-self.i)

    def _sympystr_(self):
        return str(self.v['sym']) + ">"

    #  Cross product
    def cross(self,other):
        if isinstance(other, UnitVector):
            v1xv2 = N.cross(self.v['num'],other.v['num'])

            if (v1xv2 == e1).all():
                return UnitVector(self.frame, 1)
            elif (v1xv2 == e2).all():
                return UnitVector(self.frame, 2)
            elif (v1xv2 == e3).all():
                return UnitVector(self.frame, 3)
            elif (v1xv2 == e1n).all():
                return -UnitVector(self.frame, 1)
            elif (v1xv2 == e2n).all():
                return -UnitVector(self.frame, 2)
            elif (v1xv2 == e3n).all():
                return -UnitVector(self.frame, 3)
            elif (v1xv2 == zero).all():
                return UnitVector(self.frame, 0)
        else:
            raise NotImplementedError()


class ReferenceFrame:
    """A standard reference frame with 3 dextral orthonormal vectors"""

    def __init__(self, s, matrix=None, frame=None):
        """
        ReferenceFrame('B', matrix, A)

        The matrix represents A_B, e.g. how to transform a vector from B to A.
        """
        self.name = s
        self.triad = [UnitVector(self, i) for i in (1,2,3)]
        self.transforms = {}
        self.parent = frame
        if frame is not None:
            self.append_transform(frame, matrix)
            frame.append_transform(self, matrix.T)

    def __getitem__(self, i):
        return self.triad[i-1]

    def append_transform(self, frame, matrix):
        """
        Gives us a transform to the frame "frame".
        """
        # We just append it to our "transforms" dict.
        self.transforms[frame] = matrix

    def rotate(self, name, axis, angle):
        if axis == 1:
            matrix = Matrix([
                [1, 0, 0],
                [0, cos(angle), -sin(angle)],
                [0, sin(angle), cos(angle)],
                ])
        elif axis == 2:
            matrix = Matrix([
                [cos(angle), 0, sin(angle)],
                [0, 1, 0],
                [-sin(angle), 0, cos(angle)],
                ])
        elif axis == 3:
            matrix = Matrix([
                [cos(angle), -sin(angle), 0],
                [sin(angle), cos(angle), 0],
                [0, 0, 1],
                ])
        else:
            raise ValueError("wrong axis")
        return ReferenceFrame(name, matrix, self)

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
        elif self.transforms.has_key(frame):
            return [self, frame]
        else:
            # r1, r2 contains all the frames up to the top
            r1 = [self]
            a = self.parent
            while a is not None:
                r1.append(a)
                a = a.parent
            r2 = [frame]
            a = frame.parent
            while a is not None:
                r2.append(a)
                a = a.parent

            # we strip the common right end of both lists:
            i = 1
            while not (len(r1) == i-1 or len(r2) == i-1) and r1[-i] == r2[-i]:
                i += 1
            i -= 1

            #print "r1, r2, i:", r1, r2, i
            middle = r1[-i+1]
            #print middle
            r1 = r1[:-(i-1)]
            r2 = r2[:-i]
            #print "stripped r1, r2:", r1, r2
            # e.g.: r1 == [C, B]
            # r2 == [F, E, D]
            # middle == A
            return r1 + list(reversed(r2))

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

def dot(v1,v2):
    if isinstance(v1, UnitVector) and isinstance(v2, UnitVector):
        return (v1.v['num'].T*v2.v['num'])[0]
    else:
        v1v2 = (v1*v2).expand()
        #print v1v2.args
        e = 0
        for a in v1v2.args:
            c, v1, v2 = identify(a)
            #print c, v1, v2
            if v1 is None or v2 is None:
                pass
            else:
                e += c*dot(v1, v2)
        return e


def identify(a):
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
            c = a.coeff(v1*v2)
            return c, v1, v2

    return a, None, None

def cross(v1, v2):
    #return v1.cross(v2)
    A = v1.frame
    B = v2.frame
    matrices = A.get_rot_matrices(B)
    # first vector as a list of 3 numbers in the "v2" inertial frame:
    u = Matrix(v1.v['num'])
    for m in matrices:
        #print m
        u = m*u
    #print u
    # second vector:
    v = Matrix(v2.v['num'])
    c1 = u[1]*v[2] - u[2]*v[1]
    c2 = -(u[0]*v[2] - u[2]*v[0])
    c3 = u[0]*v[1] - u[1]*v[0]
    return c1*B[1] + c2*B[2] + c3*B[3]
