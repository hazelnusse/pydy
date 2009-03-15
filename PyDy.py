import numpy as N
from sympy import Symbol, Basic, Mul, Pow, Matrix, sin, cos
e1 = N.array([1.0,0.0,0.0])
e2 = N.array([0.0,1.0,0.0])
e3 = N.array([0.0,0.0,1.0])
e1n = N.array([-1.0,0.0,0.0])
e2n = N.array([0.0,-1.0,0.0])
e3n = N.array([0.0,0.0,-1.0])
zero = N.array([0.0,0.0,0.0])


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

    def get_rot_matrices(self, frame):
        """
        Returns a list of matrices to get from self to frame.
        """
        if self.transforms.has_key(frame):
            return [self.transforms[frame]]
        else:
            raise NotImplementedError()



def dot(v1,v2):
    if isinstance(v1, UnitVector) and isinstance(v2, UnitVector):
        return N.dot(v1.v['num'],v2.v['num'])
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

def cross(v1,v2):
    #return v1.cross(v2)
    A = v1.frame
    B = v2.frame
    matrices = A.get_rot_matrices(B)
    # first vector as a list of 3 numbers in the "v2" inertial frame:
    u = Matrix(v1.v['num'])
    for m in matrices:
        u *= m

    # second vector:
    v = v2.v['num']
    c1 = u[1]*v[2] - u[2]*v[1]
    c2 = -(u[0]*v[2] - u[2]*v[0])
    c3 = u[0]*v[1] - u[1]*v[0]
    return c1*B[1] + c2*B[2] + c3*B[3]
