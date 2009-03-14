import numpy as N
from sympy import Symbol, Basic
e1 = N.array([1.0,0.0,0.0])
e2 = N.array([0.0,1.0,0.0])
e3 = N.array([0.0,0.0,1.0])
e1n = N.array([-1.0,0.0,0.0])
e2n = N.array([0.0,-1.0,0.0])
e3n = N.array([0.0,0.0,-1.0])
zero = N.array([0.0,0.0,0.0])


class UnitVector(Basic):
    """A standard unit vector  with a symbolic and a numeric representation"""    
    def __init__(self,s,i=0): #=-1,num=None):
        self.name = s    # Parent reference frame
        self.i = i
        self.v = {}
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
            

    def __neg__(self):
        return UnitVector(self.name,-self.i)

    def __repr__(self):
        return self.v['sym'].__repr__()

    #  Cross product
    def cross(self,other):
        if isinstance(other, UnitVector):
            v1xv2 = N.cross(self.v['num'],other.v['num'])

            if (v1xv2 == e1).all():
                return UnitVector(self.name,1)
            elif (v1xv2 == e2).all():
                return UnitVector(self.name,2)
            elif (v1xv2 == e3).all():
                return UnitVector(self.name,3)
            elif (v1xv2 == e1n).all():
                return -UnitVector(self.name,1)
            elif (v1xv2 == e2n).all():
                return -UnitVector(self.name,2)
            elif (v1xv2 == e3n).all():
                return -UnitVector(self.name,3)
            elif (v1xv2 == zero).all():
                return UnitVector(self.name,0)
        else:
            return NotImplemented


class ReferenceFrame:
    """A standard reference frame with 3 dextral orthonormal vectors"""
    def __init__(self,s):
        self.name = s
        self.triad = [UnitVector(s,i) for i in (1,2,3)]

    def __getitem__(self, i):
        return self.triad[i-1]
    
    
def dot(v1,v2):
    if isinstance(v1, UnitVector) and isinstance(v2, UnitVector):
        return N.dot(v1.v['num'],v2.v['num'])
    elif isinstance(other,Basic):
#            vec_expr = expand(other * Basic)
        return NotImplemeted
    

def cross(v1,v2):
    return v1.cross(v2)


#        self.name = s
#        self.i = i
#        self.v = {}
#        if i == 1:
#            self.v['sym'] = Symbol(s.lower()+str(i))
#            self.v['num'] = e1
#        elif i == 2:
#            self.v['sym'] = Symbol(s.lower()+str(i))
#            self.v['num'] = e2
#        elif i == 3:
#            self.v['sym'] = Symbol(s.lower()+str(i))
#            self.v['num'] = e3
#        if num is not None:
#            self.v['sym'] = s
#            self.v['num'] = num
            
