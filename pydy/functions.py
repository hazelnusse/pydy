from sympy import Symbol

from pydy import UnitVector, Vector, ReferenceFrame
from common import e1, e2, e3, zero, t

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

def transform_matrix(B, x, x_dependent, subs_dict=None, time=None):
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

    Returns: -inv(Bd), Bi, substitution_dict

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
    if time:
        B_dummy, d = dummy_matrix(B, 'b', time=True)
    else:
        B_dummy, d = dummy_matrix(B, 'b', time=False)

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
    #Bd_det = -factor(Bd.det().expand())
    Bd_det = -simplify(Bd.berkowitz_det().expand())
    assert Bd_det != 0, "Equations are singular."

    # Form inv(Bd)
    # inv(Bd) = adjugate(Bd) / det(Bd)
    Bd_inv = zeros((m,m))
    for i in range(m):
        for j in range(m):
            if Bd_adj[i,j] != 0:
                #Bd_inv[i,j] = factor(Bd_adj[i,j]) / Bd_det
                #Bd_inv[i,j] = Bd_adj[i,j] / Bd_det
                Bd_inv[i,j] = simplify(Bd_adj[i,j]) / Bd_det
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
    requires a special ordering, it will *NOT* be preseved by converting it to
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

def dummy_matrix(mat, char, time=None):
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
        if n == 1:
            j = 0
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
                    if time:
                        ds = Symbol(char + '%d'%i, dummy=True)(t)
                    else:
                        ds = Symbol(char + '%d'%i, dummy=True)
                    d[ds] = mij
                    dr[mij] = ds
                    new_mat[i, j] = ds
        else:
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
                        if time:
                            ds = Symbol(char + '%d%d'%(i,j), dummy=True)(t)
                        else:
                            ds = Symbol(char + '%d%d'%(i,j), dummy=True)
                        d[ds] = mij
                        dr[mij] = ds
                        new_mat[i, j] = ds
    return new_mat, d

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

