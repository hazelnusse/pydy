def convert_geometric(x):
    """Convert geometric parameters of Meijaard et al. to those of Peterson et
    al.

    Input x (list):
        w:      Wheelbase
        c:      Trail
        lmbda:  Head angle complement
        rr:     Rear wheel radius
        rrt:    Rear tire casing radius
        rf:     Front wheel radius
        rft:    Front tire casing radius

    Output (list):
        lr:     Rear offset
        lf:     Front offset
        ls:     Steer axis offset
    """
    w, c, lmbda, rr, rrt, rf, rft = x
    lr = (w + c)*cos(lmbda)-(rr + rrt)*sin(lmbda)
    lf = (rf + rft)*sin(lmbda)-c*cos(lmbda)
    ls = w*sin(lmbda) + (rr + rrt - rf - rft)*cos(lmbda)
    return lr, lf, ls

def convert_inertial(x):
    """Convert the inertial parameters of Meijaard et al. to those of Peterson
    et al.

    Input x (list/array):
        mr:     Rear wheel mass
        mb:     Frame/rider mass
        mh:     Fork/handlbar mass
        mf:     Front wheel mass
        IRxx:   Rear wheel in plane moment of inertia
        IBxx:   Frame/rider in plane moment of inertia (x)
        IBzz:   Frame/rider in plane moment of inertia (z)
        IHxx:   Fork/handlebar in plane moment of inertia (x)
        IHzz:   Fork/handlebar in plane moment of inertia (z)
        IFxx:   Front wheel in plane moment of inertia
        xb:     Frame/rider COM distance from R.W. contact in x-direction
        zb:     Frame/rider COM distance from R.W. contact in z-direction
        xh:     Fork/handlebar COM distance from R.W. contact in x-direction
        zh:     Fork/handlebar COM distance from R.W. contact in z-direction
        lmbda:  Head angle complement

    Output:

    """
    mr, mb, mh, mf, IRxx, IBxx, IBzz, IHxx, IHzz, IFxx, lmbda = x
    # Create variables for Peterson naming conventions
    mc = mr
    md = mb
    me = mh
    # Combined masses

    mcd = mr + mb
    mef = mh + mf
    # Location of rear assembly mass center relative to rear wheel center
    #l1 = ...
    #l2 = ...
    l1_t = md * l1 / mcd
    l2_t = md * l2 / mcd
    # Location of front assembly mass center relative to front wheel center
    #l3 = ...
    #l4 = ...
    l3_t = me * l3 / mef
    l4_t = me * l4 / mef



