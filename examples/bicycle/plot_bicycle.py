import bicycle_lib as bl
from convert_parameters import convert_params as convert
from numpy import pi, array, arange
from scipy.integrate import odeint

w = 1.02
c = 0.08
lmbda = pi/10.
g = 9.81
rr = 0.3
rrt = 0.0
rf = 0.35
rft = 0.0
xb = 0.3
zb = -0.9
xh = 0.9
zh = -0.7
mr = 2.0
mb = 85.0
mh = 4.0
mf = 3.0
IRxx = 0.0603
IRyy = 0.12
IBxx = 9.2
IByy = 11.0
IBzz = 2.8
IBxz = 2.4
IHxx = 0.05892
IHyy = 0.06
IHzz = 0.00708
IHxz = -0.00756
IFxx = 0.1405
IFyy = 0.28

params_mj = [rr, rrt, rf, rft, w, c, lmbda, xb, zb, xh, zh, mr, mb, mh, mf,
        IRxx, IRyy, IBxx, IByy, IBzz, IBxz, IHxx, IHyy, IHzz, IHxz, IFxx, IFyy,
        g]


params_lp = array(convert(params_mj))

param_names = ['rr', 'rrt', 'rf', 'rft', 'lr', 'ls', 'lf', 'l1', 'l2', 'l3',
        'l4', 'mcd', 'mef', 'IC22', 'ICD11', 'ICD22', 'ICD33', 'ICD13',
        'IEF11', 'IEF22', 'IEF33', 'IEF13', 'IF22', 'g']
for name, value in zip(param_names, params_lp):
    print name, '=', value
