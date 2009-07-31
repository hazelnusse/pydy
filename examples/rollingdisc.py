from sympy import symbols, Function, S, solve, simplify, \
        collect, Matrix, lambdify, trigsimp, expand, Eq, pretty_print

from pydy import *

# Constants
parameter_list = symbols("m g r t I J")
m, g, r, t, I, J = parameter_list

I = m*r**2/4  # Central moment of inertia about any diameter
J = m*r**2/2  # Central moment of inertia about normal axis

# Declare generalized coordinates and generalized speeds
(q1, q2, q3, q4, q5), q_list, qdot_list = gcs('q', 5, list=True)
(u1, u2, u3), u_list, udot_list = gcs('u', 3, list=True)

# List for holding all motion equations
eoms = []

# Create a Newtonian reference frame
N = NewtonianReferenceFrame('N')

# Assign all lists to N book keeping purposes
N.setcoords(q_list, qdot_list, u_list, udot_list)

# Intermediate reference frames
A = N.rotate("A", 3, q1)
B = A.rotate("B", 1, q2)

# Frame fixed to the torus rigid body.
C = B.rotate("C", 2, q3, I=(I, J, I, 0, 0, 0), I_frame=B)

# Locate the mass center of torus
CO = N.O.locate('CO', -r*B[3], frame=C, mass=m)

# Fixed inertial reference point
N1 = CO.locate('N1', r*B[3] - q4*N[1] - q5*N[2])

# Define the generalized speeds to be the B frame measure numbers of the angular
u_rhs = [dot(C.ang_vel(), B[i]) for i in (1, 2, 3)]


# Create the equations that define the generalized speeds, then solve them for
# the time derivatives of the generalized coordinates
u_definitions = [u - u_r for u, u_r in zip(u_list, u_rhs)]
kindiffs = solve(u_definitions, qdot_list)
print 'Kinematic differential equations'

for qd in qdot_list[:3]:
    kindiffs[qd] = expand(kindiffs[qd])
    print qd, '=', kindiffs[qd]
    eoms.append(kindiffs[qd])

# Form the expressions for q1' and q2', taken to be dependent speeds
# Completely optional, these have no influence on the dynamics.
nh = [dot(N1.vel(), N[1]), dot(N1.vel(), N[2])]
dependent_rates = solve(nh, q4.diff(t), q5.diff(t))
print 'Dependent rates:'
for qd in dependent_rates:
    dependent_rates[qd] = expand(dependent_rates[qd].subs(kindiffs))
    print qd, '=', dependent_rates[qd]
    eoms.append(dependent_rates[qd])

# Substitute the kinematic differential equations into velocity expressions,
# form partial angular velocities and partial velocites, form angular
# accelerations and accelerations
N.setkindiffs(kindiffs, dependent_rates)

# Apply gravity
N.gravity(g*A[3])

# Form Kane's equations and solve them for the udots
kanes_eqns = N.form_kanes_equations()
dyndiffs = solve(kanes_eqns, udot_list)

print 'Dynamic differential equations'
for r, ud in enumerate(udot_list):
    if r == 2:
        dyndiffs[ud] = collect(expand(dyndiffs[ud]), u1)
    else:
        dyndiffs[ud] = expand(dyndiffs[ud])
    print ud, '=', dyndiffs[ud]
    eoms.append(dyndiffs[ud])

