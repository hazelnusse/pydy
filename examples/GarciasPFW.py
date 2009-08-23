# This program calculates the full nonlinear equations of motion (EOM) of
# Garcia's simplest point foot walker model.  The complete paper can be found:
# http://ruina.tam.cornell.edu/research/topics/locomotion_and_robotics/papers/some_results_passive_dynamic.pdf

from pydy import *

# Create a Newtonian reference frame
N = NewtonianReferenceFrame('N')

# Declare parameters, coordinates, speeds
# Define Model Parameters
# Mh = hip mass, (kg)
# mf = foot mass, (kg)
# g = gravitational acceleration constant, (m/s/s)
# L = leg length
# gamma = ground slope, (rad)
Mh, mf, g, L, q3 = N.declare_parameters('Mh mf g L q3')
(q1, q2), q_list, qdot_list = N.declare_coords('q', 2)
(u1, u2), u_list, udot_list = N.declare_speeds('u', 2)

# Orient reference frames
F = N.rotate('F', 1, q3)
A = F.rotate('A', 1, q1)
B = A.rotate('B', 1, -q2)

# Locate masses
P = N.O.locate('P', L*A[2], mass=Mh)
Q = P.locate('Q', -L*B[2], mass=mf)

# Define generalized speeds
u_defs = N.define_speeds(
        [Eq(u1, dot(A.ang_vel(), A[1])),
         Eq(u2, dot(B.ang_vel(), A[1]))])
print u_defs

# Compute the kinematic differential equations
T, Tinv, kindiffs = N.form_kindiffs(u_defs, qdot_list)

print 'Kinematic differential equations'
for qd in qdot_list:
    print qd, '=', kindiffs[qd]

# Set velocities to only involve expressions with u's
A._wrel = Vector(u1*A[1])
B._wrel = Vector((u2-u1)*A[1])
P._vrel = cross(A.ang_vel(N), P.rel(N.O))
Q._vrel = cross(B.ang_vel(N), Q.rel(P))

# Set the kinematic differential equations
N.setkindiffs(kindiffs)

# Apply gravity
N.gravity(-g*N[3])

# Form Kane's equations and solve them for the udots
kanes_eqns = N.form_kanes_equations()
dyndiffs = N.solve_kanes_equations()
print 'Dynamic differential equations'
for ud in udot_list:
    dyndiffs[ud] = dyndiffs[ud].expand()
    print ud, '=', dyndiffs[ud]

N.setdyndiffs(dyndiffs)
N.output_eoms('GarciesPFW_eoms.py')
