#!/usr/bin/env python
from pydy import *
from sympy import solve, Eq

# Create a Newtonian reference frame
N = NewtonianReferenceFrame('N')

# Declare parameters, coordinates, speeds
m, g, l, b = params = N.declare_parameters('m g l b')
(q1,), (q1d,) = N.declare_coords('q', 1)
(u1,), (u1d,) = N.declare_speeds('u', 1)

# Orient the A reference frame
A = N.rotate('A', 1, q1)
# Define the position from the hinge
P = N.O.locate('P', -l*A[3], mass=m)

kd = Eq(q1d, u1)

# Set the angular velocity and velocities to only involve expressions involving
# the generalized speeds
A.abs_ang_vel = Vector(u1*A[1])
P.abs_vel = cross(A.ang_vel(N), P.rel(N.O))

A.abs_ang_acc = dt(A.abs_ang_vel, N)
P.abs_acc = dt(P.abs_vel, N)

# Apply gravity
N.gravity(-g*N[3])
# Apply a friction torque at hinge
A.apply_torque(-b*u1*A[1], N)

# Form Kane's equations and solve them for the udots
kanes_eqns = N.form_kanes_equations()
dd = Eq(u1d, solve(kanes_eqns, u1d)[u1d].expand())

# Energy calculations
ke = (m*P.abs_vel.mag_sqr)/2
pe = m*g*l*(1-cos(q1))
energy_eqs = [Eq(Symbol('ke'), ke), Eq(Symbol('pe'), pe)]

# Animation equations
anim_eqs = animate(N, ('P', P.rel(N.O)))

output_string = "from __future__ import division\n"
output_string += "from math import sin, cos\n\n"

ds = """\
Point mass pendulum equations of motion.

    _x is an array/list in the following order:
        q1:  Angle of pendulum link relative to vertical (0 downwards)
        u1:  A[1] measure number of the inertial angular velocity of the first link.

    _params is an array/list in the following order:
        m:  Mass of first pendulum point mass.
        l:  Length of first pendulum link.
        g:  Gravitational constant.
        b:  Damping coefficient at hinge.
"""

output_string += generate_function("eoms", [kd,dd], [q1, u1], params, docstring=ds,
        time=True)

ds = """\
Kinetic and Potential Energy of point mass pendulum.

    _x is an array/list in the following order:
        q1:  Angle of first pendulum link relative to vertical (0 downwards)
        u1:  A[1] measure number of the inertial angular velocity of the first link.

    _params is an array/list in the following order:
        m:  Mass of first pendulum point mass.
        l:  Length of first pendulum link.
        g:  Gravitational constant.
    Returns a list/array of kinetic energy and potential energy, respectively.
"""
output_string += generate_function("energy", energy_eqs, [q1, u1], params, docstring=ds)

ds = """\
Calculate configuration of pendulum for purposes of animation.

    _x is an array/list of the configuration coordinates of the disc:
        q1:  Angle of first pendulum link relative to vertical (0 downwards)
        u1:  A[1] measure number of the inertial angular velocity of the first link.

    _params is the radius of the disc.
        m:  Mass of first pendulum point mass.
        l:  Length of first pendulum link.
        g:  Gravitational constant.

    Output is:
          P:  Position of first point mass.
"""
output_string += generate_function("anim", anim_eqs, [q1], params, docstring=ds)

file = open('pendulum_lib.py', 'w')
file.write(output_string)
file.close()
