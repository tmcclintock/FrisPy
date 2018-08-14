# FrisPy
This is a frisbee flight simulator written (mostly) in python.

# Purpose
In order to predict how a frisbee is going to fly when thrown, 
in addition to knowing the initial conditions of the throw, 
such as the position and velocity, one must also know a small 
set of "parameters" that describe how the forces and torques scale.
There are twelve parameters in total, with 1-3 being needed to describe
each force/torque.

Obviously, the forces are supposed to cause torques, however 
in order to account for this then one would have to calculate the
center of pressure - the location that aerodynamic forces actually
act on a rigid body. Doing such a calculation is costly and very challenging,
requiring knowledge of hydrodynamics. Instead, we assume that the
torques and forces are decoupled and independent, albiet with similar
dependencies on the kinematic variables. Thus, we introduce extra parameters,
known in the literature and in the program as "flight coefficients" or
usually just "coefficients".

This decoupling has been used for a few decades, most recently in the 
Masters thesis by Sarah Hummel (2003) at UCSD. Her program was one large
Matlab program that did a calculation nearly identical to the one written here.
The motivation for writing Frispy was twofold: to write it in more "accessible"
languages, and to make it available. To elaborate, Matlab is closed source and in
general not really used outside of academic settings. On the other hand, C
and Python are very prevalent. For this project, most of the math intensive
parts are written in C, while everything else is done in Python. This lets
the Markov Chain Monte-Carlo (MCMC) sampler I use (emcee) run extremely
quickly.

In terms of availability, Hummel's code only existed in written form in the
appendix of her thesis. I have chosen to host my code on GitHub, so that in
the future as people collect data it can easily be added to the repository
for anyone to analyze.

# Modes of Operation
[In progress] As of now, there are two modes of operation, visual mode and 
flight analysis mode. In the former, the program gives a visual display of the
trajectory of a disc for a specified set of initial conditions and an
assumed set of coefficients. In the latter, the program calculates the
best set of coefficients given a set of telemetric flight data provided
by the user.

In order to run in animation mode you must first compile the C source
code that does all the heavy lifting. From the home directory:

cd src
make
cd ..

This makes the shared library src/rk4.so and takes you back to FrisPy/.
Then, in order to actually run in animation mode:

python FrisPy.py animation config_file

The configuration file config_file contains the information

In order to run in flight analysis mode, one would run:

IN PROGRESS

# Pure-Python Implementation
A slower, pure python implementation of this code has been written, but it takes
too long to return results on some machines, so I have decided to rewrite the 
math intensive parts in C. The pure python version can be found in the subdirectory
Old_Python_Animation/. IN PROGRESS