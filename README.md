# NBody-Marzagora-Tortorella-Vignolo
Implementation of "N-Body Problem" for AMSC course.

# Overview
This program implements the N-Body Problem simulation for elastic collisions between particles of the same force (gravitational or Coulomb force) with Euler's method in C++.
The simulation can be performed in 2D or 3D using the same code and the result is displayed accordingly. 
Moreover, the program runs the simulation both in parallel, implemented with OpenMP, and serially. 

# How to run

The executable is run autonomously directly by python file that displays the resulting animation.
You can find it inside the graphics folder:
`NBody-Marzagora-Tortorella-Vignolo/src/graphics`
The program can then be executed using the command 
`python3 animation.py`
or, with the same command, adding the version of python installed on the machine
(example with python 3.10)
`python3.10 animation.py`