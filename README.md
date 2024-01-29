# NBody-Marzagora-Tortorella-Vignolo
Implementation of "N-Body Problem" for AMSC course.

# Overview
This program implements two solvers for the N-Body Problem simulation for elastic collisions between particles of the same force (gravitational or Coulomb force): the first one is the standard algorithm, which takes advantage of the simmetry of the computation, Netwon's laws and the Euler's method to compute velocities and positions. The second solver, instead, is based on the Barnes and Hut algorithm to speedup the computation using a tree-like structure, a quadtree, to identify the position of the particles inside the simulation area and to take advantage of the center of mass of particles to calculate forces between particles that are far from the others. In this algorithm, we used the Verlet approximation to compute the new positions and velocities in order to obtain a better estimate during the calculations and a parallel write on file to accelerate the algorithm even more.

The simulation can be performed in 2D or 3D for the standard solver and in 2D for Barnes and Hut, using the same code and the result is displayed accordingly.
Moreover, the program runs the simulation both in parallel, implemented with OpenMP, and serially. 

# How to run the visual simulation 

The executable is run autonomously directly by python file that displays the resulting animation.

You can find it inside the graphics folder:

`src/graphics`

The program can then be executed using the command

`python3 animation.py`

or, with the same command, adding the version of python installed on the machine

(example with python 3.10)
`python3.10 animation.py`

The program will ask the user if they want to run the default version by pressing 'd'. If any other key is pressed, the program shows the helper list with 

all the parameters that can be changed and the program terminates.

You then have to run it again, adding the parameters you want to change between two "".

For example, 
`python3 animation.py "-dim 2 -simT 0"`

Default values are set on serial 2D simulation, gravitational force.


# How to run ONLY the computational part

The computational part of the code can be found inside the folder

`src/main`

It can be run alone using the following command

`g++ -fopenmp -I../utils main.cpp -o main`

The flag `-fopenmp` is needed to run the parallel code, while it will be ignored if the user intends only to run the serial simulation

The flag `-I../utils` specifies the program the path to the headers particle.hpp and force.hpp where all the information about particles and forces are implemented.

To execute it, run `./main ` 

in this way, the program will ask the user to press d to run the default function; if the user presses any other button, the helper function is shown and the program terminates. 

To change the parameters of the simulation run 

`./main -flag <value>`

where flags and the values needed are specificied inside helper function displayed by the program.

If any mistakes is made, the program terminates and it must be run again to perform the simulation.

# Documentation

You can find a detailed report about the reasoning behind the algorithms inside the "documentation folder", as well as a reference manual with all the code implemented for this project and a web page generated usind Doxygen.