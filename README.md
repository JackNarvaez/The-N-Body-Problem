# MPI-N-Body-Problem

Authors: Ochoa J. and Narváez J.

Simulation of a system of N bodies interacting with each other gravitationally. This code uses the parallel programming library: Message Passing Interface (`MPI`),
it focuses on the implementation of `ring code`.

The evolution step is performed using `Euler's method`, a first-order algorithm. By default, the code runs with the maximum number of cores.
