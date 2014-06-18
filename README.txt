A Distributed Multidimensional Knapsack Problem Solver using Island Based
Particle Swarm Optimization.


DESCRIPTION:

This package contains an efficient CPU and GPU implementations of the platform
for solving various problems using the Particle Swarm Optimization (PSO) technique.
The Multidimensional Knapsack Problem (MKP) was chosen as a benchmark
for demonstration purposes, but it is one of the most widely used combinatorial
optimization problems. This solver can also be used in multinode systems,
such as HPC clusters, where one multi-core CPU is allocated to each swarm (process).
Running in a such configuration can lead to a better solution and/or
reduced time spent in searching for a solution.


--------------------------------------------------------------------------------


External libraries:
================================================================================
* OpenMPI
* papi


Build:
================================================================================
$ cd cpu
$ mkdir build
$ cd build
$ cmake ..
$ make


Example usage (4 PSO swarms):
================================================================================
$ mpirun -n 4 ./mkp --exch 1 --stat 0 --stop 300 --epochs 3000 --particles 128 --vel "5 10 15 20" --verbose ../../data/weing8.dat

