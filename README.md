A Distributed Multidimensional Knapsack Problem Solver using Island Based Particle Swarm Optimization
=====================================================================================================

This package contains an efficient CPU and GPU implementations of the platform
for solving various problems using the Particle Swarm Optimization (PSO)
technique.  The Multidimensional Knapsack Problem (MKP) was chosen as
a benchmark for demonstration purposes, but it is one of the most widely
used combinatorial optimization problems. This solver can also be used in
parallel environment such as HPC cluster, where one or more multi-core CPUs
or GPUs are allocated for each swarm (process). Running in a such
configuration can lead to a better solution and/or reduced time spent
in searching for a solution.



### External libraries:

* OpenMPI
* papi


### Build:

> $ cd cpu
> $ mkdir build
> $ cd build
> $ cmake ..
> $ make


### Example (4 swarms):

> $ ./mkp --help
> $ mpirun -np 4 ./mkp --exch 1 --stat 0 --stop 300 --epochs 3000 --particles 128 --vel "5 10 15 20" --verbose ../../data/weing8.dat

