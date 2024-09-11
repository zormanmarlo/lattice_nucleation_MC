# Lattice Monte Carlo Model for Simulation of Nucleation

## Overview: 
This model is written in Fortran and utilizes a particle representation of ions. Coordinates are on a lattice, spacing between lattice points is user defined. We use a 5:1 ratio of translation to AVBMC-2 moves. The model outputs energy, ccceptance statistics, a trajectory, as well as bath and target cluster data. 

Right now it is built around generating FES of nucleation. To that end, it includes options to bias the size of a target cluster. The target cluster is defined as the cluster that includes the particle with an index of 1 (i.e. the first particle in output xyz trajectories). However the code should still work for unbiased simulations and bath statistics. You will simply ignore target cluster outputs and set all biases to zero.

We have built it to be easily parallelizable. Each process (known as a markov chain) is completely independent and runs on its own core(s). 

## Use:
The code requires a number of input files:

flp.in -- Specifies simulation parameters such as box size, seed, total steps, temperature, output frequency, number of particles, etc.

f1 -- Interaction potentials / force field. The first three columns specify lattice indeces, the fourth specifies electrostatic interactions, and the last three specify Type_1 - Type_1, Type_1 - Type_2, and Type_2 - Type_2 interactions, respectively
bias.txt

Kn -- Specifies values to be used for linear biasing of a target cluster. The first line specifies the size of the largest cluster, the second line specifies the bias associated with a monomer, the third line specifies the bias associated with a dimer, and so on. Right now this file is required, and if the simulation sees a cluster larger than the largest size specified here it will crash. If you want to run unbiased simulations you should be able to just leave this as is.

input-XX.xyz -- Optional input file to start the simulation from. You will need to specify that you want to use an input structure in flp.in.

While running, the code outputs several files for each markov chain:

E-XX.out -- Total system energy and total bias energy
clusters-XX.out -- First element in each line specifies largest observed cluster at that step, subsequent elements specify number of observed monomers, dimers, etc
target-cluster-XX.out -- Provides the size of the cluster with the ion of index 1 in it (if you only care about bath statistics ignore this)
stats-XX.log -- Attempts and acceptances of the different move types

In order to run locally, you can use a command similar to the one in comp.sh. In order to run on klone, you can use a submission script similar to submit.sh. 

