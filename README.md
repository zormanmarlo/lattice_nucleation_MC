# Lattice Monte Carlo Model for Simulation of Nucleation

## Overview: 
This model is written in Fortran and utilizes a particle representation of ions. Coordinates are on a lattice, spacing between lattice points is user defined. We use a 5:1 ratio of translation to AVBMC-2 moves. AVBMC-2 out->in moves utilize Rosenbluth sampling. The model outputs energy, ccceptance statistics, a trajectory, as well as bath and target cluster data.

Right now it is built around generating FES of nucleation. To that end, it includes options to bias the size of a target cluster. The target cluster is defined as the cluster that includes the particle with an index of 1 (i.e. the first particle in output xyz trajectories). However the code should still work for unbiased simulations and bath statistics. You will simply ignore target cluster outputs and set all biases to zero.

We have built it to be easily parallelizable. Each process (known as a markov chain) is completely independent and runs on its own core(s). 

## Use:
The code requires a number of input files:

* flp.in -- Specifies simulation parameters such as box size, seed, total steps, temperature, output frequency, number of particles, etc. It should take the following form:
  
_(Lattice size) (seed) (num production steps) (kT) (num equil steps) (stats output freq) (xyz output freq)_

_(num particle type 1) (num particle type 2)_

_(charge type 1) (charge type 2) (electrostatic scalar) (kT for use with electrostatics)_
 
_(bias target - ignore) (number of rosenbluth sites)_

_(1 = use input, 0 = generate initial configuration)_

* f1 -- Interaction potentials / force field. Note that since we explicit Ewald electrostatics, the particle-specific interactions will need to have electrostatics subtracted out. The columns specify:
  
_(x) (y) (z) (electrostatic potential) (Type_1-Type_1) (Type_1-Type_2) (Type_2-Type_2)_

* Kn -- Specifies values to be used for linear biasing of a target cluster. The first line specifies the size of the largest cluster, the second line specifies the bias associated with a monomer, the third line specifies the bias associated with a dimer, and so on. Right now this file is required, and if the simulation sees a cluster larger than the largest size specified here it will crash. If you want to run unbiased simulations you should be able to just leave this as is.
* input-XX.xyz -- Optional input file to start the simulation from. You will need to specify that you want to use an input structure in flp.in.

While running, the code outputs several files for each markov chain:

E-XX.out -- Total system energy and total bias energy
clusters-XX.out -- First element in each line specifies largest observed cluster at that step, subsequent elements specify number of observed monomers, dimers, etc
target-cluster-XX.out -- Provides the size of the cluster with the ion of index 1 in it (if you only care about bath statistics ignore this)
stats-XX.log -- Attempts and acceptances of the different move types

In order to run locally, you can use a command similar to the one in comp.sh. In order to run on klone, you can use a submission script similar to submit.sh.

## Notes:
Here are a few things to consider:
* The lattice spacing is set during generation of f1. For now, both f1 files you have use a spacing of 1 Angstrom. If you want / need a larger spacing you will need to generate a new f1 file and change the electrostatic scalar in flp.in accordingly (will write instructions for this out as it becomes necessary)
* You have two f1 files (f1-zncl2-cl, f1-zncl2-dft). The code will only use an f1 file with the exact name f1. So you'll need to copy/rename the one you want to use. 
* The maximum box size you can run right now is Nl = 200 (i.e. a 200x200x200 A^3 box). Note that the f1 file always goes to HALF the full lattice size.
* If you want to run smaller boxes, you'll have to trim the f1 file. I have a python script that does this. It's used like this, where input file is the original f1 file and Nl is HALF the lattice size you want:

_python3 trim_f1.py INPUT_FILE Nl

