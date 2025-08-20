#!/bin/bash
 
#SBATCH --time=24:00:00 # walltime
#SBATCH --ntasks-per-node=112 # number of processor cores (i.e. tasks)
#SBATCH --nodes=1 # number of nodes
#SBATCH --wckey moose # Project Code
#SBATCH -J "test_run" # job name
#SBATCH --mem-per-cpu=2000M # 224M? 224G?
 
module load use.moose moose-dev-openmpi
moose-dev-shell
mpiexec ../../heat_transfer-opt -i ./fine_mesh_test2.i