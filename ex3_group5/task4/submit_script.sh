#!/bin/bash
# Number of cores
#SBATCH -c 40
# Runtime of this jobs is less then 10 minutes
#            (hh:mm:ss)
#SBATCH --time=00:10:00
# Clear the environment 
module purge > /dev/null 2>&1

export OMP_NUM_THREADS=1
./solver_static static 2048 10000
./solver_static1 static1 2048 10000
./solver_dynamic dynamic 2048 10000

export OMP_NUM_THREADS=2
./solver_static static 2048 10000
./solver_static1 static1 2048 10000
./solver_dynamic dynamic 2048 10000

export OMP_NUM_THREADS=4
./solver_static static 2048 10000
./solver_static1 static1 2048 10000
./solver_dynamic dynamic 2048 10000

export OMP_NUM_THREADS=8
./solver_static static 2048 10000
./solver_static1 static1 2048 10000
./solver_dynamic dynamic 2048 10000

export OMP_NUM_THREADS=10
./solver_static static 2048 10000
./solver_static1 static1 2048 10000
./solver_dynamic dynamic 2048 10000

export OMP_NUM_THREADS=20
./solver_static static 2048 10000
./solver_static1 static1 2048 10000
./solver_dynamic dynamic 2048 10000

export OMP_NUM_THREADS=40
./solver_static static 2048 10000
./solver_static1 static1 2048 10000
./solver_dynamic dynamic 2048 10000
