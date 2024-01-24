export OMP_NUM_THREADS=1
./solver_static static 2048 1000
./solver_static1 static1 2048 1000
./solver_dynamic dynamic 2048 1000

export OMP_NUM_THREADS=2
./solver_static static 2048 1000
./solver_static1 static1 2048 1000
./solver_dynamic dynamic 2048 1000

export OMP_NUM_THREADS=4
./solver_static static 2048 1000
./solver_static1 static1 2048 1000
./solver_dynamic dynamic 2048 1000

export OMP_NUM_THREADS=8
./solver_static static 2048 1000
./solver_static1 static1 2048 1000
./solver_dynamic dynamic 2048 1000

export OMP_NUM_THREADS=10
./solver_static static 2048 1000
./solver_static1 static1 2048 1000
./solver_dynamic dynamic 2048 1000

export OMP_NUM_THREADS=12
./solver_static static 2048 1000
./solver_static1 static1 2048 1000
./solver_dynamic dynamic 2048 1000