all: solver_static solver_static1 solver_dynamic

solver_static: main_static.cpp
	g++ main_static.cpp -std=c++17 -O3 -march=native -ffast-math -fopenmp -o solver_static

solver_static1: main_static1.cpp
	g++ main_static1.cpp -std=c++17 -O3 -march=native -ffast-math -fopenmp -o solver_static1

solver_dynamic: main_dynamic.cpp
	g++ main_dynamic.cpp -std=c++17 -O3 -march=native -ffast-math -fopenmp -o solver_dynamic

clean:
	rm -f *.csv