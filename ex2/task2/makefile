output: main.o
	g++ main.o -std=c++17 -O3 -lopenblas -DEIGEN_DONT_PARALLELIZE -march=native -ffast-math -o mmm

main.o: main.cpp
	g++ -c main.cpp -std=c++17 -O3 -lopenblas -DEIGEN_DONT_PARALLELIZE -march=native -ffast-math

clean:
	rm *.o mmm
