benchmark : main.o benchmark.o
	g++ main.o benchmark.o -o benchmark

benchmark.o : Benchmark.cpp
	g++ -c Benchmark.cpp -o benchmark.o

main.o : main.cpp
	g++ -c main.cpp -o main.o