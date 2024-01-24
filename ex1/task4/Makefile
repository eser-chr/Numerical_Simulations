CC=g++
FLAGS=--std=c++17

.PHONY: all clean run

all: o_0 o_1 o_2 o_3

o_0: ex4.cpp
	$(CC) $(FLAGS) -O0 ex4.cpp -o flag_o0
o_1: ex4.cpp
	$(CC) $(FLAGS) -O1 ex4.cpp -o flag_o1
o_2: ex4.cpp
	$(CC) $(FLAGS) -O2 ex4.cpp -o flag_o2
o_3: ex4.cpp
	$(CC) $(FLAGS) -O3 ex4.cpp -o flag_o3

clean:
	rm  output.txt

run:
	./flag_o0
	./flag_o1
	./flag_o2
	./flag_o3
