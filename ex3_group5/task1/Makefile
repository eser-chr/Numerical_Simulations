CXX := clang++
CXXFLAGS := -std=c++14 -O3 -march=native -ffast-math

cg: main.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

.PHONY: clean

clean:
	rm -f cg
