CXX := clang++
CXXFLAGS := -std=c++17 -Wall -pedantic -O3 -ffast-math
INCLUDES := -I/usr/include/eigen3

TARGET := cg

all: $(TARGET)

$(TARGET): main.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $<

.PHONY: clean

clean:
	rm -f $(TARGET)
