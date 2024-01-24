output: main.o info.o
	g++ main.o info.o -o output

main.o: main.cpp
	g++ -c main.cpp

info.o : info.cpp info.hpp
	g++ -c info.cpp

clean : 
	rm -f *.o

cleanall :
	rm -f *.o *.exe output.txt