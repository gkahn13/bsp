all:
	gcc -c -std=c99 -o build/beliefMPC.o beliefMPC.c	
	g++ -c -std=c++0x -o build/timer.o timer/Timer.cpp
	g++ -I/usr/include/python2.7 -c -std=c++0x -o build/point-belief.o point-belief.cpp
	g++ -o bin/point build/beliefMPC.o build/point-belief.o build/timer.o -lpython2.7 -lboost_python

clean:
	rm -f bin/*
	rm -rf build/*
