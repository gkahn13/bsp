all:
	gcc -c -std=c99 -o build/beliefMPC.o beliefMPC.c	
	g++ -c -std=c++0x -o build/timer.o timer/Timer.cpp
	g++ -I/usr/include/python2.7 -c -std=c++0x -o build/point-belief.o point-belief.cpp
	g++ -o bin/point build/beliefMPC.o build/point-belief.o build/timer.o -lpython2.7 -lboost_python
	
state:
	gcc -c -std=gnu99 -o build/stateMPC.o stateMPC.c	
	g++ -c -std=c++0x -o build/timer.o timer/Timer.cpp
	gcc -c -std=c99 -o build/symeval.o symeval.c
	g++ -I/usr/include/python2.7 -c -std=c++0x -o build/point-state.o point-state.cpp
	g++ -o bin/point-state build/stateMPC.o build/point-state.o build/timer.o build/symeval.o -lpython2.7 -lboost_python
	

clean:
	rm -f bin/*
	rm -rf build/*
