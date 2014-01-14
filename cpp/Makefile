OBJS = build/stateMPC.o build/symeval.o build/timer.o build/point-state.o
FLAGS_C = -c -std=gnu99
BUILD_DIR = build
	
state: $(OBJS)
	g++ $(OBJS) -o bin/point-state
	
build/stateMPC.o : stateMPC.c
	gcc -c -std=gnu99 -o build/stateMPC.o stateMPC.c
	
build/timer.o : timer/Timer.h timer/Timer.cpp
	g++ -c -std=c++0x -o build/timer.o timer/Timer.cpp
	
build/symeval.o : symeval.h symeval.c
	gcc -c -std=gnu99 -o build/symeval.o symeval.c
	
build/point-state.o : timer/Timer.h symeval.h timer/Timer.cpp symeval.c stateMPC.c point-state.cpp
	g++ -I/usr/include/python2.7 -c -std=c++0x -o build/point-state.o point-state.cpp

belief:
	gcc -c -std=c99 -o build/beliefMPC.o beliefMPC.c	
	g++ -c -std=c++0x -o build/timer.o timer/Timer.cpp
	g++ -I/usr/include/python2.7 -c -std=c++0x -o build/point-belief.o point-belief.cpp
	g++ -o bin/point build/beliefMPC.o build/point-belief.o build/timer.o -lpython2.7 -lboost_python


clean:
	rm -f bin/*
	rm -rf build/*
