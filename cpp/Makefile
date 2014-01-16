OBJS_BELIEF = build/beliefMPC.o build/timer.o build/point-belief.o build/logging.o
OBJS_STATE = build/stateMPC.o build/symeval.o build/timer.o build/point-state.o build/logging.o
OBJS_LP = build/lpMPC.o build/symeval.o build/timer.o build/point-lp.o build/logging.o

build/timer.o : timer/Timer.cpp
	g++ -c -std=c++0x -o build/timer.o timer/Timer.cpp

build/logging.o : logging.cpp
	g++ -c -std=c++0x -o build/logging.o logging.cpp
	
build/symeval.o : symeval.c
	gcc -c -std=gnu99 -o build/symeval.o symeval.c

# make point-belief
belief: $(OBJS_BELIEF)
	g++ $(OBJS_BELIEF) -o bin/point-belief -lpython2.7 -lboost_python
	
build/beliefMPC.o : beliefMPC.c
	gcc -c -std=c99 -o build/beliefMPC.o beliefMPC.c
	
build/point-belief.o : point-belief.cpp
	g++ -I/usr/include/python2.7 -c -std=c++0x -o build/point-belief.o point-belief.cpp

	
# make point-state
state: $(OBJS_STATE)
	g++ $(OBJS_STATE) -o bin/point-state -lpython2.7 -lboost_python
	
build/stateMPC.o : stateMPC.c
	gcc -c -std=gnu99 -o build/stateMPC.o stateMPC.c
	
build/point-state.o : point-state.cpp
	g++ -I/usr/include/python2.7 -c -std=c++0x -o build/point-state.o point-state.cpp
	

# make point-lp
lp: $(OBJS_LP)
	g++ $(OBJS_LP) -o bin/point-lp -lpython2.7 -lboost_python
	
build/lpMPC.o : lpMPC.c
	gcc -c -std=gnu99 -o build/lpMPC.o lpMPC.c
	
build/point-lp.o : point-lp.cpp
	g++ -I/usr/include/python2.7 -c -std=c++0x -o build/point-lp.o point-lp.cpp


clean:
	rm -f bin/*
	rm -rf build/*
	rm -f *~
