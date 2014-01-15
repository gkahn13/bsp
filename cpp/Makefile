OBJS_BELIEF = build/beliefMPC.o build/timer.o build/point-belief.o
OBJS_STATE = build/stateMPC.o build/symeval.o build/timer.o build/point-state.o
OBJS_SACHINSTATE = build/sachin-stateMPC.o build/sachin-symeval.o build/timer.o build/sachin-point-state.o
OBJS_LP = build/lpMPC.o build/symeval.o build/timer.o build/point-lp.o

build/timer.o : timer/Timer.h timer/Timer.cpp
	g++ -c -std=c++0x -o build/timer.o timer/Timer.cpp
	
build/symeval.o : symeval.h symeval.c
	gcc -c -std=gnu99 -o build/symeval.o symeval.c
	

# make point-belief
belief: $(OBJS_BELIEF)
	g++ $(OBJS_BELIEF) -o bin/point-belief -lpython2.7 -lboost_python
	
build/beliefMPC.o : beliefMPC.c beliefMPC.h
	gcc -c -std=c99 -o build/beliefMPC.o beliefMPC.c
	
build/point-belief.o : timer/Timer.h timer/Timer.cpp beliefMPC.c point-belief.cpp
	g++ -I/usr/include/python2.7 -c -std=c++0x -o build/point-belief.o point-belief.cpp

	
# make point-state
state: $(OBJS_STATE)
	g++ $(OBJS_STATE) -o bin/point-state -lpython2.7 -lboost_python
	
build/stateMPC.o : stateMPC.c stateMPC.h
	gcc -c -std=gnu99 -o build/stateMPC.o stateMPC.c
	
build/point-state.o : timer/Timer.h symeval.h timer/Timer.cpp symeval.c stateMPC.c point-state.cpp
	g++ -I/usr/include/python2.7 -c -std=c++0x -o build/point-state.o point-state.cpp
	

# make sachin-point-state
sachin-state: $(OBJS_SACHINSTATE)
	g++ $(OBJS_SACHINSTATE) -o bin/sachin-point-state
	
build/sachin-symeval.o : sachin-symeval.h sachin-symeval.c
	gcc -c -std=gnu99 -o build/sachin-symeval.o sachin-symeval.c

build/sachin-stateMPC.o : sachin-stateMPC.c sachin-stateMPC.h
	gcc -c -std=gnu99 -o build/sachin-stateMPC.o sachin-stateMPC.c
	
build/sachin-point-state.o : timer/Timer.h sachin-symeval.h timer/Timer.cpp sachin-symeval.c sachin-stateMPC.c sachin-point-state.cpp
	g++ -c -std=c++0x -o build/sachin-point-state.o sachin-point-state.cpp


# make point-lp
lp: $(OBJS_LP)
	g++ $(OBJS_LP) -o bin/point-lp -lpython2.7 -lboost_python
	
build/lpMPC.o : lpMPC.c lpMPC.h
	gcc -c -std=gnu99 -o build/lpMPC.o lpMPC.c
	
build/point-lp.o : timer/Timer.h symeval.h timer/Timer.cpp symeval.c lpMPC.c point-lp.cpp
	g++ -I/usr/include/python2.7 -c -std=c++0x -o build/point-lp.o point-lp.cpp


clean:
	rm -f bin/*
	rm -rf build/*
