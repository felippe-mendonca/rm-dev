# COMPILER = arm-linux-gnueabihf-g++
COMPILER = g++
FLAGS = -std=c++11 -O3 #-Wall -Werror -Wextra
MSGPACK  = $(shell pkg-config --libs --cflags msgpack)
ZEROMQ   = $(shell pkg-config --libs --cflags libzmq)
SO_DEPS = -lpthread -lAria -ldl -lm -lrt -larmadillo -lyaml-cpp $(MSGPACK) $(ZEROMQ) -I../rm/include -I/usr/loca/include -L/usr/local/lib

BINS = trajectory formation path final_position

all: $(BINS)

clean:
	rm -rf $(BINS)

trajectory:trajectory.cpp
	$(COMPILER) $^ -o $@ $(FLAGS) $(SO_DEPS)

formation:formation.cpp
	$(COMPILER) $^ -o $@ $(FLAGS) $(SO_DEPS)

path:path.cpp
	$(COMPILER) $^ -o $@ $(FLAGS) $(SO_DEPS)

final_position:final_position.cpp
	$(COMPILER) $^ -o $@ $(FLAGS) $(SO_DEPS)