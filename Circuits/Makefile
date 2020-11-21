include ./CONFIG.mine

CC = g++
LIB = libMPC.a

CFLAGS = -Wall -std=c++11 -pedantic -Wextra -pthread -maes -mpclmul -msse4.1 -march=core2 $(FLAGS) $(OPT) -I$(OSSL)/include
CPPFLAGS = $(CFLAGS) 
LDLIBS = -L$(OSSL)/lib -lm -lssl -lcrypto -lmpirxx -lmpir -lcryptopp $(LDFLAGS)


CPPFLAGS = $(CFLAGS) 

all: 
	$(CC) $(CFLAGS) -o convert.x convert.cpp $(LIB) $(LDLIBS) -I.
	$(CC) $(CFLAGS) -o process.x process.cpp $(LIB) $(LDLIBS) -I.

clean:
	-rm *.x
	-cd VHDL; make clean

