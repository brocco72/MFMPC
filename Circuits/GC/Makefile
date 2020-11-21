OSSL = /usr/local/opt/openssl

MAX_MOD = 7
FLAGS = -DMAX_MOD_SZ=$(MAX_MOD)
OPT = -O3

CC = g++

CFLAGS = -Wall -std=c++11 -pedantic -Wextra -pthread -I$(ROOT)/src -maes -mpclmul -msse4.1 -mavx -march=core2 $(FLAGS) $(OPT) -I$(OSSL)/include
CPPFLAGS = $(CFLAGS)
LDLIBS = -L/$(OSSL)/lib -lm -lssl -lcrypto -lmpirxx -lmpir -lcryptopp $(LDFLAGS)

GC = $(patsubst %.cpp,%.o,$(wildcard ./*.cpp))
EXCEPTIONS = $(patsubst %.cpp,%.o,$(wildcard ../Exceptions/*.cpp))

LIB = libMPC.a

all: $(LIB)

%.o: %.cpp
	$(CC) $(CFLAGS) -o $@ -c $<

$(LIB): $(GC) $(EXCEPTIONS)
	$(AR) -crs libMPC.a $^
	- mv libMPC.a ../

clean:
	- rm *.o
	- cd ../Exceptions ; rm *.o
	- cd .. ; rm libMPC.a

