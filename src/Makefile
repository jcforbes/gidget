LDFLAGS		= -O3  -L/n/helmod/apps/centos7/Comp/gcc/7.1.0-fasrc01/gsl/2.4-fasrc01/lib64 -lgsl -lgslcblas -lm 
CFLAGS		= -O3  -I/n/helmod/apps/centos7/Comp/gcc/7.1.0-fasrc01/gsl/2.4-fasrc01/include -c 
CC		= g++
EXECUTABLE=../bin/gidget

SOURCES := $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE) 

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

$(EXECUTABLE2): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o $(EXECUTABLE) $(EXECUTABLE).dSYM

