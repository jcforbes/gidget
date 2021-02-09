LDFLAGS		= -O3  -L/usr/local/lib -lgsl -lgslcblas -lm 
CFLAGS		= -O3  -I/usr/local/include -c 
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

