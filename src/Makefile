LDFLAGS		=  -mmacosx-version-min=10.4 -O3 -g libgsl.a libgslcblas.a -static-libgcc -static-libstdc++
CFLAGS		=  -O3 -g -I../inc -c 
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

