LDFLAGS		=  -mmacosx-version-min=10.4 -O3 libgsl.a libgslcblas.a -static-libgcc -static-libstdc++ 
CFLAGS		= -O3 -I../inc -c
CC		= g++ 
EXECUTABLE=../bin/gidget

SOURCES := $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o $(EXECUTABLE)
