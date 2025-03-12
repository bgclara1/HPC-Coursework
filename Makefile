compiler      = g++
flags  = -std=c++11
target   = simulation
sources  = simulation.cpp
objects  = $(SOURCES:.cpp=.o)

all: md

md: $(sources)
	$(compiler) $(flags) -o $(target) $(sources)

clean:
	rm -f $(target) $(objects) output.txt energy.txt positions.txt

.PHONY: all md clean

