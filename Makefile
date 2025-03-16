CXX = g++
CXXFLAGS = -std=c++11 -O2 -ftree-vectorize

#give actual names instead of using default symbols for ease
TARGET = md
UNITTESTS = unit_tests

SIM_SOURCES = serialSim.cpp
UNITTESTS_SOURCES = unit_tests.cpp

all: $(TARGET)

$(TARGET): $(SIM_SOURCES)
	$(CXX) $(CXXFLAGS) $(SIM_SOURCES) -o $(TARGET)

unittests: $(TARGET) $(UNITTESTS_SOURCES)
	$(CXX) $(CXXFLAGS) -o $(UNITTESTS) $(UNITTESTS_SOURCES)
	

.PHONY: unit_tests
unit_tests: unittests
	./$(UNITTESTS)

doc: Doxyfile
	doxygen Doxyfile

.PHONY: clean 

clean:
	rm -f $(TARGET) $(UNITTESTS) *.o *.txt



