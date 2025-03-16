CXX = g++
CXXFLAGS = -std=c++11 

#give actual names instead of using default symbols for ease
TARGET = md
UNITTESTS = unit_tests

SIM_SOURCES = serialSim.cpp
UNITTESTS_SOURCES = unit_tests.cpp

all: $(TARGET)

$(TARGET): $(SIM_SOURCES)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SIM_SOURCES)

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



