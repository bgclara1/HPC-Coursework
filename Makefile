CXX = g++
CXXFLAGS = -std=c++11 -O3

NVCC = nvcc

MD = md
UNITTESTS = unit_tests
MDPAR  = mdpar
MDCUDA = mdcuda

SERIAL_SOURCES    = serialSim.cpp
UNITTESTS_SOURCES = unit_tests.cpp
PAR_SOURCES       = parallelSim.cpp
CUDA_SOURCES      = cuda.cu

all: $(MD) $(UNITTESTS) $(MDPAR) $(MDCUDA)

$(MD): $(SERIAL_SOURCES)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(UNITTESTS): $(UNITTESTS_SOURCES)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(MDPAR): $(PAR_SOURCES)
	$(CXX) $(CXXFLAGS) -fopenmp -o $@ $^

$(MDCUDA): $(CUDA_SOURCES)
	$(NVCC) $(CXXFLAGS) -o $@ $^

.PHONY: unittests
unittests: $(UNITTESTS)
	./$(UNITTESTS)

doc: Doxyfile
	doxygen Doxyfile

.PHONY: clean
clean:
	rm -f $(MD) $(UNITTESTS) $(MDPAR) $(MDCUDA) *.o *.txt
