# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -O2 -fPIC `root-config --cflags`
LDFLAGS = -shared `root-config --libs` -lRooFit -lRooFitCore

# Files
SRC = PhiBinnedFitRunner.cpp
HDR = PhiBinnedFitRunner.h
DICT = PhiBinnedFitRunnerDict.cxx
DICT_HDR = LinkDef.h

OBJ = PhiBinnedFitRunner.o PhiBinnedFitRunnerDict.o
LIB = libPhiBinnedFitRunner.so

# Default target
all: $(LIB)

# Generate dictionary source with rootcling
$(DICT): $(HDR) $(DICT_HDR)
	rootcling -f $(DICT) -c $(HDR) $(DICT_HDR)

# Compile object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
	
%.o: %.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Build shared library
$(LIB): $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@

# Clean generated files
clean:
	rm -f $(OBJ) $(DICT) $(DICT:.cxx=.pcm) $(LIB)

.PHONY: all clean
