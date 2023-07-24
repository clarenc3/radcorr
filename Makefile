# Shared library
LIB=libRadCorr.so
# Executables
EXE=Calculate.exe CheckQ2.exe
# Input file
SRC=RadCorrCalculator.cpp
OBJ=$(subst .cpp,.o,$(SRC))

CXX=g++
CXXFLAGS=-g -O3 -Wall -Wextra -Wl,-rpath=./
LDFLAGS=-shared

ROOT_INCLUDES=$(shell root-config --cflags)
ROOT_LIBS=$(shell root-config --libs)
INCLUDES=$(ROOT_INCLUDES) -I.

all: $(LIB) $(EXE)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -fpic -c $^ $(ROOT_LIBS)

$(LIB): $(OBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

%.exe: %.cpp $(LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(ROOT_LIBS) 

clean:
	rm -f *.exe *.so *.o
