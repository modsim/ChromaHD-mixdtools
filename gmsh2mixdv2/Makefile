# INCLUDES=-I../mixdclass
# CXX=g++-9.2
CXX=g++
CXXFLAGS=-O3 -fopenmp

SRC=gmsh2mixdv2.cpp endian.cpp
OBJ := $(addsuffix .o, $(basename $(SRC)))
EXE=gmsh2mixdv2

all: $(EXE)

%.o: %cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

$(EXE): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ 
	# cp $@ gmsh2mixd.v2

clean:
	$(RM) $(OBJ)
	$(RM) $(EXE)
	# $(RM) gmsh2mixd.v2

install: $(EXE)
	cp $< ~/local/bin/

.PHONY: all clean

