INCLUDES=-I../mixdclass
CXX=g++
CXXFLAGS=-g $(INCLUDES)

SRC=gmsh2mixd_3D.cpp
OBJ := $(addsuffix .o, $(basename $(SRC)))
EXE=$(basename $(SRC))

all: $(EXE)

%.o: %cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(EXE): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^
	cp $@ gmsh2mixd

clean:
	$(RM) $(OBJ)
	$(RM) $(EXE)

install: $(EXE)
	cp $< ~/local/bin/

.PHONY: all clean

