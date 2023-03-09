INCLUDES=-I../mixdclass
CXX=g++
CXXFLAGS=-g $(INCLUDES)

SRC=$(wildcard *.cpp)
OBJ := $(addsuffix .o, $(basename $(SRC)))
EXE=$(basename $(SRC))

all: $(EXE)

%.o: %cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(EXE): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	$(RM) $(OBJ)
	$(RM) $(EXE)

install: $(EXE)
	cp $< ~/local/bin

.PHONY: all clean

