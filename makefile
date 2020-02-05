# make           # compile all binary
# make clean     # remove ALL objects
# make clean_all # remove ALL objects and binary

.PHONY = all install clean rmproper

CC:=g++
EXE:=ss3d
SOURCE:=$(wildcard *.cpp) $(wildcard Algorithm/*.cpp) $(wildcard Module/*.cpp)
OBJ:=$(SOURCE:%.cpp=%.o)
FLAG:=-O2
LDLIBS:=-lm
LDFLAGS:=
STD:=-std=c++14

all: $(EXE)
	@echo "DBinary is in the latest version!"

$(EXE) : $(OBJ)
	@echo "DLinking Objects ... "
	$(CC) $^ -o $(EXE) $(LDFLAGS) $(LDLIBS)

%.o: %.cpp
	@echo "Compiling Source File... "
	$(CC) $< -o $@ $(FLAG) $(STD) -c

clean:
	@echo "Cleaning Object Files ..."
	rm -rvf $(OBJ)

clean_all: clean
	@echo "DCleaning Binary file ..."
	rm -rvf $(EXE)
