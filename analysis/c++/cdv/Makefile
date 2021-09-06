CXX=g++
CXXFLAGS=-std=c++11 -O2 -m64

SRC=$(shell ls *.cxx)
HED=$(shell ls *.h)
OBJ=$(SRC:.cxx=.o)

.SUFFIXES: .cxx .h. .o

.cxx.o:	$(SRC) $(HED)
	$(CXX) $(CXXFLAGS) -c $< -o $@

all:	$(OBJ)
	$(CXX) $(CXXFLAGS) $^ -o define_inner.out
	$(CXX) $(CXXFLAGS) -DOPEN -c main.cxx -o main.o
	$(CXX) $(CXXFLAGS) $^ -o open_cylinder.out

clean:
	rm $(OBJ) $(TARGET)
