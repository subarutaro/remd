CXX=g++
CXXFLAGS=-std=c++11 -Wall

SRC=resize.cxx
OBJ=resize.o
EXE=resize.out

$(EXE):	$(OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@
$(OBJ):	$(SRC)
	$(CXX) $(CXXFLAGS) $^ -c -o $@
clean:
	rm $(OBJ) $(EXE)
