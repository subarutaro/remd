OMPPPATH=$(HOME)/software/ompp/0.8.
OMPP=kinst-ompp
OMPPFLAGS=-I

CUDAPATH ?= /usr/local/cuda/6.5

CC = gcc
CXX = g++
NVCC = $(CUDAPATH)/bin/nvcc
#CC = kinst-ompp gcc
#CXX = kinst-ompp g++

#CFLAGS += -g -O0
CFLAGS += -O2
#CFLAGS += -O3
#CFLAGS +=  -Wall
CFLAGS += -m64
#CFLAGS += -DDEBUG

OMPFLAGS += -fopenmp
#CFLAGS +=  -lpthread -lgomp -lcudart

CXXFLAGS = -std=c++11

NVCCFLAGS  = --ptxas-options=-v
#for Fermi GPUs
#NVCCFLAGS += -gencode arch=compute_20,code=sm_20
#for Kepler 6xx series
#NVCCFLAGS += -gencode arch=compute_30,code=sm_30
#for Kepler titan and titan black
#NVCCFLAGS += -gencode arch=compute_35,code=sm_35
#for Maxwell 9xx
NVCCFLAGS += -gencode arch=compute_52,code=sm_52
#for Pascal 10xx
NVCCFLAGS += -gencode arch=compute_61,code=sm_61
NVCCFLAGS += -rdc=true

NVCCFLAGS += -Xcompiler $(OMPFLAGS)

#NVCCFLAGS += -G

PLUGIN_DIR = ../plugin

INC = -I$(CUDAPATH)/include  -I$(PLUGIN_DIR) -I../inc

#SRC=$(shell ls *.cxx)
#CUSRC=$(shell ls *.cu)
#HED=$(shell ls *.h)
#CUHED=$(shell ls *.cuh)
SRC=force_calculator.cxx integrator.cxx  main.cxx molecules.cxx remdinfo.cxx gro.cxx io_manager.cxx molecular_dynamics.cxx remd.cxx wham.cxx
CUSRC=remdgpu.cu
CUSRC+= bitonic_sort.cu  quicksort.cu sort.cu
HED=force_calculator.h io_manager.h molecules.h remd.h rotation.h unit.h wham.h gro.h mdunit.h parameter.h remdgpu.h timer.h variables.h integrator.h molecular_dynamics.h product.h remdinfo.h typedef.h vector.h
CUHED= reduce.cuh vector.cuh
CUHED+= 4site.cuh bitonic_sort.cuh quicksort.cuh scan.cuh  sort.cuh

DEP =$(SRC:.cxx=.d)
#DEP+=$(CUSRC:.cu=.d)

OBJ =$(SRC:.cxx=.o)
OBJ+=$(CUSRC:.cu=.o)

TARGET=md.out

.SUFFIXES:	.cxx .C .cu .cuh .h. .o
.cxx.o:	$(SRC) $(HED)
	$(CXX) $(CFLAGS) $(CXXFLAGS) $(OMPFLAGS) $(INC) -c -MMD -MP $< -o $@
.c.o:	$(SRC) $(HED)
	$(CC) $(CFLAGS) -c $< -o $@
.cu.o:	$(CUSRC) $(HED) $(CUHED)
	$(NVCC) $(NVCCFLAGS) $(CFLAGS) $(CXXFLAGS) $(INC) -c $< -o $@


 all:	$(TARGET)

-include $(DEP)

 $(TARGET):	$(OBJ)
#	$(CXX) $(CFLAGS) $(CXXFLAGS) -o $(TARGET) $(OBJ)
	$(NVCC) $(CFLAGS) $(CXXFLAGS) $(NVCCFLAGS) $^ -o $(TARGET)
test:	$(OBJ)
#	$(CXX) $(CFLAGS) $(CXXFLAGS) -o $(TARGET) $(OBJ)
	$(NVCC) $(CFLAGS) $(CXXFLAGS) $(NVCCFLAGS) -o test.out $(OBJ)
takaiwa:	$(OBJ)
	$(CXX) $(CFLAGS) $(CXXFLAGS) $(OMPFLAGS) $(INC) -c -MMD -MP -DTAKAIWA force_calculator.cxx -o force_calculator.o
	$(NVCC) $(CFLAGS) $(CXXFLAGS) $(NVCCFLAGS) $(INC) -DTAKAIWA -c remdgpu.cu -o remdgpu.o
	$(NVCC) $(CFLAGS) $(CXXFLAGS) $(NVCCFLAGS) -DTAKAIWA -o takaiwa.out $(OBJ)
	rm remdgpu.o
steepest:	$(OBJ)
	$(NVCC) $(CFLAGS) $(CXXFLAGS) $(NVCCFLAGS) $(INC) -DSTEEPEST -c remdgpu.cu -o remdgpu.o
	$(NVCC) $(CFLAGS) $(CXXFLAGS) $(NVCCFLAGS) -DSTEEPEST -o steepest.out $(OBJ)
	rm remdgpu.o
tuning:	$(OBJ)
	$(NVCC) $(CFLAGS) $(CXXFLAGS) $(NVCCFLAGS) $(INC) -DTUNING -DTIP3P -c remdgpu.cu -o remdgpu.o
	$(NVCC) $(CFLAGS) $(CXXFLAGS) $(NVCCFLAGS) -DTUNING -DTIP3P -o tip3p.out $(OBJ)
	$(NVCC) $(CFLAGS) $(CXXFLAGS) $(NVCCFLAGS) $(INC) -DTUNING -DTIP4P -c remdgpu.cu -o remdgpu.o
	$(NVCC) $(CFLAGS) $(CXXFLAGS) $(NVCCFLAGS) -DTUNING -DTIP4P -o tip4p.out $(OBJ)
	$(NVCC) $(CFLAGS) $(CXXFLAGS) $(NVCCFLAGS) $(INC) -DTUNING -DTIP4PICE -c remdgpu.cu -o remdgpu.o
	$(NVCC) $(CFLAGS) $(CXXFLAGS) $(NVCCFLAGS) -DTUNING -DTIP4PICE -o tip4pice.out $(OBJ)
	$(NVCC) $(CFLAGS) $(CXXFLAGS) $(NVCCFLAGS) $(INC) -DTUNING -DTIP5P -c remdgpu.cu -o remdgpu.o
	$(NVCC) $(CFLAGS) $(CXXFLAGS) $(NVCCFLAGS) -DTUNING -DTIP5P  -o tip5p.out $(OBJ)
clean:
	rm $(OBJ) $(TARGET) $(DEP)
