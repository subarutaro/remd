#ifndef H_GRO
#define H_GRO

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>

#include "vector.h"

class AtomGro{
public:
/*
  std::string AtomName;
  std::string MolName;
//*/
  char  AtomName[6];
  char  MolName[6];
  int   AtomID;
  int   MolID;
  dvec3 coor;
  dvec3 vel;

  friend std::ostream &operator<<(std::ostream &s, AtomGro a){//output atom info in gro file format
    s << std::setw(5) << a.MolID
      << std::setw(5) << a.MolName
      << std::setw(5) << a.AtomName
      << std::setw(5) << a.AtomID;
    s << std::fixed
      << std::setprecision(3);
    for(int d=0;d<3;d++){
      s << std::setw(8) << a.coor[d];
    }
    for(int d=0;d<3;d++){
      s << std::setw(8) << a.vel[d];
    }
    return s;
  };

};

void read_gro_natom(char*,int&);
void read_gro_atom(char*,AtomGro*,int);
void read_gro_cellsize(char*r,dvec3&);

void read_gro_file(char*,int&,AtomGro*,dvec3&,int&);
void write_gro_file(const char*,const int,const AtomGro*,const dvec3,const char *);

#endif
