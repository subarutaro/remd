#ifndef H_FORCE_CALCULATOR
#define H_FORCE_CALCULATOR

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <string>

#include <iostream>
#include <fstream>

#include "vector.h"
//#include "molecules.h"
#include "variables.h"
#include "rotation.h"
#include "unit.h"
#include "timer.h"
#include "parameter.h"

//#define NO_PBC
//#define DEBUG


class ForceCalculator{
 public:
  //for vdw
  double **sgm,**eps;
  double rcut;
  //for ewald sum
  double *q;
  double kcut;
  double alpha;
  int    nwave;
  ivec3  *kvec;
  //for switching function
  double rswitch;
  //for confined
  int     confined_dim;
  double  wall_length;
  int     nfwall;
  dvec2 *fwall;
  double  sgm_wall;
  double  eps_wall;
  double  rho_wall;

  int nthreads;

  std::string prefix;
  //setting function
  void SetKvec();
  void CalcWallForce();

  void EwaldDirect(const int,const int,Molecule*,MolTypeList,Property&,const dvec3,double*);
  void EwaldWave(const int,const int,Molecule*,MolTypeList,Property&,const dvec3,double*);
  void EwaldSelf(const int,const int,Molecule*,MolTypeList,Property&,double*);
  void EwaldIntra(const int,const int,Molecule*,MolTypeList,Property&,const dvec3,double*);

 public:
  ForceCalculator(const double _rcut,const double _kmax,const double _alpha,const int _nthreads);
  ForceCalculator(const Parameter _p);
  ~ForceCalculator();
  void MakeParamList(MolTypeList);

  void LJ(Molecule*,MolTypeList,Property&,const int,const dvec3,double*);
  void Direct(Molecule*,const long,const dvec3);
  void Ewald(const int,const int,Molecule*,MolTypeList,Property&,const dvec3,double*);

  void SR(const Molecule*,Atom*,const dvec3,const int);
  void LR(Atom*,const dvec3,const int);
  void Confined(Molecule*,Atom*,const MolTypeList,const dvec3,const int,Property&);

  void Switching(Molecule*,Atom*,const MolTypeList,const dvec3,const int,Property&);
};

#endif
