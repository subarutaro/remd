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
#include "profiler.h"

//#define NO_PBC
//#define DEBUG
#ifdef USE_FP32
#define FP float
#else
#define FP double
#endif
//#define FP float

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

  Profiler prof;

  int nthreads;

#ifdef SWITCHING
  __declspec(align(64)) FP *gx = nullptr;
  __declspec(align(64)) FP *gy = nullptr;
  __declspec(align(64)) FP *gz = nullptr;
  __declspec(align(64)) FP *gfx = nullptr;
  __declspec(align(64)) FP *gfy = nullptr;
  __declspec(align(64)) FP *gfz = nullptr;
  __declspec(align(64)) FP *gvx = nullptr;
  __declspec(align(64)) FP *gvy = nullptr;
  __declspec(align(64)) FP *gvz = nullptr;
  __declspec(align(64)) FP *glj = nullptr;
  __declspec(align(64)) FP *gcl = nullptr;

  __declspec(align(64)) FP *ax = nullptr;
  __declspec(align(64)) FP *ay = nullptr;
  __declspec(align(64)) FP *az = nullptr;
  __declspec(align(64)) FP *afx = nullptr;
  __declspec(align(64)) FP *afy = nullptr;
  __declspec(align(64)) FP *afz = nullptr;
  __declspec(align(64)) FP *as = nullptr;
  __declspec(align(64)) FP *ae = nullptr;
  __declspec(align(64)) FP *aq = nullptr;
#endif

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
  template<int>
  void SwitchingTuning(Molecule*,Atom*,const MolTypeList,const dvec3,const int,Property&,const int*,const int *,const int);
  void Switching3site (Molecule*,Atom*,const MolTypeList,const dvec3,const int,Property&,const int*,const int *,const int);
  void Switching4site (Molecule*,Atom*,const MolTypeList,const dvec3,const int,Property&,const int*,const int *,const int);
  void Switching5site (Molecule*,Atom*,const MolTypeList,const dvec3,const int,Property&,const int*,const int *,const int);
};

#endif
