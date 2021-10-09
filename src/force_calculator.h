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
#define FPvec fvec3
#else
#define FP double
#define FPvec dvec3
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

  // tmp for calc force
  bool first_call = true;
  FP* lj_omp;
  FP* clmb_omp;
  FPvec*  vir_omp;
  FP* wall_omp;

#ifdef SWITCHING
#ifdef __INTEL_COMPILER
#define ATTR_ALIGN __declspec(align(64))
#else
#define ATTR_ALIGN
#endif
  ATTR_ALIGN FP *gx  = nullptr;
  ATTR_ALIGN FP *gy  = nullptr;
  ATTR_ALIGN FP *gz  = nullptr;
  ATTR_ALIGN FP *gfx = nullptr;
  ATTR_ALIGN FP *gfy = nullptr;
  ATTR_ALIGN FP *gfz = nullptr;
  ATTR_ALIGN FP *gvx = nullptr;
  ATTR_ALIGN FP *gvy = nullptr;
  ATTR_ALIGN FP *gvz = nullptr;
  ATTR_ALIGN FP *glj = nullptr;
  ATTR_ALIGN FP *gcl = nullptr;

  ATTR_ALIGN FP *ax  = nullptr;
  ATTR_ALIGN FP *ay  = nullptr;
  ATTR_ALIGN FP *az  = nullptr;
  ATTR_ALIGN FP *afx = nullptr;
  ATTR_ALIGN FP *afy = nullptr;
  ATTR_ALIGN FP *afz = nullptr;
  ATTR_ALIGN FP *as  = nullptr;
  ATTR_ALIGN FP *ae  = nullptr;
  ATTR_ALIGN FP *aq  = nullptr;
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
  void Confined(Molecule*,Atom*,const MolTypeList,const dvec3,const int,Property&,const int*,const int*);

  void Switching(Molecule*,Atom*,const MolTypeList,const dvec3,const int,Property&);
  template<int>
  void SwitchingTuning(Molecule*,Atom*,const MolTypeList,const dvec3,const int,Property&,const int*,const int *,const int,const int*,const int*);
  void Switching3site (Molecule*,Atom*,const MolTypeList,const dvec3,const int,Property&,const int*,const int *,const int,const int*,const int*);
  void Switching4site (Molecule*,Atom*,const MolTypeList,const dvec3,const int,Property&,const int*,const int *,const int,const int*,const int*);
  void Switching5site (Molecule*,Atom*,const MolTypeList,const dvec3,const int,Property&,const int*,const int *,const int,const int*,const int*);
};

#endif
