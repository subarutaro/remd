#ifndef H_MOLECULES
#define H_MOLECULES

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <vector>
#include <string>

#include "remdinfo.h"
#include "vector.h"
#include "unit.h"
#include "rotation.h"
#include "force_calculator.h"
//#include "io_manager.h"
#include "parameter.h"
#include "product.h"
#include "variables.h"
#include "integrator.h"
#include "profiler.h"

#ifdef _OPENMP
#include <omp.h>
#endif

class Molecules{
 private:

 public:
  const int nmol;
  const int natom;
  const double dt,dthalf;
  double T,P;

  Molecule   *mlcl;
#ifdef ENABLE_AOS_TO_SOA_CONVERSION
  bool isAoS;
  double *mi;
  double *rx, *ry, *rz;
  double *vx, *vy, *vz;
  double *fx, *fy, *fz;
  double *qx, *qy, *qz, *qw;
  double *px, *py, *pz, *pw;
  double *nx, *ny, *nz, *nw;
  double *ix, *iy, *iz;
#endif

  Atom       *atom;
  Thermostat *tst;
  Barostat   *bst;
  MolTypeList mlist;
  dvec3       L;
  Property    prop;
  Average     ave;
  Parameter   param;
  ForceCalculator *fclcltr;

  Profiler prof;

  int mode;

  int nthreads;
  int* is;
  int* ie;
  // tmp for integrator
  double* sum_Ps_omp;
  dvec3*  sum_Pv_omp;

  //constexpr int nlane = nsimd*2;
  int* jstart;
  int* jend;

#if defined(__AVX512F__)
  const int nsimd = 64 / sizeof(FP);
#elif defined(__AVX2__)
  const int nsimd = 32 / sizeof(FP);
#else
  const int nsimd = 1;
#endif
  const int nlane = nsimd;

  Molecules
  (const Parameter _param)
    :nmol(_param.nmol),natom(_param.natom),
    dt(_param.dt*1e-15/unit_time), dthalf(_param.dt*1e-15/unit_time*0.5)
  {
    param = _param;

    T = param.temperature/unit_temp;
    P = param.pressure / unit_press;

    mlcl = new Molecule[nmol];
#if ENABLE_AOS_TO_SOA_CONVERSION
#ifdef __INTEL_COMPILER
    mi = (double*)_mm_malloc(sizeof(double)*nmol,64);
    rx = (double*)_mm_malloc(sizeof(double)*nmol,64);
    ry = (double*)_mm_malloc(sizeof(double)*nmol,64);
    rz = (double*)_mm_malloc(sizeof(double)*nmol,64);
    vx = (double*)_mm_malloc(sizeof(double)*nmol,64);
    vy = (double*)_mm_malloc(sizeof(double)*nmol,64);
    vz = (double*)_mm_malloc(sizeof(double)*nmol,64);
    fx = (double*)_mm_malloc(sizeof(double)*nmol,64);
    fy = (double*)_mm_malloc(sizeof(double)*nmol,64);
    fz = (double*)_mm_malloc(sizeof(double)*nmol,64);
    qx = (double*)_mm_malloc(sizeof(double)*nmol,64);
    qy = (double*)_mm_malloc(sizeof(double)*nmol,64);
    qz = (double*)_mm_malloc(sizeof(double)*nmol,64);
    qw = (double*)_mm_malloc(sizeof(double)*nmol,64);
    px = (double*)_mm_malloc(sizeof(double)*nmol,64);
    py = (double*)_mm_malloc(sizeof(double)*nmol,64);
    pz = (double*)_mm_malloc(sizeof(double)*nmol,64);
    pw = (double*)_mm_malloc(sizeof(double)*nmol,64);
    nx = (double*)_mm_malloc(sizeof(double)*nmol,64);
    ny = (double*)_mm_malloc(sizeof(double)*nmol,64);
    nz = (double*)_mm_malloc(sizeof(double)*nmol,64);
    nw = (double*)_mm_malloc(sizeof(double)*nmol,64);
    ix = (double*)_mm_malloc(sizeof(double)*nmol,64);
    iy = (double*)_mm_malloc(sizeof(double)*nmol,64);
    iz = (double*)_mm_malloc(sizeof(double)*nmol,64);
#else
    mi = new double[nmol];
    rx = new double[nmol];
    ry = new double[nmol];
    rz = new double[nmol];
    vx = new double[nmol];
    vy = new double[nmol];
    vz = new double[nmol];
    fx = new double[nmol];
    fy = new double[nmol];
    fz = new double[nmol];
    qx = new double[nmol];
    qy = new double[nmol];
    qz = new double[nmol];
    qw = new double[nmol];
    px = new double[nmol];
    py = new double[nmol];
    pz = new double[nmol];
    pw = new double[nmol];
    nx = new double[nmol];
    ny = new double[nmol];
    nz = new double[nmol];
    nw = new double[nmol];
    ix = new double[nmol];
    iy = new double[nmol];
    iz = new double[nmol];
#endif
    isAoS = true;
#endif
    atom = new Atom[natom];
    tst  = new Thermostat(param.tstat_mass);
    bst  = new Barostat(param.bstat_mass);

    nthreads = param.nthreads;

    //fclcltr = new ForceCalculator(param.rcut,param.kcut,param.alpha,param.nthreads);
    fclcltr = new ForceCalculator(param);
    //MakeMolTypeList();
    MakeMolTypeList(param.input_prefix+param.mtl_in);
    fclcltr->MakeParamList(mlist);

    mode = (param.confined<<CSHIFT) + (param.pconstant<<PSHIFT) + (param.tconstant<<TSHIFT);
    std::cout << "  mode: " << mode << std::endl;
    std::cout << " Tmode: " << ((mode>>TSHIFT)&MASK);
    std::cout << " Pmode: " << ((mode>>PSHIFT)&MASK);
    std::cout << " Cmode: " << ((mode>>CSHIFT)&MASK);
    std::cout << std::endl;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#else
    nthreads = 1;
#endif
    is = new int[nthreads];
    ie = new int[nthreads];
    is[0] = 0;
    ie[0] = (nmol%nthreads == 0) ? nmol/nthreads : nmol/nthreads + 1;

    for(int i=1;i<nthreads;i++){
      is[i] = ie[i-1];
      ie[i] = is[i] + ((nmol/nthreads)/nlane)*nlane;
    }
    ie[nthreads-1] = nmol;
    for(int i=0;i<nthreads;i++) printf("%d: (is,ie)= (%d,%d)\n",i,is[i],ie[i]);
    sum_Ps_omp  = new double[nlane*nthreads];
    sum_Pv_omp  = new  dvec3[nlane*nthreads];

    jstart = new int[(nmol+nlane-1)/nlane];;
    jend   = new int[(nmol+nlane-1)/nlane];;
  };

  ~Molecules(){
    delete fclcltr;
    delete bst;
    delete tst;
    delete[] atom;
#ifdef ENABLE_AOS_TO_SOA_CONVERSION
#ifdef __INTEL_COMPILER
    _mm_free(mi);
    _mm_free(rx);
    _mm_free(ry);
    _mm_free(rz);
    _mm_free(vx);
    _mm_free(vy);
    _mm_free(vz);
    _mm_free(fx);
    _mm_free(fy);
    _mm_free(fz);
    _mm_free(qx);
    _mm_free(qy);
    _mm_free(qz);
    _mm_free(qw);
    _mm_free(px);
    _mm_free(py);
    _mm_free(pz);
    _mm_free(pw);
    _mm_free(nx);
    _mm_free(ny);
    _mm_free(nz);
    _mm_free(nw);
    _mm_free(ix);
    _mm_free(iy);
    _mm_free(iz);
#else
    delete [] mi;
    delete [] rx;
    delete [] ry;
    delete [] rz;
    delete [] vx;
    delete [] vy;
    delete [] vz;
    delete [] fx;
    delete [] fy;
    delete [] fz;
    delete [] qx;
    delete [] qy;
    delete [] qz;
    delete [] qw;
    delete [] px;
    delete [] py;
    delete [] pz;
    delete [] pw;
    delete [] nx;
    delete [] ny;
    delete [] nz;
    delete [] nw;
    delete [] ix;
    delete [] iy;
    delete [] iz;
#endif
#endif
    delete[] mlcl;
  }

  void InitializeProperty(){
    #pragma omp parallel
    {
      CalcForcePot();
    }
    prop.gkT = 6.0*(double)nmol*T;
    CalcProperties();
    prop.H0  = prop.tot;
    prop.ham = prop.tot - prop.H0;
    std::cout << prop << std::endl;
  };

  void MakeMolTypeList();
  void MakeMolTypeList(std::string);
  void SetCubicFCC();
  void KillMomentum();
  void KillAngularMomentumInTube();
  void SetMassCenterToSystemCenter();
  void RandomVel();
  void VelocityScaling();
  void RandomAngularVel();
  void ZeroAngularVel();
  void AngVelocityScaling();
  void ResetAngularMomentum();

  void Sort();

  void ConvertToAtoms();
  void ConvertFromAtoms();

  void CalcForcePot();
#ifdef ENABLE_AOS_TO_SOA_CONVERSION
  void AoStoSoA();
  void SoAtoAoS();
#endif
  template <int> void init_D3D2D1D2D3_PsPv();
  template <int> void D3D2D1D2D3_PsPv();
  template <int> void D1();
  template <int> void D3();
  template <int> void D2();
  template <int> void D4();
  template <int> void D5();
  template <int> void D6();
  template<int MODE>
  void ExecuteStep();
  void ExecuteSteps();

  dvec3  TranslationalEnergy();
  double RotationalEnergy();
  double ThermostatEnergy(const double);
  double BarostatEnergy(const double);

  dvec3 Momentum();
  dvec3 RotationalMomentum();

  void CalcProperties();
  void UpdateHamiltonian(){
    prop.H0  = prop.tot;
    prop.ham = 0.0;
  }
  void PrintProperties(std::ostream &s){s << prop; s << " " << L[0]*L[1]*L[2] << std::endl;;};

  //void OutputCDV(const std::string,const long);
  //void OutputGRO(const std::string,const long,const dvec3);
  //void UpdateOutputs(){iomngr->UpdateOutputs(mlcl,tst,bst,prop,L,mlist);};
  //void WriteCDV(const long step){iomngr->WriteCDV(mlcl,mlist,L,step);};
  //void WriteCheckPoint(){iomngr->WriteCheckPoint(mlcl,tst,bst,L,prop);};

  //for remd
  double GetPressure(){return P;};
  double GetTemperature(){return T;};
  void   SetPressure(double _P){P = _P;};
  void   SetTemperature(double _T){T = _T; prop.gkT=6.0*(double)nmol*T;};
  void   ChangeTemperature(double _T){
    const double coef = sqrt(_T / T);
    for(int i=0;i<nmol;i++){
      mlcl[i].v *= coef / tst->s;
      mlcl[i].p *= coef / tst->s;
    }
    tst->Ps *= coef;
    if(param.pconstant == 1) bst->Pv[0] *= coef * tst->s;
    if(param.pconstant == 2) bst->Pv[2] *= coef * tst->s;
    if(param.pconstant == 3){
      bst->Pv[0] *= coef * tst->s;
      bst->Pv[1] *= coef * tst->s;
    }
    SetTemperature(_T);
    //reset s
    //without reset, s becomes smaller/larger for low/high temperature, and
    //s sometimes becomes smaller/larger and smaller/larger with replica exchanges
    tst->s = 1.0;
    //tst->Ps = 0.0;
    //bst->Pv = 0.0;
  };

  double GetPotential() const {return prop.pot;};
  double GetKinetic() const {return prop.kin;};
  double GetVolume() const {
    if(param.confined == 1){
      return GetBottomArea() * L[2];
    }else if (param.confined == 2){
      const double wl = param.wall_length * 1e-10 / unit_length;
      return L[0] * L[1] * wl;
    }else{
      return L[0]*L[1]*L[2];
    }
  };
  double GetBottomArea() const {
    if(param.confined == 1){
      const double sigma_c = 3.4; // angstrom
      const double wl = (param.wall_length - 0.5*sigma_c) * 1e-10 / unit_length;
      return M_PI * wl * wl;
    }else{
      fprintf(stderr,"error: do not use GetBottomArea for not-1D confined system");
      exit(EXIT_FAILURE);
    }
  };
  double GetVirial()      const {return sum(prop.vir);}
  double GetHamiltonian() const {return prop.tot;};
  Property GetProperty()  const {return prop;}
  Average  GetAverage()   const {
    Average ave;
    ave.sum[Average::LJ]   = prop.lj;
    ave.sum[Average::CLMB] = prop.clmb;
    ave.sum[Average::WALL] = prop.wall;
    ave.sum[Average::VOL]  = GetVolume();//prop.vol;
    ave.sum[Average::TRA]  = prop.tra;
    ave.sum[Average::ROT]  = prop.rot;
    ave.sum[Average::T]    = prop.T;
    ave.sum[Average::P]    = prop.prs[3];
    ave.sum[Average::Px]   = prop.prs[0];
    ave.sum[Average::Py]   = prop.prs[1];
    ave.sum[Average::Pz]   = prop.prs[2];
    ave.sum[Average::TOT]  = prop.tot;
    ave.sum[Average::DRF]  = prop.ham;
    ave.sum[Average::TSTs] = tst->s;
    ave.sum[Average::TSTv] = tst->Ps;
    ave.sum[Average::BSTx] = bst->Pv[0];
    ave.sum[Average::BSTy] = bst->Pv[1];
    ave.sum[Average::BSTz] = bst->Pv[2];

    return ave;
  }

  void PrintAll(std::ostream&);
};
#endif

