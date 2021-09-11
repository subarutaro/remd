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

class Molecules{
 private:

 public:
  const int nmol;
  const int natom;
  const double dt,dthalf;
  double T,P;

  Molecule   *mlcl;
  Atom       *atom;
  Thermostat *tst;
  Barostat   *bst;
  MolTypeList mlist;
  dvec3       L;
  Property    prop;
  Parameter   param;
  ForceCalculator *fclcltr;

  Profiler prof;

  int mode;

  int nthreads;

  Molecules
  (const Parameter _param)
    :nmol(_param.nmol),natom(_param.natom),
    dt(_param.dt*1e-15/unit_time), dthalf(_param.dt*1e-15/unit_time*0.5)
  {
    param = _param;

    T = param.temperature/unit_temp;
    P = param.pressure / unit_press;

    mlcl = new Molecule[nmol];
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
  };

  ~Molecules(){
    delete fclcltr;
    delete bst;
    delete tst;
    delete[] atom;
    delete[] mlcl;
  }

  void InitializeProperty(){
    CalcForcePot();
    prop.gkT = 6.0*(double)nmol*T;
    CalcProperties();
    prop.H0  = prop.tot;
    prop.ham = prop.tot - prop.H0;
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
  template <int> void D1();
  template <int> void D3();
  template <int> void D2();
  template <int> void D4();
  template <int> void D5();
  template <int> void D6();

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

  double GetPotential(){return prop.pot;};
  double GetKinetic(){return prop.kin;};
  double GetVolume(){
    if(param.confined == 1){
      return GetBottomArea() * L[2];
    }else if (param.confined == 2){
      const double wl = param.wall_length * 1e-10 / unit_length;
      return L[0] * L[1] * wl;
    }else{
      return L[0]*L[1]*L[2];
    }
  };
  double GetBottomArea(){
    if(param.confined == 1){
      const double sigma_c = 3.4; // angstrom
      const double wl = (param.wall_length * 1e-10 / unit_length - 0.5*sigma_c);
      return M_PI * wl * wl;
    }else{
      fprintf(stderr,"error: do not use GetBottomArea for not-1D confined system");
      exit(EXIT_FAILURE);
    }
  };
  double GetVirial(){return sum(prop.vir);}
  double GetHamiltonian(){return prop.tot;};
  Property GetProperty(){return prop;}

  void PrintAll(std::ostream&);
};
#endif

