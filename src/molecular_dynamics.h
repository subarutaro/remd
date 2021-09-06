#ifndef H_MOLECULAR_DYNAMICS
#define H_MOLECULAR_DYNAMICS

#include "io_manager.h"
#include "molecules.h"
#include "force_calculator.h"

struct MolecularDynamics{
  // parameters
  //double T;   // temperature
  //double P;   // pressure
  //dvec3  L;   // length of cell edge
  //double gkT;
  //double H0;

  //double dt;   // delta time
  //long   nstep;// total step
  Molecules  *mlcls;
  IOManager  *iomngr;

  Property pry;
  //double   Tave;
  //dvec3    Pave;
  //int      nave;

  MolecularDynamics();
  ~MolecularDynamics();

  MolecularDynamics(int argc, char **argv);

  void CalcProperties();
  void OutputProperty(std::ostream &s);
  //void Test();
};

#endif
