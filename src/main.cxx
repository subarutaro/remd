#include "molecular_dynamics.h"
#include "remd.h"
#if 1
int main(int argc,char **argv){
  //MolecularDynamics md(argc,argv);
  REMD md;
  md.ExecuteREMD();
  return 0;
};
#else //test
int main(){
  long nmol = 108;
  ForceCalculator *f;
  Molecules *mlcls;
  dvec3 L;
  L[0] = L[1] = L[2] = 10.;

  std::ostream *s;
  s = new std::iostream(std::cout.rdbuf());

  mlcls = new Molecules(nmol);
  mlcls->SetCubicFCC();

  mlcls->PrintAll(*s,L);

  //f->SPCE_Ewald(mlcls->mlcl,nmol,L);
  f->SPCE_Direct(mlcls->mlcl,nmol,L);

  mlcls->PrintAll(*s,L);
  delete s;

  return 0;
}
#endif
