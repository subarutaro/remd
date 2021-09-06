#include "molecular_dynamics.h"
#include "timer.h"
#include "parameter.h"

MolecularDynamics::MolecularDynamics(int argc,char **argv){

  Parameter param;
  param.read(std::cin);
  param.print(std::cout);
  mlcls = new Molecules(param);
  iomngr = new IOManager(param);

  //mlcls->SetCubicFCC();
  //mlcls->RandomVel();
  //mlcls->RandomAngularVel();
  //mlcls->VelocityScaling();
  //mlcls->AngVelocityScaling();
  //mlcls->CalcForcePot();

  std::ostream *o_eng;
  o_eng = new std::ofstream(param.ene_out.c_str(),std::ios_base::out);
  for(int s=1;s<=param.ninterval;s++){
    for(int i=0;i<param.interval;i++){
      mlcls->ExecuteSteps();
    }
    if(s%param.energy_interval==0){
      mlcls->CalcProperties();
      mlcls->PrintProperties(*o_eng);
    }
    if(s%param.snap_interval==0){
      iomngr->UpdateOutputs(mlcls);
      //mlcls->WriteCDV(s/param.snap_interval);
    }
    if(s%param.mom_cancel_interval==0){
      mlcls->KillMomentum();
    }
  }
  iomngr->UpdateOutputs(mlcls);
#if 0
  CalcProperties();
  H0  = pry.ham;

  /*
  const double rcut  = 9.0;
  const double rtol  = 1.0e-5;
  const double fs    = 6.0; // fourierspacing
  const double kcut  = L[0] / fs;
  const double alpha = sqrt(-log(rtol)/(rcut*rcut));
  fclcltr = new ForceCalculator(rcut,kcut,alpha);
  //std::cout << "alpha and kcut = " << alpha << " " << kcut << std::endl;
  //*/

  nstep = 0;
  int nstepout = 0;
  std::ostream *o_eng, *o_dbg;
  o_eng = new std::ofstream(iomngr.iprm->oene,std::ios_base::out); //file output
  //o_eng = new std::iostream(std::cout.rdbuf()); //stdio
  o_dbg = new std::iostream(std::cout.rdbuf()); //stdio
  std::string cdv_prefix = "cdv/";

  int nstep_precalc = 0;
  for(int s=0;s<nstep_precalc;s++){
    mlcls->VelocityScaling(T);
    mlcls->AngVelocityScaling(T);
    CalcProperties();
    ExecuteSteps();
  }
  pry.time = 0.0;
  CalcProperties();
  H0  = pry.ham;
  for(int s=0;s<nstep;s++){
    CalcProperties();
    if(s%nstepout == 0){
      OutputProperty(*o_eng);
      mlcls->OutputCDV(s/nstepout);
      iomngr.UpdateOutputs(mlcls->mlcl,mlcls->mlist,mlcls->L);
      //mlcls->PrintAll(*o_dbg,mlcls->L);
    }
    ExecuteSteps();
  }
  iomngr.UpdateOutputs(mlcls->mlcl,mlcls->mlist,mlcls->L);
  //iomngr.PrintOptions();
#endif
};

MolecularDynamics::~MolecularDynamics(){
  //delete mlcls;
  std::cout << "MD simulation successfully finished!" << std::endl;
};


void MolecularDynamics::OutputProperty(std::ostream &s){
  if(pry.time == 0.){
    s << "#     1.time";
#ifdef DEBUG
    s << "        2.vdw";
    s << "        3.clm";
#endif
    s << "        4.pot";
    s << "        5.tra";
    s << "        6.rot";
    s << "        7.kin";
    s << "        8.tsm";
    s << "        9.tsp";
    s << "       10.tst";
    s << "       11.bst";
    s << "       12.tot";
    s << "       13.ham";
    s << "     13.tmo_x";
    s << "     14.tmo_y";
    s << "     15.tmo_z";
    /*
    s << "   16.rmo_x";
    s << "   17.rmo_y";
    s << "   18.rmo_z";
    s << "   19.vir_x";
    s << "   20.vir_y";
    s << "   21.vir_z";
    s << "   22.prs_x";
    s << "   23.prs_y";
    s << "   24.prs_z";
    //*/
    s << std::endl;
  }
  s << pry << std::endl;
}
/*
void MolecularDynamics::CalcForce(){
  double* forcepot;
  if((forcepot = (double*)malloc((unsigned int)(sizeof(double)*4*mlcls->natom)))==NULL){
    fprintf(stderr,"error: malloc forcepot failed\n");
    exit(EXIT_FAILURE);
  }
  for (int i = 0; i < 4*mlcls->natom; i++){
    forcepot[i] = 0.e0;
  }
  for(int j=0;j<mlcls->nmol;j++){
    mlcls->mlcl[j].f = 0.0;
    mlcls->mlcl[j].n = 0.0;
  }
  Time lj("LJ"),ew("EW");
  lj.beg();
  fclcltr->LJ(mlcls->mlcl,mlcls->mlist,pry,mlcls->nmol,mlcls->L,forcepot);
  lj.end();
  ew.beg();
  fclcltr->Ewald(mlcls->nmol,mlcls->natom,mlcls->mlcl,mlcls->mlist,pry,mlcls->L,forcepot);
  ew.end();
  lj.print();ew.print();

  //get force, torque  and pot from forcepot
  double tmp = 0.;
  //#pragma omp parallel reduction (+: tmp)
  for (int j = 0; j < mlcls->nmol; j++){
    Molecule m = mlcls->mlcl[j];
    MolType mt = mlcls->mlist[m.type];
    for(unsigned int d = 0;d<mt.a.size();d++){
      dvec3 f;
      f[0] = forcepot[4*(m.id + d)+0];
      f[1] = forcepot[4*(m.id + d)+1];
      f[2] = forcepot[4*(m.id + d)+2];
      //std::cout << (m.id + d) << " " << f << std::endl;
      m.f += f;

      const dvec3 fbody = space_to_body(m.q,f);
      m.n[0] += scalar_prod(mt.a[d].r,fbody);
      m.n[1] += mt.a[d].r[1]*fbody[2] - mt.a[d].r[2]*fbody[1];
      m.n[2] += mt.a[d].r[2]*fbody[0] - mt.a[d].r[0]*fbody[2];
      m.n[3] += mt.a[d].r[0]*fbody[1] - mt.a[d].r[1]*fbody[0];

      tmp += forcepot[4*(m.id + d)+3];
    }
    mlcls->mlcl[j] = m;
  }
  pry.pot = tmp;
#ifdef DEBUG
  pry.clm = pry.drct + pry.wave + pry.self + pry.intra;
  std::cout << "potential of van der waals:" << std::endl;
  std::cout << "vdw=   " << pry.vdw << std::endl;
  std::cout << "potential of ewald:"<< std::endl;
  std::cout << "drct= " << pry.drct << std::endl;
  std::cout << "wave= " << pry.wave << std::endl;
  std::cout << "self= " << pry.self << std::endl;
  std::cout << "intra=" << pry.intra<< std::endl;
  std::cout << "clm=  " << pry.clm  << std::endl;
  //pry.pot = pry.vdw + pry.clm;//uptdate potential energy for D5
#endif
  free(forcepot);
}
//*/
/*
void MolecularDynamics::ExecuteSteps(){
  Time t("total"),f("force"),i("integrator");
  t.beg();i.beg();
  mlcls->D6(0.5*dt);
  mlcls->D5(0.5*dt,P,pry.pot,pry.vir);
  mlcls->D4(0.5*dt);
  mlcls->D3(0.5*dt);
  mlcls->D2(0.5*dt);
  mlcls->D1(dt,gkT,H0);
  mlcls->D2(0.5*dt);
  mlcls->D3(0.5*dt);
  mlcls->D4(0.5*dt);
  i.end();

  f.beg();
  CalcForce();
  f.end();

  i.beg();
  mlcls->D5(0.5*dt,P,pry.pot,pry.vir);
  mlcls->D6(0.5*dt);
  pry.time += unit_time * dt;
  i.end();
  t.end();
  //i.print();f.print();
  t.print();
}
//*/

/*
void MolecularDynamics::Test(){
  SPCE spce;
  I[0] = spce.m[1]*spce.r[1][1]*spce.r[1][1] + spce.m[2]*spce.r[2][1]*spce.r[2][1];
  I[1] = spce.m[0]*spce.r[0][0]*spce.r[0][0] + spce.m[1]*spce.r[1][0]*spce.r[1][0] + spce.m[2]*spce.r[2][0]*spce.r[2][0];;
  I[2] = spce.m[0]*scalar_prod(spce.r[0],spce.r[0]) + spce.m[1]*scalar_prod(spce.r[1],spce.r[1]) + spce.m[2]*scalar_prod(spce.r[2],spce.r[2]);
  std::cout << "I = " << I[0] << " " << I[1] << " " << I[2] << std::endl;

  nmol = 32;
  const double mass_water = spce.m[0] + spce.m[1] + spce.m[2];
  const double density = 0.997e3 / (unit_mass / (unit_length*unit_length*unit_length));
  L[0] = L[1] = L[2] = pow(mass_water*(double)nmol/density,1./3.);
  //L[0] = L[1] = L[2] = 10.;
  T = 200 / unit_temp;
  //P = 101.3e3 * unit_length*unit_length*unit_length / unit_energy;
  P = 0.0;
  //gkT = (6.0*(double)nmol+2.0)*T;
  gkT = (6.0*(double)nmol)*T;

  dt = 2.0e-15 / unit_time;
  rc = 0.5*L[0];
  //rc = 10.;

  const double rtol  = 1.0e-5;
  const double kcut  = -log(rtol) * L[0] / (M_PI*rc);
  const double alpha = sqrt(M_PI*kcut/(L[0]*rc));
  std::cout << "kcut = " << kcut << ", alpha = " << alpha << std::endl;

  mlcls = new Molecules(nmol);
#if 0
  //set intitial coordinate
  mlcls->SetMass(mass_water);
  mlcls->SetCubicFCC();
  //set intitial velocity
  mlcls->RandomVel(L);
#else
  
#endif
  mlcls->VelocityScaling(T,L,I);
  //set intitial angular velocity
  mlcls->RandomAngularVel();
  mlcls->AngVelocityScaling(T,L,I);
  //calclulate initial force and potential
  fclcltr = new ForceCalculator(rc,kcut,alpha);
  fclcltr->SPCE_LJ(mlcls->mlcl,pry,nmol,L);
  fclcltr->SPCE_Ewald(nmol,natom,mlcls->mlcl,mlcls->mlist,pry,L);

  //calculate initial properties
#ifdef NPT
  CalcProperties();
  P = norm(pry.prs);
#endif
  CalcProperties();
  H0 = pry.tot;

  std::cout << "L= " << L[0] << " " << L[1] << " " << L[2] << std::endl;
  std::cout << "dt= " << dt << ", rc= " << rc << std::endl;
  std::cout << "eps= " << spce.eps << " sgm = " << spce.sgm << std::endl;
  std::cout << "T= " << T << " P = " << P << std::endl;
  std::cout << "gkT= " << gkT << std::endl;

  //set output source
  std::ostream *o_eng, *o_dbg;
  o_eng = new std::ofstream("eng.txt",std::ios_base::out); //file output
  //o_eng = new std::iostream(std::cout.rdbuf()); //stdio
  o_dbg = new std::iostream(std::cout.rdbuf()); //stdio
  std::string cdv_prefix = "cdv/test";

  nstep = 1;
  const long out_interval = 1;
  bool scaling_flg = false;

  std::cout << "roop start" << std::endl;
  for(long s=0;s<nstep;s++){
    if(s<nstep/10){
      mlcls->tst->s  = 1.0;
      mlcls->tst->Ps = 0.0;

      mlcls->VelocityScaling(T,L,I);
      mlcls->AngVelocityScaling(T,L,I);
#ifdef NPT
      CalcProperties();
      P = norm(pry.prs);
#endif
      CalcProperties();
      H0 = pry.tot;
    }

    if(s%out_interval==0){
      std::cout << "output " << s << std::endl;
      mlcls->PrintAll(*o_dbg,L);
      mlcls->OutputCDV(cdv_prefix,s/out_interval,L);
      mlcls->OutputGRO(cdv_prefix,s/out_interval,L);
      CalcProperties();
      OutputProperty(*o_eng);
    }
    ExecuteSteps();
  }

  std::cout << "Tave= " << Tave/(double)nave << "Pave= " << Pave/(double)nave << std::endl;

  delete o_eng;
  delete o_dbg;

  delete mlcls;
}
//*/
