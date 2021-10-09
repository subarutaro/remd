
#include<cstdio>
#include<cstdlib>
#include<cmath>

#include <iostream>

#include "mdunit.h"
#include "remd.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

//===================== class REMD ======================
//===================== (de)constructor of REMD ======================

void Abort(const int ret = 0){
#ifdef ENABLE_MPI
  MPI_Finalize();
  exit(ret);
#endif
}

REMD::REMD(std::string input){
  std::ifstream ifs(input);
  param.read(ifs);
  param.print(std::cout);

  remdinfo = new REMDInfo(param);

  const int nreplica = remdinfo->GetNumReplica();
  const int offset   = remdinfo->GetReplicaOffset(); // for MPI

  assert(!(param.confined == 1 && (param.pconstant == 1 || param.pconstant == 3)));
  assert(!(param.confined == 2 && (param.pconstant == 1 || param.pconstant == 2)));

  if(param.ninterval > 0){
    //struct pointer of pointer *md
    md = new Molecules*[nreplica];
    iomngr = new IOManager*[nreplica];
    for(int i=0;i<nreplica;i++){
      std::stringstream strs;
      strs << std::setw(4) << std::setfill('0') << i+offset;
      Parameter p = param;
      if(p.gro_in  != "null") p.gro_in  += strs.str();
      if(p.trr_in  != "null") p.trr_in  += strs.str();
      if(p.chk_in  != "null") p.chk_in  += strs.str();
      if(p.tkw_in  != "null") p.tkw_in  += strs.str();
      if(p.xyz_in  != "null") p.xyz_in  += strs.str();
      if(p.cdv_in  != "null") p.cdv_in  += strs.str();

      if(p.gro_out != "null") p.gro_out += strs.str();
      if(p.trr_out != "null") p.trr_out += strs.str();
      //if(p.chk_out != "null") p.chk_out += strs.str();
      if(p.ene_out != "null") p.ene_out += strs.str();
      if(p.cdv_out != "null") p.cdv_out += strs.str();

      iomngr[i] = new IOManager(p);
      md[i]     = new Molecules(p);
      if(p.gro_in != "null" || p.chk_in != "null"){
	iomngr[i]->ReadInputs(md[i]);
      }
      if(p.adjust_center==1)
	md[i]->SetMassCenterToSystemCenter();

      if(p.tconstant > 0){
	md[i]->SetTemperature(remdinfo->GetTemperature(i));
      }
      if(p.pconstant > 0){
	md[i]->SetPressure(remdinfo->GetPressure(i));
      }

      if(p.gen_vel == 1){
	md[i]->RandomVel();
	md[i]->VelocityScaling();
	//md[i]->ZeroAngularVel();
	md[i]->RandomAngularVel();
	md[i]->AngVelocityScaling();
      }

      if(p.init_scaling == 1){
	md[i]->VelocityScaling();
	md[i]->AngVelocityScaling();
      }

      if(p.bstat_mass != 0.0){
	md[i]->bst->W = p.bstat_mass;
	md[i]->bst->flush();
      }
      if(p.tstat_mass != 0.0){
	md[i]->tst->Q = p.tstat_mass;
	md[i]->tst->flush();
      }
      md[i]->InitializeProperty();

      //md[i]->Sort();
    }

    /*
      for(int i=0;i<nreplica;i++){
      md[i]->VelocityScaling();
      md[i]->AngVelocityScaling();
      md[i]->InitializeProperty();
      md[i]->PrintProperties(std::cout);
      }
    //*/
#ifdef ENABLE_GPU
    ngpus = param.ngpus;
    if(ngpus>0){
      fprintf(stdout,"GPU is used. # of gpu is %d, offest is %d\n",param.ngpus,param.gpu_offset);
      mdgpu = new REMDGPU(md,nreplica);
    }
#else
    ngpus = 0;
#endif
    snap_interval   = param.snap_interval;
    chk_interval    = param.chk_interval;
    energy_interval = param.energy_interval;
    step = 0;

    const int rank = remdinfo->GetRank();
    if(rank == 0){
      const int nreplica_global = remdinfo->GetNumReplicaGlobal();
      pot = new D[nreplica_global];
      vol = new D[nreplica_global];
      ave = new Average[nreplica_global];
    }else{
      pot = new D[nreplica];
      vol = new D[nreplica];
      ave = new Average[nreplica];
    }
    for(int i=0;i<nreplica;++i){
      pot[i] = md[i]->GetPotential();
    }
    for(int i=0;i<nreplica;++i){
      vol[i] = md[i]->GetVolume();
    }
    for(int i=0;i<nreplica;++i){
      ave[i].flush();
    }
    GatherResult();
    if(remdinfo->GetRank() == 0){
      const int nreplica_global = remdinfo->GetNumReplicaGlobal();
      std::stringstream strs;
      strs << remdinfo->GetOutputDir() << "/energy.dat";
      os_ene = new std::ofstream(strs.str().c_str(),std::ios_base::out);
      if(os_ene->fail()){
	std::cout << "error: open " << strs.str() << "failed" << std::endl;
	exit(EXIT_FAILURE);
      }
      binary = new std::fstream[nreplica_global];
      for(int rep=0;rep<nreplica_global;rep++){
	std::stringstream strs;
	strs << remdinfo->GetOutputDir() << "/P" << remdinfo->GetPressure(rep) << 'T' << remdinfo->GetTemperature(rep) << ".log";
	const std::string filename = strs.str();
	std::cout << "opening " << filename << std::endl;
	if(param.restart == 1){
	  binary[rep].open(filename,std::ios::in  | std::ios::binary);
	  if(binary[rep].fail()) {
	    std::cerr << "error: open " << filename << " failed" << std::endl;
	    exit(EXIT_FAILURE);
	  }
	  double tmp[3];
	  binary[rep].read((char*)tmp,sizeof(double)*3);
	  if(tmp[0] != remdinfo->GetTemperature(rep) || tmp[1] != remdinfo->GetPressure(rep)){
	    std::cerr << "error: open " << filename << " failed" << std::endl;
	    exit(EXIT_FAILURE);
	  }
	  binary[rep].close();

	  binary[rep].open(filename,std::ios::out | std::ios::binary | std::ios::app);
	  if(binary[rep].fail()){
	    std::cout << "error: open " << filename << "failed" << std::endl;
	    exit(EXIT_FAILURE);
	  }
	}else{
	  binary[rep].open(filename,std::ios::out | std::ios::binary);
	  if(binary[rep].fail()){
	    std::cout << "error: open " << filename << "failed" << std::endl;
	    exit(EXIT_FAILURE);
	  }
	  const double conditions[3]
	    = {remdinfo->GetTemperature(rep),remdinfo->GetPressure(rep),1.0};
	  binary[rep].write((char*)conditions,3*sizeof(double));
	  const int nelem = Average::NUM_ELEM;
	  binary[rep].write((char*)&nelem,sizeof(int));
	}
      }
    }
  } // if (ninterval > 0)
};

REMD::~REMD(){
  const int nreplica = remdinfo->GetNumReplica();
  for(int rep=0;rep<nreplica;rep++){
    Profiler p = md[rep]->prof + md[rep]->fclcltr->prof;
    p.print();
  }

  if(os_ene != nullptr) delete os_ene;
  if(wham != nullptr){
    std::cout << "delete wham" << std::endl;
    delete wham;
  }
  if(binary != nullptr){
    for(int rep=0;rep<param.nreplica;rep++) if(binary[rep].is_open()) binary[rep].close();
    delete[] binary;
  }
  if(pot != nullptr){
    std::cout << "delete pot" << std::endl;
    delete[] pot;
  }
  if(vol != nullptr){
    std::cout << "delete vol" << std::endl;
    delete[] vol;
  }
  if(vir != nullptr){
    std::cout << "delete vir" << std::endl;
    delete[] vir;
  }
  if(md != nullptr){
    const int nreplica = remdinfo->GetNumReplica();
    for(int i=0;i<nreplica;++i){
      std::cout << "delete md " << i << std::endl;
      delete md[i];
    }
    delete[] md;
  }
#ifdef ENABLE_GPU
  if(mdgpu != nullptr) delete mdgpu;
#endif
};

void REMD::ExecuteMDs(){
  const int nreplica = remdinfo->GetNumReplica();
  const int rank     = remdinfo->GetRank();
  if(ngpus>0){
#ifdef ENABLE_GPU
    mdgpu->ExecuteSteps(remdinfo->temperature,remdinfo->pressure);
    for(int rep=0;rep<nreplica;++rep){
      ave[rep] = mdgpu->GetAverages(rep);
      pot[rep] = ave[rep].sum[Average::LJ] + ave[rep].sum[Average::CLMB] + ave[rep].sum[Average::WALL];
      vol[rep] = ave[rep].sum[Average::VOL];
    }
#else
    fprintf(stderr,"error: gpu is not enabled\n");
    exit(EXIT_FAILURE);
#endif
  }else{
    //#pragma omp parallel for num_threads(param.nthreads)
    for(int rep=0;rep<nreplica;++rep){
      if(param.tconstant == 2){
	md[rep]->VelocityScaling();
	md[rep]->AngVelocityScaling();
      }
      md[rep]->ExecuteSteps();//interval is included in md class
      pot[rep] = md[rep]->GetPotential();
      vol[rep] = md[rep]->GetVolume();
      ave[rep] = md[rep]->GetAverage();
    }
  }
}

//====================== functions of REMD ====================
void REMD::ExecuteREMD(){
  std::cout << "start REMD" << std::endl;
  const unsigned long step_max = remdinfo->GetStepMax();
  const unsigned long interval = remdinfo->GetInterval();

  const int nreplica = remdinfo->GetNumReplica();
  //const int nprocs   = remdinfo->GetNumProc();
  const int rank = remdinfo->GetRank();
  //const int num_bkup = remdinfo->GetNumBackup();
  //const int num_ene  = remdinfo->GetNumPointEnergy();

  if(param.adjustor == 1){
    std::cout << "adjusting conditon" << std::endl;
    do{
      for(int s=0;s<param.adjust_interval;s++){
	if(s%100==0) std::cout << s << "th step" << std::endl;
	ExecuteMDs();
	//IncrementHistogram();
	ReplicaExchange();
	step++;
      }
    }while(remdinfo->ConditionAdjustor());
    step = 0;
  }
  
  for(unsigned long s=0;s<step_max/interval;++s){
    ExecuteMDs();
#if 0
    if(s%1000) std::cerr << "=== " << rank << " " << s << " intervals done ===" << std::endl;
#else
    for(int rep=0;rep<nreplica;rep++){
      std::cout << " rep " << rep;
      std::cout << " " << s;
      std::cout << " " << pot[rep];
      std::cout << " " << vol[rep];
      if(ngpus>0){
#ifdef ENABLE_GPU
 	Average a = mdgpu->GetAverages(rep);
	std::cout << " " << a.sum[Average::LJ];
	std::cout << " " << a.sum[Average::CLMB];
	std::cout << " " << a.sum[Average::WALL];
	std::cout << " " << a.sum[Average::T];
	if(param.confined == 1){
	  std::cout << " " << a.sum[Average::Pz];
	}else if(param.confined == 2){
	  std::cout << " " << a.sum[Average::Px];
	  std::cout << " " << a.sum[Average::Py];
	}
	std::cout << " " << a.sum[Average::TRA];
	std::cout << " " << a.sum[Average::ROT];
	std::cout << " " << a.sum[Average::TSTs];
	std::cout << " " << a.sum[Average::TSTv];
	if(param.confined == 1){
	  std::cout << " " << a.sum[Average::BSTz];
	}else if(param.confined == 2){
	  std::cout << " " << a.sum[Average::BSTx];
	  std::cout << " " << a.sum[Average::BSTy];
	}
#if 0
	if(param.tconstant == 1){
	  std::cout << " " << a.sum[Average::TST]
	}
	if(param.pconstant > 0){
	  std::cout << " " << a.sum[Average::BST]
	}
#endif
	std::cout << " " << a.sum[Average::DRF];
	std::cout << " " << a.sum[Average::TOT];
#else
        fprintf(stderr,"error: gpu is not enabled!");
	exit(EXIT_FAILURE);
#endif
      }else{
        md[rep]->CalcProperties();
        const Property prop = md[rep]->GetProperty();
        std::cout << " " << prop.lj;
        std::cout << " " << prop.clmb;
        std::cout << " " << prop.wall;
        std::cout << " " << prop.T;
        if(param.confined == 1){
	  std::cout << " " << prop.prs[2];
        }else if(param.confined == 2){
	  std::cout << " " << prop.prs[0];
	  std::cout << " " << prop.prs[1];
        }
        std::cout << " " << prop.tra;
        std::cout << " " << prop.rot;
        //std::cout << " " << prop.kin;
        //std::cout << " " << prop.tsm;
        //std::cout << " " << prop.tsp;
        //std::cout << " " << prop.tst;
        std::cout << " " << md[rep]->tst->s;
        std::cout << " " << md[rep]->tst->Ps;
        if(param.confined == 1){
	  std::cout << " " << md[rep]->bst->Pv[2];
        }else if(param.confined == 2){
	  std::cout << " " << md[rep]->bst->Pv[0];
	  std::cout << " " << md[rep]->bst->Pv[1];
        }
#if 0
	if(param.tconstant == 1){
	  std::cout << " " << prop.tst;
	}
	if(param.pconstant > 0){
	  std::cout << " " << prop.bst;
	}
#endif
        std::cout << " " << prop.ham;
        std::cout << " " << prop.tot;
#if 0
	std::cout << " " << prop.Tave / (double)prop.nave;
	if(param.confined == 1){
	  std::cout << " " << prop.Pave[2] / (double)prop.nave;
	}else if(param.confined == 2){
	  std::cout << " " << prop.Pave[0] / (double)prop.nave;
	  std::cout << " " << prop.Pave[1] / (double)prop.nave;
	}
#endif
      } // gpu or not
      std::cout << std::endl;
    }
#endif
    GatherResult(); //for mpi
    if(rank == 0){
      //OutputDebugInfo();
      //IncrementHistogram();
      OutputBinary();
    }
    if(param.rem_type==0){
      ReplicaExchange();
    }else if(param.rem_type==1){
      DesignedWalkReplicaExchange();
    }else if(param.rem_type < 0){
      // no replica exchange
    }else{
      std::cerr << "error: unsupported replica exhcange type" << std::endl;
      exit(EXIT_FAILURE);
    }
    //OutputBinary();
    if(step%snap_interval==0){
      if(ngpus>0){
#ifdef ENABLE_GPU
	mdgpu->OutputCDV((param.output_prefix+param.cdv_out).c_str(),step/snap_interval,remdinfo->index);
#endif
      }else{
	for(int rep=0;rep<nreplica;++rep){
	  iomngr[rep]->WriteCDV(md[rep],step/snap_interval);
	}
      }
    }
    if(step%chk_interval==0){
#ifdef ENABLE_GPU
      //mdgpu->OutputCheckPoint((param.output_prefix+param.chk_out).c_str(),remdinfo->index);
      if(ngpus>0) mdgpu->ConvertVariablesToHost(md);
#endif
      for(int rep=0;rep<nreplica;++rep){
	const int ind = remdinfo->GetIndex(rep);
	iomngr[rep]->UpdateOutputs(md[rep],ind);
      }
      //remdinfo->OutputHistogram();
    }
#if 0
    if(step%energy_interval==0){
      remdinfo->OutputAcceptRate();
      remdinfo->OutputTunnelCount();
      OutputEnergy();
    }
#endif
    step++;
  } // loop of s
  if(step_max > 0){
#ifdef ENABLE_GPU
    if(ngpus>0) mdgpu->ConvertVariablesToHost(md);
#endif
    for(int rep=0;rep<nreplica;++rep){
      const int ind = remdinfo->GetIndex(rep);
      iomngr[rep]->UpdateOutputs(md[rep],ind);

      iomngr[rep]->WriteCDV(md[rep],step/snap_interval);
    }
    if(rank  == 0){
      //output result
      //remdinfo->OutputHistogram();
      remdinfo->OutputAcceptRate();
      remdinfo->OutputTunnelCount();
      OutputEnergy();
      OutputIndex();
    }
  } // if step_max > 0
  //ExecuteWHAM();
}

static const D kb = 1.0;
static D random_number(){return ((D)rand() / (D)RAND_MAX);};

D REMD::CalcExchangeProb(int type,int ind){
  const Pair pair = remdinfo->GetPair(type,ind);
  const int ind1    = remdinfo->GetIndex(pair.a);
  const int ind2    = remdinfo->GetIndex(pair.b);

  const D beta1  = 1./remdinfo->GetTemperature(ind1);
  const D beta2  = 1./remdinfo->GetTemperature(ind2);
  const D press1 = remdinfo->GetPressure(ind1);
  const D press2 = remdinfo->GetPressure(ind2);

  const int ensemble = remdinfo->GetEnsemble();
  D delta;
  delta = (beta2 - beta1)*(pot[ind1]-pot[ind2]);
  if(((ensemble>>PSHIFT) & MASK) > 0){
    delta += (beta2*press2 - beta1*press1)*(vol[ind1]-vol[ind2]);
  }
  return exp(-delta);
};

void REMD::ReplicaExchange(){
  const int nreplica    = remdinfo->GetNumReplica();

  if(false) if(remdinfo->GetRank() == 0){
    printf("(debug:result)");
    for(int i=0;i<remdinfo->GetNumReplicaGlobal();i++){
      printf(" (%d,%lf,%lf)",remdinfo->GetIndex(i),pot[i],vol[i]);
    }
    printf("\n");
  }

  //exchange temperature and pressure
  if(remdinfo->GetRank() == 0){
    remdinfo->SetAllIsExchangedGlobal(false);
    const int ec_type     = step%4;//set exchange pair list number(0 to 3 for 2D rem,0 or 1 for 1D)
    const int list_length = remdinfo->GetListLength(ec_type);//set exhcange list length of ec_type
    for(int i=0;i<list_length;++i){
      const Pair pair = remdinfo->GetPair(ec_type,i);
      if(CalcExchangeProb(ec_type,i)>random_number()){
	remdinfo->ExchangeAccepted(ec_type,i);
      }else{
	remdinfo->ExchangeRejected(ec_type,i);
      }
    }
    remdinfo->CheckTunneling();
  }
  remdinfo->BroadcastConditionAndIndex();
#ifdef ENABLE_MPI
  MPI_Barrier(MPI_COMM_WORLD);

  if(false) if(remdinfo->GetRank() == 0){
    printf("(debug:condition)");
    for(int i=0;i<remdinfo->GetNumReplicaGlobal();i++){
      printf(" (%d,%lf,%lf,%d)",remdinfo->GetIndex(i),remdinfo->GetTemperature(i),remdinfo->GetPressure(i),remdinfo->GetIsExchanged(i));
    }
    printf("\n");
  }
#endif
  //set temperature of all replica from remdinfo
  //remdinfo has conditions of i th replica in order
  if(remdinfo->GetNumGPU()<=0){
    for(int i=0;i<nreplica;++i){
      if(remdinfo->GetIsExchanged(i)){
	md[i]->ChangeTemperature(remdinfo->GetTemperature(i));
	md[i]->KillMomentum();
	md[i]->SetPressure(remdinfo->GetPressure(i));
	// update hamiltonian for integration
	md[i]->CalcProperties();
	md[i]->UpdateHamiltonian();
      }
    }
  }
};

void REMD::DesignedWalkReplicaExchange(){
#ifdef ENABLE_MPI
  fprintf(stderr,"DWREM is not supported with MPI\n");
  Abort();
#endif
  const int nreplica    = remdinfo->GetNumReplica();
  const int ec_type     = remdinfo->GetListType();
  const int list_length = remdinfo->GetListLength(ec_type);//set exhcange list length of ec_type

  //exchange temperature and pressure
  bool isAllExchanged = true;
  for(int i=0;i<list_length;++i){
    const Pair pair  = remdinfo->GetPair(ec_type,i);
    const int  ind1  = remdinfo->GetIndex(pair.a);
    const int  ind2  = remdinfo->GetIndex(pair.b);
    const bool isex1 = remdinfo->GetIsExchanged(ind1);
    const bool isex2 = remdinfo->GetIsExchanged(ind2);
    if(!isex1 && !isex2){
      if(CalcExchangeProb(ec_type,i)>random_number()){
	remdinfo->ExchangeAccepted(ec_type,i);
	if(remdinfo->GetNumGPU()<=0){
	  md[ind1]->ChangeTemperature(remdinfo->GetTemperature(pair.a));
	  md[ind2]->ChangeTemperature(remdinfo->GetTemperature(pair.b));
	  md[ind1]->SetPressure(remdinfo->GetPressure(pair.a));
	  md[ind2]->SetPressure(remdinfo->GetPressure(pair.b));
	}
      }else{
	remdinfo->ExchangeRejected(ec_type,i);
	isAllExchanged = false;
      }
    }
  }
  if(isAllExchanged){
    printf("index:");for(int i=0;i<nreplica;i++) printf(" %2d",remdinfo->GetIndex(i));printf("\n");
    remdinfo->MoveToNextList();
    remdinfo->SetAllIsExchanged(false);
  }
  remdinfo->CheckTunneling();
};

template <class T>
static inline T max(T a,T b){
  if(a>b) return a;
  return b;
};

template <class T>
static inline T min(T a,T b){
  if(a<b) return a;
  return b;
};

template <class T>
static inline T min(T a,T b,T c){
  return min( min(a,b),c );
};

template <class T>
static inline T min(T a,T b,T c,T d){
  return min( min( min(a,b),c ),d );
};

static inline double sum_log(double a,double b){
  if(a==0. && b==0.) return a+b;
  if(a>b) return a + log(1. + exp(b-a));
  return b + log(1. + exp(a-b));
};

static const
int factorial(int N){
  if(N<1) return 1;
  return N * factorial(N-1);
};

static inline void CalcWeight(double *t,double *u,double *w,int **ind,int n){
  double b[n];
  for(int i=0;i<n;i++){
    b[i] = 1e0 / t[i];
  }
  for(int i=0;i<factorial(n);i++){
    w[i] = 0e0;
    for(int j=0;j<n;j++){
      w[i] -= u[ind[i][j]]*b[j];
    }
  }

  /* for N=4
  w[0] = -(u[0]*b[0] + u[1]*b[1] + u[2]*b[2] + u[3]*b[3]);
  w[1] = -(u[0]*b[0] + u[1]*b[1] + u[3]*b[2] + u[2]*b[3]);
  w[2] = -(u[0]*b[0] + u[2]*b[1] + u[1]*b[2] + u[3]*b[3]);
  w[3] = -(u[0]*b[0] + u[2]*b[1] + u[3]*b[2] + u[1]*b[3]);
  w[4] = -(u[0]*b[0] + u[3]*b[1] + u[1]*b[2] + u[2]*b[3]);
  w[5] = -(u[0]*b[0] + u[3]*b[1] + u[2]*b[2] + u[1]*b[3]);
  w[6] = -(u[1]*b[0] + u[0]*b[1] + u[2]*b[2] + u[3]*b[3]);
  w[7] = -(u[1]*b[0] + u[0]*b[1] + u[3]*b[2] + u[2]*b[3]);
  w[8] = -(u[1]*b[0] + u[2]*b[1] + u[0]*b[2] + u[3]*b[3]);
  w[9] = -(u[1]*b[0] + u[2]*b[1] + u[3]*b[2] + u[0]*b[3]);
  w[10]= -(u[1]*b[0] + u[3]*b[1] + u[0]*b[2] + u[2]*b[3]);
  w[11]= -(u[1]*b[0] + u[3]*b[1] + u[2]*b[2] + u[0]*b[3]);
  w[12]= -(u[2]*b[0] + u[0]*b[1] + u[1]*b[2] + u[3]*b[3]);
  w[13]= -(u[2]*b[0] + u[0]*b[1] + u[3]*b[2] + u[1]*b[3]);
  w[14]= -(u[2]*b[0] + u[1]*b[1] + u[0]*b[2] + u[3]*b[3]);
  w[15]= -(u[2]*b[0] + u[1]*b[1] + u[3]*b[2] + u[0]*b[3]);
  w[16]= -(u[2]*b[0] + u[3]*b[1] + u[0]*b[2] + u[2]*b[3]);
  w[17]= -(u[2]*b[0] + u[3]*b[1] + u[2]*b[2] + u[0]*b[3]);
  w[18]= -(u[3]*b[0] + u[0]*b[1] + u[1]*b[2] + u[2]*b[3]);
  w[19]= -(u[3]*b[0] + u[0]*b[1] + u[2]*b[2] + u[1]*b[3]);
  w[20]= -(u[3]*b[0] + u[1]*b[1] + u[0]*b[2] + u[2]*b[3]);
  w[21]= -(u[3]*b[0] + u[1]*b[1] + u[2]*b[2] + u[0]*b[3]);
  w[22]= -(u[3]*b[0] + u[2]*b[1] + u[0]*b[2] + u[1]*b[3]);
  w[23]= -(u[3]*b[0] + u[2]*b[1] + u[1]*b[2] + u[0]*b[3]);
  //*/
};

static inline
double CalcS(double *w,int ind,int imax,int n){
  double tmp = 0e0;

  if(ind == -1){
    for(int i=imax;i<n;i++){
      tmp = sum_log(tmp,w[i]);
    }
    return tmp;
  }
  if(ind >= imax){
    for(int i=imax;i<ind+1;i++){
      tmp = sum_log(tmp,w[i]);
    }
    return tmp;
  }
  if(ind < imax){
    for(int i=imax;i<n;i++){
      tmp = sum_log(tmp,w[i]);
    }
    for(int i=0;i<ind+1;i++){
      tmp = sum_log(tmp,w[i]);
    }
    return tmp;
  }
};

static inline
void CalcTransProb(double *w,double *p,int imax,int n){
  int nf = factorial(n);
  double s_a = CalcS(w,0,imax,nf);
  for(int i=0;i<nf;i++){
    double s_b = CalcS(w,i-1,imax,nf);
    double delta;// = s_a  + log(1. - exp(s_b - s_a) + exp(w[imax] - s_a));
    double tmp[4];
    if(sum_log(s_a,w[imax]) > s_b){
      delta = s_a + log( 1. - exp(s_b-s_a) + exp(w[imax]-s_a) );
    }else{//if w(X_a)+w(X_)-delta < 0, p will be 0
      delta =  -1e10;
    }
    tmp[0] = exp(delta-w[0]);
    if(sum_log(w[0],w[i]) > delta){
      tmp[1]= w[0] + log( 1. + exp(w[i]-w[0]) - exp(delta-w[0]) );
    }else{//if w(X_a)+w(X_)-delta < 0, p will be 0
      tmp[1]= -1e10;
    }
    tmp[1] = exp(tmp[1]-w[0]);
    tmp[2] = exp(w[0]-w[0]);
    tmp[3] = exp(w[i]-w[0]);
    p[i]= max( 0., min( tmp[0], tmp[1], tmp[2], tmp[3] ) );
  }
};


void fx(int n,int k,int *num,int *frag,int **ind,int &count){
  int i,f,j;
  for(i=0;i<n;++i){
    f=frag[i];
    if(f == 0){
      frag[i]=1;
      num[k] = i+1;
      if(k==1){
	for(j=n;j > 0;--j){
	  ind[count/n][count%n] = num[j]-1;
	  count++;
	}
      }else{
	fx(n,k-1,num,frag,ind,count);
      }
      frag[i] = 0;
    }
  }
}

static void GetPermutation(int **ind,int n){
  int num[n],frag[n];
  static int count;
  for(int i=0;i<n;i++){
    num[i]=frag[i]=0;
  }
  count = 0;
  fx(n,n,num,frag,ind,count);
};

//template<uint N>
void REMD::ReplicaPermutation(){// now only for N=4
  //printf("replica permutation method\n");
  const int nreplica = remdinfo->GetNumReplica();
  //const int N = nreplica;
  const int N = 4;
  const int Nfac = factorial(N);
  //printf("N= %d,Nfac= %d,step= %d\n",N,Nfac,step);
  //calc weight factor
  int **indexlist;
  SAFE_MALLOC(indexlist,Nfac*sizeof(int*),int*);
  for(int i=0;i<Nfac;i++){
    SAFE_MALLOC(indexlist[i],N*sizeof(int),int);
  }
  GetPermutation(indexlist,N);
  double w[Nfac],p[Nfac];
  int EvenOrOdd = step%2;

  for(int i=0;i<(nreplica-EvenOrOdd)/N;i++){
    int ind[N];

    for(int n=0;n<N;n++){
      ind[n] = remdinfo->GetIndex(i*N+n+2*EvenOrOdd);
    }
    double t[N];

    for(int n=0;n<N;n++){
      t[n] = remdinfo->GetTemperature(i*N+n+2*EvenOrOdd);
    }
    double u[N];
    for(int n=0;n<N;n++){
      u[n] = pot[i*N+n+2*EvenOrOdd];
    }

    CalcWeight(t,u,w,indexlist,N);
    
    double wmax  = w[0];
    int    iwmax = 0;
    for(int j=0;j<Nfac;j++){
      if(wmax<w[j]){
	wmax = w[j];iwmax=j;
      }
    }
    CalcTransProb(w,p,iwmax,N);
    double sum=0e0,rndm=random_number();
    for(int j=0;j<Nfac;j++){
      //printf("p[%d]= %e,rndm=%e\n",j,p[j],rndm);
      sum += p[j];
      if(sum>rndm){
	remdinfo->Permutate(ind,indexlist[j],t,N,N*i+2*EvenOrOdd);
	break;
      }
    }
    //printf("sum= %e\n",sum);
  }

  for(int i=0;i<nreplica;++i){
    md[i]->ChangeTemperature(remdinfo->GetTemperature(i));
    //md[i]->SetPressure(remdinfo->GetPressure(i));
  }
  for(int i=0;i<nreplica;i++){
    //printf("%d %lf\n",i,md[i]->GetTemperature());
  }
};

void REMD::IncrementHistogram(){
  const int nreplica = remdinfo->GetNumReplica();

  for(int i=0;i<nreplica;++i){
    const int ind = remdinfo->GetIndex(i);
    const D p = pot[ind]/(D)remdinfo->GetNmol();
    const D v = vol[ind]/(D)remdinfo->GetNmol();
    remdinfo->IncrementHistogram(v,p,i);
    const Average a = ave[ind];
    remdinfo->KeepAverages(v,p,i,a);
  }
};

void REMD::OutputBackUp(){
  printf("===OutputBackUp() of %ld REMD step start===\n",step);
  const int nreplica = remdinfo->GetNumReplica();
  const char *output_dir = remdinfo->GetOutputDir();

  for(int i=0;i<nreplica;++i){
    const int ind   = i;
    //const D temp  = md[i]->GetTemperature();
    //const D press = md[i]->GetPressure();
    //iomanager[ind]->UpdateOutputs(md[ind]);
  }
  char indexfile[256];
  sprintf(indexfile,"%s/index.dat",output_dir);
  printf("====making %s====\n",indexfile);
  FILE *fp;
  SAFE_FILEOPEN(fp,indexfile,"w");
  for(int i=0;i<nreplica;++i){
    fprintf(fp,"%4d %4d\n",i,remdinfo->GetIndex(i));
  }
  fclose(fp);
};

void REMD::OutputEnergy(){
  *os_ene << std::setw(8) << step;
  for(int i=0;i<remdinfo->GetNumReplica();++i){
    double p = pot[i];
    if(param.pconstant > 0){
      //p += vol[i]*md[i]->GetPressure();
      p += vol[i]*remdinfo->GetPressure(i);
    }
    *os_ene << " " << std::setw(16) << std::setprecision(16) << p;
  }
  *os_ene << std::endl;
};

void REMD::OutputIndex(){
  FILE *fp_t,*fp_p;
  char *output_dir = remdinfo->GetOutputDir();
  char filename1[256],filename2[256];
  if(step==0){
    printf("===%s/index.dat is output every 100 interval===\n",output_dir);
    sprintf(filename1,"%s/index_temp.dat",output_dir);
    sprintf(filename2,"%s/index_press.dat",output_dir);
    SAFE_FILEOPEN(fp_t,filename1,"w");
    SAFE_FILEOPEN(fp_p,filename2,"w");
  }else{
    sprintf(filename1,"%s/index_temp.dat",output_dir);
    sprintf(filename2,"%s/index_press.dat",output_dir);
    SAFE_FILEOPEN(fp_t,filename1,"a");
    SAFE_FILEOPEN(fp_p,filename2,"a");
  }

  fprintf(fp_t," %7ld",step);
  fprintf(fp_p," %7ld",step);
  for(int i=0;i<remdinfo->GetNumReplica();++i){
    fprintf(fp_t," %lf",remdinfo->GetTemperature(i));
    fprintf(fp_p," %lf",remdinfo->GetPressure(i));
  }
  fprintf(fp_t,"\n");
  fprintf(fp_p,"\n");

  fclose(fp_t);
  fclose(fp_p);
};

void REMD::ExecuteWHAM(){
  const int iter = 1000;
  for(int i=0;i<iter;++i){
    wham->CalcDensState(remdinfo);
    wham->CalcG(remdinfo);
  }
  wham->CalcPhysValue(remdinfo);
  wham->Output(remdinfo);
};

void REMD::OutputBinary(){
  if(remdinfo->GetRank() == 0){
    for(int rep=0;rep<remdinfo->GetNumReplicaGlobal();rep++){
#if 0
      if(ngpus > 0){
#ifdef ENABLE_GPU
	const Average a = mdgpu->GetAverages(remdinfo->GetIndex(rep));
	binary[rep].write((char*)&a.sum,Average::NUM_ELEM*sizeof(double));
#else
	fprintf(stderr,"%s : %s : error: gpu is not enabled!",__FILE__,__LINE__);
	exit(EXIT_FAILURE);
#endif
      }else{
	const Molecules* m = md[remdinfo->GetIndex(rep)];
	const Property prop = m->GetProperty();
	Average a;
	a.sum[Average::LJ]   = prop.lj;
	a.sum[Average::CLMB] = prop.clmb;
	a.sum[Average::WALL] = prop.wall;
	a.sum[Average::TRA]  = prop.tra;
	a.sum[Average::ROT]  = prop.rot;
	a.sum[Average::T]    = prop.T;
	a.sum[Average::P]    = prop.prs[3];
	a.sum[Average::Px]   = prop.prs[0];
	a.sum[Average::Py]   = prop.prs[1];
	a.sum[Average::Pz]   = prop.prs[2];
	a.sum[Average::TOT]  = prop.tot;
	a.sum[Average::DRF]  = prop.ham;
	a.sum[Average::TSTs] = m->tst->s;
	a.sum[Average::TSTv] = m->tst->Ps;;
	a.sum[Average::BSTx] = m->bst->Pv[0];
	a.sum[Average::BSTy] = m->bst->Pv[1];
	a.sum[Average::BSTz] = m->bst->Pv[2];
	binary[rep].write((char*)&a.sum,Average::NUM_ELEM*sizeof(double));
      }
#else
      Average a = ave[remdinfo->GetIndex(rep)];
      binary[rep].write((char*)&a.sum,Average::NUM_ELEM*sizeof(double));
#endif
    }
  }
}

void REMD::OutputDebugInfo(){
  int nmol = remdinfo->GetNumReplica();
  char *output_dir = remdinfo->GetOutputDir();
  char filename[256];

  FILE *fp_pot,*fp_kin,*fp_vol;
  //FILE *fp_pre;
  sprintf(filename,"%s/debug_pot.dat",output_dir);
  SAFE_FILEOPEN(fp_pot,filename,"a");
  sprintf(filename,"%s/debug_kin.dat",output_dir);
  SAFE_FILEOPEN(fp_kin,filename,"a");
  sprintf(filename,"%s/debug_vol.dat",output_dir);
  SAFE_FILEOPEN(fp_vol,filename,"a");
  //sprintf(filename,"%s/debug_pre.dat",output_dir);
  //SAFE_FILEOPEN(fp_pre,filename,"a");

  fprintf(fp_pot," %ld",step);
  for(int i=0;i<nmol;++i) fprintf(fp_pot," %lf",md[i]->GetPotential());
  fprintf(fp_pot,"\n");

  fprintf(fp_kin," %ld",step);
  for(int i=0;i<nmol;++i) fprintf(fp_kin," %lf",md[i]->GetKinetic());
  fprintf(fp_kin,"\n");

  fprintf(fp_vol," %ld",step);
  for(int i=0;i<nmol;++i) fprintf(fp_vol," %lf",md[i]->GetVolume());
  fprintf(fp_vol,"\n");
  /*
  fprintf(fp_pre," %ld",step);
  for(int i=0;i<nmol;++i) fprintf(fp_pre," %lf",md[i]->GetPressureTmp());
  fprintf(fp_pre,"\n");
  //*/
  fclose(fp_pot);
  fclose(fp_kin);
  fclose(fp_vol);
  //fclose(fp_pre);
};

void REMD::GatherResult(){
  const int rank     = remdinfo->GetRank();
  const int nproc    = remdinfo->GetNumProc();
  const int nreplica = remdinfo->GetNumReplica();
#ifdef ENABLE_MPI
  if(rank == 0){
    // gather pot, vol, and ae
    MPI_Request req[nproc];
    MPI_Status  stat[nproc];
    for(int i=1; i<nproc;i++) MPI_Irecv(pot+i*nreplica,nreplica,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&req[i-1]);
    MPI_Waitall(nproc-1,req,stat);
    for(int i=1; i<nproc;i++) MPI_Irecv(vol+i*nreplica,nreplica,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&req[i-1]);
    MPI_Waitall(nproc-1,req,stat);
    for(int i=1; i<nproc;i++) MPI_Irecv((double*)(ave+i*nreplica),nreplica*sizeof(Average)/sizeof(double),MPI_DOUBLE,i,i,MPI_COMM_WORLD,&req[i-1]);
    MPI_Waitall(nproc-1,req,stat);
  }else{
    MPI_Send(pot,nreplica,MPI_DOUBLE,0,rank,MPI_COMM_WORLD);
    MPI_Send(vol,nreplica,MPI_DOUBLE,0,rank,MPI_COMM_WORLD);
    MPI_Send((double*)ave,nreplica*sizeof(Average)/sizeof(double),MPI_DOUBLE,0,rank,MPI_COMM_WORLD);
  }
#endif
}
