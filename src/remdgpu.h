#ifndef H_REMDGPU
#define H_REMDGPU

#include "molecules.h"
#include "io_manager.h"
#include "integrator.h"
#include "remdinfo.h"

#include "typedef.h"
#include <vector_types.h>

#include <vgg_common.h>

//#define DEBUG
//#define TUNING
//#define HIGHACC
#define FASTMATH

//#define DIRECT
#define SWITCHING

#define tid threadIdx.x
#define bid blockIdx.x
#define bdm blockDim.x

#define NGPUMAX 12
#define NTHREADS 192 // for kepler series, a SMX has 192 cuda cores
#define NMOLTYPE 1

//#define NUMSHARED 96 // NUMSHARED must be equal to # of threads!!!
//#define NUMSHARED 128 // NUMSHARED must be equal to # of threads!!!
//#define NUMSHARED 192 // NUMSHARED must be equal to # of threads!!!
//#define NUMSHARED 256 // NUMSHARED must be equal to # of threads!!!
//#define NUMSHARED 768 // NUMSHARED must be equal to # of threads!!!
#define NUMSHARED NTHREADS

#if 1
typedef int     rtype;
typedef int4    rtype4;
typedef double  vtype;
typedef double4 vtype4;
typedef double  ftype;
typedef double4 ftype4;
typedef double  qtype;
typedef double4 qtype4;
typedef double  ptype;
typedef double4 ptype4;
typedef double  ntype;
typedef double4 ntype4;
#else
typedef int    rtype;
typedef int4   rtype4;
typedef float  vtype;
typedef float4 vtype4;
typedef float  ftype;
typedef float4 ftype4;
typedef float  qtype;
typedef float4 qtype4;
typedef float  ptype;
typedef float4 ptype4;
typedef float  ntype;
typedef float4 ntype4;
#endif

typedef struct{
  float x;
  float y;
  float z;
  int   w;
}f3i1;

typedef struct{
  float  lj;
  float  clmb;
  float  wall;
  float  pot;
  float  vol;
  float  tra;
  float  rot;
  float  kin;
  float2 tst;
  float2 bst;
  float  tot;
  float  drf;
  float3 vir;
  float3 pres;
  float  temp;
  float2 t;
  float3 b;
}Result;

class REMDGPU{
private:
  //parameters
  unsigned int  nthreads;
  unsigned int  ngpus;
  unsigned int  gpu_offset;

  unsigned int  nmol;
  unsigned int  natom;
  unsigned int  nsite;
  unsigned int  nreplica;
  float         dt;
  int           mode;

  float         rcut;
  float         alpha;
  unsigned long interval;
  unsigned int  nwave;
  float         rswitch;
  // host pointer
  int4   *gcoor;
  float4 *vel;
  float4 *force;
  qtype4 *angle;
  ptype4 *angvel;
  float4 *torque;

  float3 *tst;
  float4 *bst;

  float3 *L;

  float4 *potvir;
  float2 *cond0;
  float2 *cond1;
  Result *rslt;
  float  *H0;

  int4   *acoor;
  float4 *aforce;

  int3   *kvec;
  float2 *kvecr;

  // host pinter for constant memory
  float3 *iner;
  float3 *atom;
  float3 *atomf;
  float  **epsilon;
  float  **sigma;
  float  **ljc12;
  float  **ljc06;
  float  **lje12;
  float  **lje06;
  float  *charge;

  // device pointer for gpus
  int4   *gcoor_dev[NGPUMAX];
  float4 *vel_dev[NGPUMAX];
  float4 *force_dev[NGPUMAX];
  qtype4 *angle_dev[NGPUMAX];
  ptype4 *angvel_dev[NGPUMAX];
  float4 *torque_dev[NGPUMAX];

  float3 *tst_dev[NGPUMAX];
  float4 *bst_dev[NGPUMAX];

  int4   *acoor_dev[NGPUMAX];
  float4 *aforce_dev[NGPUMAX];
  float2 *kvecr_dev[NGPUMAX];

  int3   *kvc_dev[NGPUMAX];

  float3 *L_dev[NGPUMAX];
  float4 *potvir_dev[NGPUMAX];

  float2 *cond0_dev[NGPUMAX];
  float2 *cond1_dev[NGPUMAX];
  Result *rslt_dev[NGPUMAX];
  float  *H0_dev[NGPUMAX];

  //for confined
  cudaArray *cu_array[NGPUMAX];
  cudaArray *cu_array2[NGPUMAX];
  //cudaArray *switching_array[NGPUMAX];
  int       nfwall;
  float     wall_length;

  void ExecuteKernel();

public:
  REMDGPU(Molecules**, const unsigned int);
  ~REMDGPU();
  void ExecuteSteps(double*,double*);
  void GetBackAllVariablesFromGPU();

  // output
  void OutputCheckPoint(const char*,const int*);
  void OutputCDV(const char*,const unsigned int,const int*);
  void ConvertVariablesToHost(Molecules**);

  // getter for remd
  float GetPotential     (unsigned int rep){return rslt[rep].pot;};
  float GetVolume        (unsigned int rep){return rslt[rep].vol;};
  /*
  float GetTmpPressure   (unsigned int rep){return rslt[rep].pres;};
  float GetTmpTemperature(unsigned int rep){return rslt[rep].temp;};
  float GetLJ            (unsigned int rep){return rslt[rep].lj;};
  float GetCoulomb       (unsigned int rep){return rslt[rep].clmb;};
  float GetWall          (unsigned int rep){return rslt[rep].wall;};
  float GetTra           (unsigned int rep){return rslt[rep].tra;};
  float GetRot           (unsigned int rep){return rslt[rep].rot;};
  float GetKin           (unsigned int rep){return rslt[rep].kin;};
  float GetTot           (unsigned int rep){return rslt[rep].tot;};
  float3 GetVirial       (unsigned int rep){return rslt[rep].vir;};
  //*/
  Average GetAverages    (unsigned int rep){
    Average tmp;
    tmp.sum[Average::LJ]   = rslt[rep].lj;
    tmp.sum[Average::CLMB] = rslt[rep].clmb;
    tmp.sum[Average::WALL] = rslt[rep].wall;

    tmp.sum[Average::VOL]  = rslt[rep].vol;
    tmp.sum[Average::T]    = rslt[rep].temp;
    tmp.sum[Average::Px]   = rslt[rep].pres.x;
    tmp.sum[Average::Py]   = rslt[rep].pres.y;
    tmp.sum[Average::Pz]   = rslt[rep].pres.z;

    tmp.sum[Average::TRA]  = rslt[rep].tra;
    tmp.sum[Average::ROT]  = rslt[rep].rot;
    tmp.sum[Average::TOT]  = rslt[rep].tot;
    tmp.sum[Average::DRF]  = rslt[rep].drf;

    tmp.sum[Average::TSTs] = rslt[rep].t.x;
    tmp.sum[Average::TSTv] = rslt[rep].t.y;

    tmp.sum[Average::BSTx] = rslt[rep].b.x;
    tmp.sum[Average::BSTy] = rslt[rep].b.y;
    tmp.sum[Average::BSTz] = rslt[rep].b.z;

    return tmp;
  };

  //debug function
  void PrintAll(std::ostream&);
  void Diff(Molecules**,std::ostream&);
};

#endif
