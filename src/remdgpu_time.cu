#include "remdgpu.h"
#include "reduce.cuh"

#include "vector.cuh"

#include <string>

const int MAXNATOM = 5;

__device__ __constant__ float3 inertia[NMOLTYPE];
__device__ __constant__ float3 atm[MAXNATOM*NMOLTYPE];
__device__ __constant__ float  eps[NMOLTYPE*MAXNATOM+1][NMOLTYPE*MAXNATOM+1];
__device__ __constant__ float  sgm02[NMOLTYPE*MAXNATOM+1][NMOLTYPE*MAXNATOM+1];
__device__ __constant__ float  chg[NMOLTYPE*MAXNATOM+1];

texture<float2,cudaTextureType1D,cudaReadModeElementType> tex_wall;
texture<float2,cudaTextureType1D,cudaReadModeElementType> tex_switch;

#include "sort.cuh"

#ifdef TUNING
#include "4site.cuh"
#endif

class GPUTimer{
  cudaEvent_t begin,end;
  float kernel;
  float memcpy;
public:
  GPUTimer(){
    kernel = memcpy = 0.f;
    cudaEventCreate(&begin);
    cudaEventCreate(&end);
  }
  ~GPUTimer(){
    const float total = kernel + memcpy;
    printf("kernel:  %e usec, %4.2f%\n",kernel,kernel/total*100);
    printf("memcpy:  %e usec, %4.2f%\n",memcpy,memcpy/total*100);
  }
  void start(){
    cudaEventRecord(begin,0);
  }
  void stop(std::string type){
    cudaEventRecord(end,0);
    cudaEventSynchronize(end);
    cudaThreadSynchronize();

    float elapsed_time;
    cudaEventElapsedTime(&elapsed_time,begin,end);
    if(type == "kernel")      kernel += elapsed_time;
    else if(type == "memcpy") memcpy += elapsed_time;
    else fprintf(stderr,"error: wrong type set in time profiler in REMDGPU\n");
  }
};

class KernelTimer{
  cudaEvent_t begin,end;
  float velscale;
  float sort;
  float remd;
public:
  KernelTimer(){
    velscale = remd = sort = 0.f;
    cudaEventCreate(&begin);
    cudaEventCreate(&end);
  }
  ~KernelTimer(){
    const float total = remd + sort + velscale;
    printf("remd:     %e usec, %4.2f%\n",remd,remd/total*100);
    printf("sort:     %e usec, %4.2f%\n",sort,sort/total*100);
    printf("velscale: %e usec, %4.2f%\n",velscale,velscale/total*100);

    cudaEventDestroy(begin);
    cudaEventDestroy(end);
  }
  void start(){
    cudaEventRecord(begin,0);
  }
  void stop(std::string type){
    cudaEventRecord(end,0);
    cudaEventSynchronize(end);
    cudaThreadSynchronize();

    float elapsed_time;
    cudaEventElapsedTime(&elapsed_time,begin,end);
    if(type == "remd")          remd   += elapsed_time;
    else if(type == "sort")     sort     += elapsed_time;
    else if(type == "velscale") velscale += elapsed_time;
    else fprintf(stderr,"error: wrong type set in time profiler in KernelTimer\n");
  }
};

template<int>
__global__
void vel_scale
(
 int4   *,
 float4 *,
 qtype4 *,
 ptype4 *,
 float3 *,
 float3 *,
 float2 *,
 const int
 );

template<int>
__global__
void remd_kernel
(int4*,
 float4*,
 float4*,
 qtype4*,
 ptype4*,
 float4*,
 int4*,
 float4*,
 float2*,
 int3*,
 float4*,
 float3*,
 float4*,
 float3*,
 float2*,
 float2*,
 Result*,
 float*,
 const unsigned int,
 const unsigned int,
 const float,
 const float,
 const float,
 const unsigned int,
 const float,
 const unsigned long,
 const float
 );

template<int>
__global__
void calc_direct_potvir
(int4*,
 float4*,
 float4*,
 qtype4*,
 ptype4*,
 float4*,
 int4*,
 float4*,
 float2*,
 int3*,
 float4*,
 float3*,
 float4*,
 float3*,
 float2*,
 float2*,
 float4*,
 float*,
 const unsigned int,
 const float,
 const float,
 const float,
 const unsigned int,
 const unsigned long,
 const float
 );


REMDGPU::REMDGPU(Molecules** mlcls, const unsigned int _nreplica){
  nthreads   = mlcls[0]->param.nthreads;
  ngpus      = mlcls[0]->param.ngpus;
  gpu_offset = mlcls[0]->param.gpu_offset;

  nreplica = mlcls[0]->param.nreplica;
  nmol     = mlcls[0]->nmol;
  nsite    = mlcls[0]->param.nsite;
  /*
  if(nmol%32!=0){
    std::cerr << "error: unsupported number of molecules!" << std::endl;
    exit(EXIT_FAILURE);
  }
  //*/
  nwave    = mlcls[0]->fclcltr->nwave;
  mode     = mlcls[0]->mode;

  dt       = mlcls[0]->dt;

  rcut     = mlcls[0]->param.rcut;
  alpha    = mlcls[0]->param.alpha;
  interval = mlcls[0]->param.interval;
  rswitch  = mlcls[0]->param.rswitch;

  if(mlcls[0]->param.confined > 0){
    nfwall = mlcls[0]->param.nfwall;
    wall_length = mlcls[0]->param.wall_length;
  }

  natom = nsite*nmol;
  if(natom!=mlcls[0]->natom){
    std::cerr << "error: natom in REMDGPU (=" << natom << ") is not equal to that of Molecules (=" << mlcls[0]->natom << ")!" << std::endl;
    std::cerr << "error: You should change the value of nsite in input file" << std::endl;
    exit(EXIT_FAILURE);
  }

  // malloc host variables
#define HostMalloc(type,ptr,size)					\
  if( (ptr = (type)malloc(size)) == NULL){				\
    fprintf(stderr,"error: malloc "#ptr" failed\n");			\
    exit(EXIT_FAILURE);							\
  }

  HostMalloc(int4*,   gcoor, nreplica * nmol  * sizeof(int4)  );
  HostMalloc(float4*, vel,   nreplica * nmol  * sizeof(float4));
  HostMalloc(float4*, force, nreplica * nmol  * sizeof(float4));
  HostMalloc(qtype4*, angle, nreplica * nmol  * sizeof(qtype4));
  HostMalloc(ptype4*, angvel,nreplica * nmol  * sizeof(ptype4));
  HostMalloc(float4*, torque,nreplica * nmol  * sizeof(float4));

  HostMalloc(float3*, tst,   nreplica *         sizeof(float3));
  HostMalloc(float4*, bst,   nreplica *         sizeof(float4));

  HostMalloc(float3*, L,     nreplica *         sizeof(float3));

  HostMalloc(float4*, potvir,nreplica *         sizeof(float4));
  HostMalloc(float2*, cond0, nreplica *         sizeof(float2));
  HostMalloc(float2*, cond1, nreplica *         sizeof(float2));
  HostMalloc(Result*, rslt,  nreplica *         sizeof(Result));
  HostMalloc(float*,    H0,  nreplica *         sizeof(float) );

  HostMalloc(int4*,   acoor, nreplica * natom * sizeof(int4) );
  HostMalloc(float4*, aforce,nreplica * natom * sizeof(float4));

  HostMalloc(int3*,   kvec,  nreplica * nwave * sizeof(int3)  );
  HostMalloc(float2*, kvecr, nreplica * nwave * sizeof(float2));

  //set switching function to host
  float2 *switching;
  HostMalloc(float2*, switching, nfwall * sizeof(float2));
  const float rlc = rswitch - rcut;
  const float dr  = (rcut-rswitch) / (float)(nfwall-1);
  const float coef_sw = 1.f / (rlc*rlc*rlc*rlc*rlc);
  for(int s=0;s<nfwall;s++){
    const float rl = dr * s;
    const float rc = rl + rlc;
    switching[s].x = coef_sw * rc*rc*rc*(10.f*rl*rl - 5.f*rl*rc + rc*rc);
    switching[s].y = coef_sw * 30.f*rc*rc*rl*rl;
  }

  // set host variables
  for(unsigned int rep=0;rep<nreplica;rep++){
    for(unsigned int i=0;i<nmol;i++){
      gcoor[nmol*rep + i].x  = (int)((double)IntMax * (mlcls[rep]->mlcl[i].r[0] - 0.5));
      gcoor[nmol*rep + i].y  = (int)((double)IntMax * (mlcls[rep]->mlcl[i].r[1] - 0.5));
      gcoor[nmol*rep + i].z  = (int)((double)IntMax * (mlcls[rep]->mlcl[i].r[2] - 0.5));
      gcoor[nmol*rep + i].w  = mlcls[rep]->mlcl[i].type;

      vel[nmol*rep + i].x    = (float)mlcls[rep]->mlcl[i].v[0];
      vel[nmol*rep + i].y    = (float)mlcls[rep]->mlcl[i].v[1];
      vel[nmol*rep + i].z    = (float)mlcls[rep]->mlcl[i].v[2];
      vel[nmol*rep + i].w    = (float)mlcls[rep]->mlcl[i].m;

      force[nmol*rep + i].x  = (float)mlcls[rep]->mlcl[i].f[0];
      force[nmol*rep + i].y  = (float)mlcls[rep]->mlcl[i].f[1];
      force[nmol*rep + i].z  = (float)mlcls[rep]->mlcl[i].f[2];

      angle[nmol*rep + i].x  = (qtype)mlcls[rep]->mlcl[i].q[0];
      angle[nmol*rep + i].y  = (qtype)mlcls[rep]->mlcl[i].q[1];
      angle[nmol*rep + i].z  = (qtype)mlcls[rep]->mlcl[i].q[2];
      angle[nmol*rep + i].w  = (qtype)mlcls[rep]->mlcl[i].q[3];

      angvel[nmol*rep + i].x = (ptype)mlcls[rep]->mlcl[i].p[0];
      angvel[nmol*rep + i].y = (ptype)mlcls[rep]->mlcl[i].p[1];
      angvel[nmol*rep + i].z = (ptype)mlcls[rep]->mlcl[i].p[2];
      angvel[nmol*rep + i].w = (ptype)mlcls[rep]->mlcl[i].p[3];

      torque[nmol*rep + i].x = (float)mlcls[rep]->mlcl[i].n[0];
      torque[nmol*rep + i].y = (float)mlcls[rep]->mlcl[i].n[1];
      torque[nmol*rep + i].z = (float)mlcls[rep]->mlcl[i].n[2];
      torque[nmol*rep + i].w = (float)mlcls[rep]->mlcl[i].n[3];
    }

    tst[rep].x = mlcls[rep]->tst->s;
    tst[rep].y = mlcls[rep]->tst->Ps;
    tst[rep].z = 0.5f / mlcls[rep]->tst->Q;

    bst[rep].x = mlcls[rep]->bst->Pv[0];
    bst[rep].y = mlcls[rep]->bst->Pv[1];
    bst[rep].z = mlcls[rep]->bst->Pv[2];
    bst[rep].w = 0.5f / mlcls[rep]->bst->W;

    L[rep].x = mlcls[rep]->L[0];
    L[rep].y = mlcls[rep]->L[1];
    L[rep].z = mlcls[rep]->L[2];

    potvir[rep].x = mlcls[rep]->prop.vir[0];
    potvir[rep].y = mlcls[rep]->prop.vir[1];
    potvir[rep].z = mlcls[rep]->prop.vir[2];
    potvir[rep].w = mlcls[rep]->prop.pot;

    cond0[rep].x = cond1[rep].x = mlcls[rep]->T;
    cond0[rep].y = cond1[rep].y = mlcls[rep]->P;

    rslt[rep].lj  = mlcls[rep]->prop.pot;
    H0[rep]      = mlcls[rep]->prop.H0;
    std::cout << rep <<  " H0:" << H0[rep] << std::endl;

    for(unsigned int i=0;i<nwave;i++){
      kvec[nwave*rep+i].x = mlcls[rep]->fclcltr->kvec[i][0];
      kvec[nwave*rep+i].y = mlcls[rep]->fclcltr->kvec[i][1];
      kvec[nwave*rep+i].z = mlcls[rep]->fclcltr->kvec[i][2];
    }
  }

  //set inertia and atom coordinate and force in body space
  iner  = (float3*)malloc(NMOLTYPE * sizeof(float3));
  atom  = (float3*)malloc(NMOLTYPE * nsite * sizeof(float3));
  for(unsigned int i=0;i<mlcls[0]->mlist.size();i++){
    iner[i].x = mlcls[0]->mlist[i].i[0];
    iner[i].y = mlcls[0]->mlist[i].i[1];
    iner[i].z = mlcls[0]->mlist[i].i[2];
    for(unsigned int j=0;j<mlcls[0]->mlist[i].a.size();j++){
      atom[nsite*i+j].x = mlcls[0]->mlist[i].a[j].r[0];
      atom[nsite*i+j].y = mlcls[0]->mlist[i].a[j].r[1];
      atom[nsite*i+j].z = mlcls[0]->mlist[i].a[j].r[2];
    }
  }

  // set force parameter
  epsilon = (float**)malloc((NMOLTYPE*nsite+1)*sizeof(float*));
  sigma   = (float**)malloc((NMOLTYPE*nsite+1)*sizeof(float*));
  for(int i=0;i<NMOLTYPE*nsite+1;i++){
    epsilon[i] = (float*)malloc((NMOLTYPE*nsite+1)*sizeof(float));
    sigma[i]   = (float*)malloc((NMOLTYPE*nsite+1)*sizeof(float));
  }
  charge = (float*)malloc((NMOLTYPE*nsite+1)*sizeof(float));

  unsigned int count_i = 0,count_j;
  for(unsigned int ii=0;ii<mlcls[0]->mlist.size();ii++){
    for(unsigned int ij=0;ij<mlcls[0]->mlist[ii].a.size();ij++){
      count_j = 0;
      for(unsigned int ji=0;ji<mlcls[0]->mlist.size();ji++){
	for(unsigned int jj=0;jj<mlcls[0]->mlist[ji].a.size();jj++){
	  sigma[count_i][count_j] = 0.5f*(mlcls[0]->mlist[ii].a[ij].s + mlcls[0]->mlist[ji].a[jj].s);
	  sigma[count_i][count_j] = sigma[count_i][count_j]*sigma[count_i][count_j];
	  epsilon[count_i][count_j] = 4.f * sqrt(mlcls[0]->mlist[ii].a[ij].e * mlcls[0]->mlist[ji].a[jj].e);
	  count_j++;
      }}
      charge[count_i] = mlcls[0]->mlist[ii].a[ij].q;
      count_i++;
  }}
  for(unsigned int i=0;i<NMOLTYPE*nsite+1;i++){
    sigma[NMOLTYPE*nsite][i]   = 0.f;
    sigma[i][NMOLTYPE*nsite]   = 0.f;
    epsilon[NMOLTYPE*nsite][i] = 0.f;
    epsilon[i][NMOLTYPE*nsite] = 0.f;
  }
  charge[NMOLTYPE*nsite] = 0.f;

  //set wall interaction for confined system
  if(mlcls[0]->param.confined > 0){
    // wall interaction for confined water
    float2 fwall[nfwall];
    for(int i=0;i<nfwall;i++){
      fwall[i].x = mlcls[0]->fclcltr->fwall[i][0];
      fwall[i].y = mlcls[0]->fclcltr->fwall[i][1];
    }

    // set to texture memory
    for(int dev=0;dev<ngpus;dev++){
      cudaSetDevice(dev+gpu_offset);

      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float2>();
      cudaMallocArray(&cu_array[dev], &channelDesc, nfwall);
      cudaMemcpyToArray(cu_array[dev], 0, 0, fwall, sizeof(float2)*nfwall, cudaMemcpyHostToDevice);
      tex_wall.normalized = true;
      tex_wall.addressMode[0] = cudaAddressModeClamp;
      tex_wall.filterMode = cudaFilterModeLinear;
      cudaBindTextureToArray(tex_wall, cu_array[dev], channelDesc);
      //cudaBindTexture(0, tex_wall, cu_array[dev], channelDesc, sizeof(float2)*nfwall);
    }

    for(int dev=0;dev<ngpus;dev++){
      cudaSetDevice(dev+gpu_offset);

      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float2>();
      cudaMallocArray(&cu_array2[dev], &channelDesc, nfwall);
      cudaMemcpyToArray(cu_array2[dev], 0, 0, switching, sizeof(float2)*nfwall, cudaMemcpyHostToDevice);
      tex_switch.normalized = true;
      tex_switch.addressMode[0] = cudaAddressModeClamp;
      tex_switch.filterMode = cudaFilterModeLinear;
      cudaBindTextureToArray(tex_switch, cu_array2[dev], channelDesc);
    }
  }

  std::cout << "malloc device memory and memcpy" << std::endl;
  for(unsigned int dev=0;dev<ngpus;dev++){
    cudaSetDevice(dev+gpu_offset);
    
    cudaMalloc((void**) &gcoor_dev[dev],  nmol*nreplica*sizeof(  int4)/ngpus);
    cudaMalloc((void**)   &vel_dev[dev],  nmol*nreplica*sizeof(float4)/ngpus);
    cudaMalloc((void**) &force_dev[dev],  nmol*nreplica*sizeof(float4)/ngpus);
    cudaMalloc((void**) &angle_dev[dev],  nmol*nreplica*sizeof(qtype4)/ngpus);
    cudaMalloc((void**)&angvel_dev[dev],  nmol*nreplica*sizeof(ptype4)/ngpus);
    cudaMalloc((void**)&torque_dev[dev],  nmol*nreplica*sizeof(float4)/ngpus);
    cudaMalloc((void**) &acoor_dev[dev], natom*nreplica*sizeof(  int4)/ngpus);
    cudaMalloc((void**)&aforce_dev[dev], natom*nreplica*sizeof(float4)/ngpus);
    cudaMalloc((void**)   &tst_dev[dev],       nreplica*sizeof(float3)/ngpus);
    cudaMalloc((void**)     &L_dev[dev],       nreplica*sizeof(float3)/ngpus);
    cudaMalloc((void**)   &bst_dev[dev],       nreplica*sizeof(float4)/ngpus);
    cudaMalloc((void**)&potvir_dev[dev],       nreplica*sizeof(float4)/ngpus);
    cudaMalloc((void**) &cond0_dev[dev],       nreplica*sizeof(float2)/ngpus);
    cudaMalloc((void**) &cond1_dev[dev],       nreplica*sizeof(float2)/ngpus);
    cudaMalloc((void**)  &rslt_dev[dev],       nreplica*sizeof(Result)/ngpus);
    cudaMalloc((void**)    &H0_dev[dev],       nreplica*sizeof( float)/ngpus);
    cudaMalloc((void**)   &kvc_dev[dev], nwave*nreplica*sizeof(  int3)/ngpus);
    cudaMalloc((void**) &kvecr_dev[dev], nwave*nreplica*sizeof(float2)/ngpus);

    cudaMemcpy( gcoor_dev[dev], gcoor + dev *  nmol*nreplica/ngpus,  nmol*nreplica*sizeof(  int4)/ngpus, cudaMemcpyHostToDevice);
    cudaMemcpy(   vel_dev[dev],   vel + dev *  nmol*nreplica/ngpus,  nmol*nreplica*sizeof(float4)/ngpus, cudaMemcpyHostToDevice);
    cudaMemcpy( force_dev[dev], force + dev *  nmol*nreplica/ngpus,  nmol*nreplica*sizeof(float4)/ngpus, cudaMemcpyHostToDevice);
    cudaMemcpy( angle_dev[dev], angle + dev *  nmol*nreplica/ngpus,  nmol*nreplica*sizeof(qtype4)/ngpus, cudaMemcpyHostToDevice);
    cudaMemcpy(angvel_dev[dev],angvel + dev *  nmol*nreplica/ngpus,  nmol*nreplica*sizeof(ptype4)/ngpus, cudaMemcpyHostToDevice);
    cudaMemcpy(torque_dev[dev],torque + dev *  nmol*nreplica/ngpus,  nmol*nreplica*sizeof(float4)/ngpus, cudaMemcpyHostToDevice);
    cudaMemcpy(   kvc_dev[dev],  kvec + dev * nwave*nreplica/ngpus, nwave*nreplica*sizeof(  int3)/ngpus, cudaMemcpyHostToDevice);
    cudaMemcpy(   tst_dev[dev],   tst + dev *       nreplica/ngpus,       nreplica*sizeof(float3)/ngpus, cudaMemcpyHostToDevice);
    cudaMemcpy(   bst_dev[dev],   bst + dev *       nreplica/ngpus,       nreplica*sizeof(float4)/ngpus, cudaMemcpyHostToDevice);
    cudaMemcpy(     L_dev[dev],     L + dev *       nreplica/ngpus,       nreplica*sizeof(float3)/ngpus, cudaMemcpyHostToDevice);
    cudaMemcpy(potvir_dev[dev],potvir + dev *       nreplica/ngpus,       nreplica*sizeof(float4)/ngpus, cudaMemcpyHostToDevice);
    cudaMemcpy( cond0_dev[dev], cond0 + dev *       nreplica/ngpus,       nreplica*sizeof(float2)/ngpus, cudaMemcpyHostToDevice);
    cudaMemcpy( cond1_dev[dev], cond1 + dev *       nreplica/ngpus,       nreplica*sizeof(float2)/ngpus, cudaMemcpyHostToDevice);
    cudaMemcpy(  rslt_dev[dev],  rslt + dev *       nreplica/ngpus,       nreplica*sizeof(Result)/ngpus, cudaMemcpyHostToDevice);
    cudaMemcpy(    H0_dev[dev],    H0 + dev *       nreplica/ngpus,       nreplica*sizeof( float)/ngpus, cudaMemcpyHostToDevice);

    cudaMemcpyToSymbol(inertia, iner,         NMOLTYPE*sizeof(float3));
    cudaMemcpyToSymbol(    atm, atom,nsite*NMOLTYPE*sizeof(float3));
    for(unsigned int i=0;i<NMOLTYPE*nsite+1;i++){
      cudaMemcpyToSymbol(sgm02[i],sigma[i],  (NMOLTYPE*nsite+1)*sizeof(float));
      cudaMemcpyToSymbol(eps[i],epsilon[i],(NMOLTYPE*nsite+1)*sizeof(float));
    }
    cudaMemcpyToSymbol(chg,charge,(NMOLTYPE*nsite+1)*sizeof(float));
  }

  // for initial force calculation and hamiltonian
  std::cout << "calc initial force" << std::endl;
  float dt_tmp = dt;
  dt = 0.f;
  ExecuteKernel();
  for(unsigned int dev=0;dev<ngpus;dev++){
    cudaSetDevice(dev+gpu_offset);
    cudaMemcpy(rslt + dev*nreplica/ngpus,rslt_dev[dev],nreplica*sizeof(Result)/ngpus, cudaMemcpyDeviceToHost);
  }
  dt = dt_tmp;

  for(int rep=0;rep<nreplica;rep++){
    H0[rep] = rslt[rep].tot;
  }
  for(unsigned int dev=0;dev<ngpus;dev++){
    cudaSetDevice(dev+gpu_offset);
    cudaMemcpy(    H0_dev[dev],    H0 + dev *       nreplica/ngpus,       nreplica*sizeof( float)/ngpus, cudaMemcpyHostToDevice);
  }

#ifdef DIRECT
  std::cout << "direct coulomb summation begin" << std::endl;
  for(int dev=0;dev<ngpus;dev++){
    cudaSetDevice(dev+gpu_offset);
    calc_direct_potvir<NPT1D><<<nreplica/ngpus,nmol>>>
      ( gcoor_dev[dev],    vel_dev[dev],  force_dev[dev],
	angle_dev[dev], angvel_dev[dev], torque_dev[dev],
	acoor_dev[dev], aforce_dev[dev],
	kvecr_dev[dev],    kvc_dev[dev], potvir_dev[dev],
	tst_dev[dev],    bst_dev[dev],      L_dev[dev],
	cond0_dev[dev],  cond1_dev[dev],   rslt_dev[dev],	H0_dev[dev],
	nmol, dt, rcut*rcut, alpha, nwave, interval,
	wall_length
	);
      }
  for(unsigned int dev=0;dev<ngpus;dev++){
    cudaSetDevice(dev+gpu_offset);
    cudaMemcpy(rslt + dev*nreplica/ngpus,rslt_dev[dev],nreplica*sizeof(Result)/ngpus, cudaMemcpyDeviceToHost);
  }
  for(int rep=0;rep<nreplica;rep++) std::cout << "rep " << rep << ": " << rslt[rep].x << std::endl;
  std::cout << "direct coulomb summation end" << std::endl;
#endif

}

REMDGPU::~REMDGPU(){
  //*
  for(unsigned int dev=0;dev<ngpus;dev++){
    cudaSetDevice(dev+gpu_offset);

    cudaMemcpy( gcoor + dev *  nmol*nreplica/ngpus, gcoor_dev[dev],  nmol*nreplica*sizeof(  int4)/ngpus, cudaMemcpyDeviceToHost);
    cudaMemcpy(   vel + dev *  nmol*nreplica/ngpus,   vel_dev[dev],  nmol*nreplica*sizeof(float4)/ngpus, cudaMemcpyDeviceToHost);
    cudaMemcpy( force + dev *  nmol*nreplica/ngpus, force_dev[dev],  nmol*nreplica*sizeof(float4)/ngpus, cudaMemcpyDeviceToHost);
    cudaMemcpy( angle + dev *  nmol*nreplica/ngpus, angle_dev[dev],  nmol*nreplica*sizeof(qtype4)/ngpus, cudaMemcpyDeviceToHost);
    cudaMemcpy(angvel + dev *  nmol*nreplica/ngpus,angvel_dev[dev],  nmol*nreplica*sizeof(ptype4)/ngpus, cudaMemcpyDeviceToHost);
    cudaMemcpy(torque + dev *  nmol*nreplica/ngpus,torque_dev[dev],  nmol*nreplica*sizeof(float4)/ngpus, cudaMemcpyDeviceToHost);
    cudaMemcpy(   tst + dev *       nreplica/ngpus,   tst_dev[dev],       nreplica*sizeof(float3)/ngpus, cudaMemcpyDeviceToHost);
    cudaMemcpy(   bst + dev *       nreplica/ngpus,   bst_dev[dev],       nreplica*sizeof(float4)/ngpus, cudaMemcpyDeviceToHost);
    cudaMemcpy(     L + dev *       nreplica/ngpus,     L_dev[dev],       nreplica*sizeof(float4)/ngpus, cudaMemcpyDeviceToHost);
    cudaMemcpy( acoor + dev * natom*nreplica/ngpus, acoor_dev[dev], natom*nreplica*sizeof(  int4)/ngpus, cudaMemcpyDeviceToHost);
    cudaMemcpy(aforce + dev * natom*nreplica/ngpus,aforce_dev[dev], natom*nreplica*sizeof(float4)/ngpus, cudaMemcpyDeviceToHost);
    cudaMemcpy( kvecr + dev * nwave*nreplica/ngpus, kvecr_dev[dev], nwave*nreplica*sizeof(float4)/ngpus, cudaMemcpyDeviceToHost);
    cudaMemcpy(  kvec + dev * nwave*nreplica/ngpus,   kvc_dev[dev], nwave*nreplica*sizeof(  int3)/ngpus, cudaMemcpyDeviceToHost);
  //std::cout << "(gpu)recvec 0: " << kvecr[0].x << " " << kvecr[0].y << std::endl;

    cudaMemcpy(cond0 + dev * nreplica/ngpus,cond0_dev[dev], nreplica*sizeof(float2)/ngpus, cudaMemcpyDeviceToHost);
    cudaMemcpy(cond1 + dev * nreplica/ngpus,cond1_dev[dev], nreplica*sizeof(float2)/ngpus, cudaMemcpyDeviceToHost);
    cudaMemcpy( rslt + dev * nreplica/ngpus, rslt_dev[dev], nreplica*sizeof(Result)/ngpus, cudaMemcpyDeviceToHost);
  }
  //*/
  //PrintAll(std::cout);

  free( gcoor);
  free(   vel);
  free( force);
  free( angle);
  free(angvel);
  free(torque);
  free(   tst);
  free(   bst);
  free( acoor);
  free(aforce);
  free( kvecr);
  free(  kvec);
  free(     L);
  free( cond0);
  free( cond1);
  free(  rslt);
  for(unsigned int dev=0;dev<ngpus;dev++){
    cudaSetDevice(dev+gpu_offset);

    cudaFree( gcoor_dev[dev]);
    cudaFree(   vel_dev[dev]);
    cudaFree( force_dev[dev]);
    cudaFree( angle_dev[dev]);
    cudaFree(angvel_dev[dev]);
    cudaFree(torque_dev[dev]);
    cudaFree(   tst_dev[dev]);
    cudaFree(   bst_dev[dev]);
    cudaFree( acoor_dev[dev]);
    cudaFree(aforce_dev[dev]);
    cudaFree( kvecr_dev[dev]);
    cudaFree(   kvc_dev[dev]);
    cudaFree(     L_dev[dev]);
    cudaFree(potvir_dev[dev]);
    cudaFree( cond0_dev[dev]);
    cudaFree( cond1_dev[dev]);
    cudaFree(  rslt_dev[dev]);
  }
}

void REMDGPU::ExecuteKernel(){
  static KernelTimer timer;
#define KERNEL(x)							\
  case(x):								\
    if(((x>>TSHIFT)&MASK)==2){						\
      timer.start();							\
    vel_scale< x ><<<nreplica/ngpus,NTHREADS>>>				\
	( gcoor_dev[dev],vel_dev[dev],					\
	  angle_dev[dev],angvel_dev[dev],				\
	  tst_dev[dev],L_dev[dev],					\
	  cond0_dev[dev],nmol						\
	  );								\
      timer.stop("velscale");						\
    }									\
    timer.start();							\
    sort<<<nreplica/ngpus,nmol>>>					\
      (gcoor_dev[dev],    vel_dev[dev],  force_dev[dev],		\
       angle_dev[dev], angvel_dev[dev], torque_dev[dev],512);		\
    timer.stop("sort");							\
    timer.start();							\
    remd_kernel< x ><<<nreplica/ngpus,NTHREADS>>>			\
      ( gcoor_dev[dev],    vel_dev[dev],  force_dev[dev],		\
	angle_dev[dev], angvel_dev[dev], torque_dev[dev],		\
	acoor_dev[dev], aforce_dev[dev],				\
	kvecr_dev[dev], kvc_dev[dev], potvir_dev[dev],			\
	tst_dev[dev],   bst_dev[dev],  L_dev[dev],			\
	cond0_dev[dev], cond1_dev[dev],rslt_dev[dev],H0_dev[dev],	\
	nmol, nsite, dt, rcut*rcut, alpha, nwave, rswitch, interval,	\
	wall_length							\
	);								\
    timer.stop("remd");							\
    break;

#pragma omp parallel for num_threads(nthreads)
  for(unsigned int dev=0;dev<ngpus;dev++){
    cudaSetDevice(dev+gpu_offset);
    switch(mode){
#ifndef TUNING
      KERNEL(NVE);
      KERNEL(NVT);
      KERNEL(NVTSCALE);
      KERNEL(NPH);
      KERNEL(NPT);
      KERNEL(NPTSCALE);
      KERNEL(NAxyPzH);
      KERNEL(NAxyPzT);
      KERNEL(NAxyPzTSCALE);
      KERNEL(NPxyLzH);
      KERNEL(NPxyLzT);
      KERNEL(NPxyLzTSCALE);
#endif
      KERNEL(NVE1D);
      KERNEL(NVT1D);
      KERNEL(NVTSCALE1D);
      KERNEL(NPH1D);
      KERNEL(NPT1D);
      KERNEL(NPTSCALE1D);
      KERNEL(NAxyPzH1D);
      KERNEL(NAxyPzT1D);
      KERNEL(NAxyPzTSCALE1D);
      KERNEL(NPxyLzH1D);
      KERNEL(NPxyLzT1D);
      KERNEL(NPxyLzTSCALE1D);
      /*
      KERNEL(NVE2D);
      KERNEL(NVT2D);
      KERNEL(NVTSCALE2D);
      KERNEL(NPH2D);
      KERNEL(NPT2D);
      KERNEL(NPTSCALE2D);
      KERNEL(NAxyPzH2D);
      KERNEL(NAxyPzT2D);
      KERNEL(NAxyPzTSCALE2D);
      KERNEL(NPxyLzH2D);
      KERNEL(NPxyLzT2D);
      KERNEL(NPxyLzTSCALE2D);
      //*/
    default:
      std::cerr << "error: undefined ensemble or confined dimention" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
//cudaDeviceSynchronize();
#undef KERNEL
}

void REMDGPU::ExecuteSteps(double* temperature,double *pressure){
  static GPUTimer timer;
  timer.start();
  for(int rep=0;rep<nreplica;rep++){
    cond0[rep].x = (float)temperature[rep];
    cond0[rep].y = (float)pressure[rep];
  }
#pragma omp parallel for num_threads(nthreads)
  for(unsigned int dev=0;dev<ngpus;dev++){
    cudaSetDevice(dev+gpu_offset);
    cudaMemcpy(cond0_dev[dev],cond0 + dev*nreplica/ngpus,nreplica*sizeof(float2)/ngpus,cudaMemcpyHostToDevice);
  }
  timer.stop("memcpy");

  timer.start();
  ExecuteKernel();
  timer.stop("kernel");

  timer.start();
#pragma omp parallel for num_threads(nthreads)
  for(unsigned int dev=0;dev<ngpus;dev++){
    cudaSetDevice(dev+gpu_offset);
    cudaMemcpy(rslt + dev*nreplica/ngpus,rslt_dev[dev],nreplica*sizeof(Result)/ngpus, cudaMemcpyDeviceToHost);
  }
  timer.stop("memcpy");
}

void REMDGPU::GetBackAllVariablesFromGPU(){
  for(int dev=0;dev<ngpus;dev++){
    cudaSetDevice(dev+gpu_offset);
    cudaMemcpy( gcoor + dev* nmol*nreplica/ngpus, gcoor_dev[dev],  nmol*nreplica*sizeof(  int4)/ngpus, cudaMemcpyDeviceToHost);
    cudaMemcpy(   vel + dev* nmol*nreplica/ngpus,   vel_dev[dev],  nmol*nreplica*sizeof(float4)/ngpus, cudaMemcpyDeviceToHost);
    cudaMemcpy( force + dev* nmol*nreplica/ngpus, force_dev[dev],  nmol*nreplica*sizeof(float4)/ngpus, cudaMemcpyDeviceToHost);
    cudaMemcpy( angle + dev* nmol*nreplica/ngpus, angle_dev[dev],  nmol*nreplica*sizeof(qtype4)/ngpus, cudaMemcpyDeviceToHost);
    cudaMemcpy(angvel + dev* nmol*nreplica/ngpus,angvel_dev[dev],  nmol*nreplica*sizeof(ptype4)/ngpus, cudaMemcpyDeviceToHost);
    cudaMemcpy(torque + dev* nmol*nreplica/ngpus,torque_dev[dev],  nmol*nreplica*sizeof(float4)/ngpus, cudaMemcpyDeviceToHost);
    cudaMemcpy( acoor + dev*natom*nreplica/ngpus, acoor_dev[dev], natom*nreplica*sizeof(  int4)/ngpus, cudaMemcpyDeviceToHost);
    cudaMemcpy(aforce + dev*natom*nreplica/ngpus,aforce_dev[dev], natom*nreplica*sizeof(float4)/ngpus, cudaMemcpyDeviceToHost);
    cudaMemcpy(   tst + dev*      nreplica/ngpus,   tst_dev[dev],       nreplica*sizeof(float3)/ngpus, cudaMemcpyDeviceToHost);
    cudaMemcpy(   bst + dev*      nreplica/ngpus,   bst_dev[dev],       nreplica*sizeof(float4)/ngpus, cudaMemcpyDeviceToHost);
    cudaMemcpy(     L + dev*      nreplica/ngpus,     L_dev[dev],       nreplica*sizeof(float3)/ngpus, cudaMemcpyDeviceToHost);
  }
}

void REMDGPU::OutputCDV(const char *prefix,const unsigned int step,const int *index){
  for(int dev=0;dev<ngpus;dev++){
    cudaSetDevice(dev+gpu_offset);

#ifdef SWITCHING
    cudaMemcpy(gcoor+dev*nmol*nreplica/ngpus, gcoor_dev[dev], nreplica*nmol*sizeof(int4)/ngpus,  cudaMemcpyDeviceToHost);
#endif
    cudaMemcpy(acoor+dev*natom*nreplica/ngpus, acoor_dev[dev], nreplica*natom*sizeof(int4)/ngpus,  cudaMemcpyDeviceToHost);
    cudaMemcpy(L+dev*nreplica/ngpus,L_dev[dev],nreplica*sizeof(float3)/ngpus,cudaMemcpyDeviceToHost);
  }

  const float ImaxInv = 2.3283064365387e-10f;
  for(int rep=0;rep<nreplica;rep++){
    //int idx = index[rep];
    int idx = rep;
    const float3 ItoF = {L[idx].x*ImaxInv,L[idx].y*ImaxInv,L[idx].z*ImaxInv};

    const double T = cond0[rep].x * unit_temp;
    const double P = cond0[rep].y * unit_press;

    std::stringstream strs;
    strs << prefix;
    strs << "R" << std::setw(4) << std::setfill('0') << rep;
    //strs << "P" << (int)P << "T" << (int)T;
    strs << "s" << std::setw(7) << std::setfill('0') << step;
    strs << ".cdv";

    std::ostream *s;
    s = new std::ofstream(strs.str().c_str(),std::ofstream::out);
    if(s->fail()){
      std::cerr << "error: open " << strs.str() << " failed" << std::endl;
      return;
    }
    *s << "'";
    const int cmode = ((mode>>CSHIFT)&MASK);
    if(cmode==0){
      *s << " box_sx=" << -0.5*L[idx].x; //start of box for x axis
      *s << " box_sy=" << -0.5*L[idx].y; //start of box for y axis
      *s << " box_sz=" << -0.5*L[idx].z; //start of box for z axis
      *s << " box_ex=" <<  0.5*L[idx].x; //end of box for x axis
      *s << " box_ey=" <<  0.5*L[idx].y; //end of box for y axis
      *s << " box_ez=" <<  0.5*L[idx].z; //end of box for z axis
    }else if(cmode == 1){
      *s << " box_sx=" << -wall_length;  //start of box for x axis
      *s << " box_sy=" << -wall_length;  //start of box for y axis
      *s << " box_sz=" << -0.5*L[idx].z; //start of box for z axis
      *s << " box_ex=" <<  wall_length;  //end of box for x axis
      *s << " box_ey=" <<  wall_length;  //end of box for y axis
      *s << " box_ez=" <<  0.5*L[idx].z; //end of box for z axis
    }else if(cmode == 2){
      *s << " box_sx=" << -0.5*L[idx].x;    //start of box for x axis
      *s << " box_sy=" << -0.5*L[idx].y;    //start of box for y axis
      *s << " box_sz=" << -0.5*wall_length; //start of box for z axis
      *s << " box_ex=" <<  0.5*L[idx].x;    //end of box for x axis
      *s << " box_ey=" <<  0.5*L[idx].y;    //end of box for y axis
      *s << " box_ez=" <<  0.5*wall_length; //end of box for z axis
    }

    // for water color
    *s << " r0=1.0 r1=0.8 r2=0.8";
    *s << " c0=(1.0,0.0,0.0)";
    *s << " c1=(1.0,1.0,1.0)";
    *s << " c2=(1.0,1.0,1.0)";
    *s << " bond0_c=(1.0,1.0,1.0)";
    *s << " st0=\"T=" << T << " K, P=" << P/1e6 << " MPa\"";
    *s << " st0_pos=(-2.0,-2.0)";
    *s << std::endl;

    // writing coordinates
#ifdef SWITCHING
    f3i1 *acoorf = (f3i1*)acoor;
    for(unsigned int i=0;i<nmol;i++){
      for(unsigned int j=0;j<nsite;j++){
	float3 ar = {(float)gcoor[nmol*idx+i].x*ItoF.x + acoorf[natom*idx + j*nmol + i].x,
		     (float)gcoor[nmol*idx+i].y*ItoF.y + acoorf[natom*idx + j*nmol + i].y,
		     (float)gcoor[nmol*idx+i].z*ItoF.z + acoorf[natom*idx + j*nmol + i].z};
	*s << i << " " << j;
	*s << " " << ar.x;
	*s << " " << ar.y;
	*s << " " << ar.z;
	*s << std::endl;
      }

#else
    for(unsigned int i=0;i<natom;i++){
      float3 ar = {(float)acoor[natom*idx+i].x*ItoF.x,
		   (float)acoor[natom*idx+i].y*ItoF.y,
		   (float)acoor[natom*idx+i].z*ItoF.z};
      *s << i << " " << i/nmol;
      *s << " " << ar.x;
      *s << " " << ar.y;
      *s << " " << ar.z;
      *s << std::endl;
#endif
    }
    // writing bonds for water
    /*
    for(unsigned int i=0;i<nmol;i++){
      *s << "CDVIEW_BOND " << i << " " << i+nmol   << " 0" << std::endl;
      *s << "CDVIEW_BOND " << i << " " << i+2*nmol << " 0" << std::endl;
    }
    //*/
    delete s;
  }
};

void REMDGPU::ConvertVariablesToHost(Molecules **mlcls){
  GetBackAllVariablesFromGPU();
  const float ImaxInv = 2.3283064365387e-10f;
  for(int rep=0;rep<nreplica;rep++){
    for(int i=0;i<nmol;i++){
      mlcls[rep]->mlcl[i].r[0] = (double)(ImaxInv*(float)gcoor[nmol*rep + i].x) + 0.5;
      mlcls[rep]->mlcl[i].r[1] = (double)(ImaxInv*(float)gcoor[nmol*rep + i].y) + 0.5;
      mlcls[rep]->mlcl[i].r[2] = (double)(ImaxInv*(float)gcoor[nmol*rep + i].z) + 0.5;

      mlcls[rep]->mlcl[i].v[0] = (double)vel[nmol*rep + i].x;
      mlcls[rep]->mlcl[i].v[1] = (double)vel[nmol*rep + i].y;
      mlcls[rep]->mlcl[i].v[2] = (double)vel[nmol*rep + i].z;

      mlcls[rep]->mlcl[i].f[0] = (double)force[nmol*rep + i].x;
      mlcls[rep]->mlcl[i].f[1] = (double)force[nmol*rep + i].y;
      mlcls[rep]->mlcl[i].f[2] = (double)force[nmol*rep + i].z;

      mlcls[rep]->mlcl[i].q[0] = (double)angle[nmol*rep + i].x;
      mlcls[rep]->mlcl[i].q[1] = (double)angle[nmol*rep + i].y;
      mlcls[rep]->mlcl[i].q[2] = (double)angle[nmol*rep + i].z;
      mlcls[rep]->mlcl[i].q[3] = (double)angle[nmol*rep + i].w;

      mlcls[rep]->mlcl[i].p[0] = (double)angvel[nmol*rep + i].x;
      mlcls[rep]->mlcl[i].p[1] = (double)angvel[nmol*rep + i].y;
      mlcls[rep]->mlcl[i].p[2] = (double)angvel[nmol*rep + i].z;
      mlcls[rep]->mlcl[i].p[3] = (double)angvel[nmol*rep + i].w;

      mlcls[rep]->mlcl[i].n[0] = (double)torque[nmol*rep + i].x;
      mlcls[rep]->mlcl[i].n[1] = (double)torque[nmol*rep + i].y;
      mlcls[rep]->mlcl[i].n[2] = (double)torque[nmol*rep + i].z;
      mlcls[rep]->mlcl[i].n[3] = (double)torque[nmol*rep + i].w;
    }
    mlcls[rep]->tst->s  = (double)tst[rep].x;
    mlcls[rep]->tst->Ps = (double)tst[rep].y;

    mlcls[rep]->bst->Pv[0] = (double)bst[rep].x;
    mlcls[rep]->bst->Pv[1] = (double)bst[rep].y;
    mlcls[rep]->bst->Pv[2] = (double)bst[rep].z;

    mlcls[rep]->L[0] = (double)L[rep].x;
    mlcls[rep]->L[1] = (double)L[rep].y;
    mlcls[rep]->L[2] = (double)L[rep].z;

    mlcls[rep]->prop.vir[0] = (double)potvir[rep].x;
    mlcls[rep]->prop.vir[1] = (double)potvir[rep].y;
    mlcls[rep]->prop.vir[2] = (double)potvir[rep].z;
    mlcls[rep]->prop.pot    = (double)potvir[rep].w;
  }
 }

void REMDGPU::OutputCheckPoint(const char *prefix,const int *index){
  GetBackAllVariablesFromGPU();

  for(unsigned int rep=0;rep<nreplica;rep++){
    std::stringstream strs;
    strs << prefix << std::setw(4) << std::setfill('0') << rep << ".chk";
    
    std::ostream *os;
    os = new std::ofstream(strs.str().c_str(),std::ofstream::out);
    if(os->fail()){
      std::cerr << "error: open " << strs.str() << " failed" << std::endl;
      return;
    }

    const float ImaxInv = 2.3283064365387e-10f;
    int idx = index[rep];

    *os << nmol << " " << natom << std::endl;
    for(unsigned int i=0;i<nmol;i++){
      float  s = tst[rep].x;
      float3 r = {((float)gcoor[i+nmol*idx].x*ImaxInv+0.5f)*L[idx].x,
		  ((float)gcoor[i+nmol*idx].y*ImaxInv+0.5f)*L[idx].y,
		  ((float)gcoor[i+nmol*idx].z*ImaxInv+0.5f)*L[idx].z};
      float3 v = {vel[i+nmol*idx].x/(L[idx].x*s),vel[i+nmol*idx].y/(L[idx].y*s),vel[i+nmol*idx].z/(L[idx].z*s)};
      qtype4 q = { angle[i+nmol*idx].x, angle[i+nmol*idx].y, angle[i+nmol*idx].z, angle[i+nmol*idx].w};
      ptype4 p = {angvel[i+nmol*idx].x/s,angvel[i+nmol*idx].y/s,angvel[i+nmol*idx].z/s,angvel[i+nmol*idx].w/s};

      float  m = vel[i+nmol*idx].w;
      int  type= gcoor[i+nmol*idx].w;
      float3 in= {iner[type].x,iner[type].y,iner[type].z};
      int    id= i * nsite;// you have to fix here for multi molecule type

      *os << std::scientific;
      *os << " " << std::setw(12) << r.x;
      *os << " " << std::setw(12) << r.y;
      *os << " " << std::setw(12) << r.z;

      *os << " " << std::setw(12) << v.x;
      *os << " " << std::setw(12) << v.y;
      *os << " " << std::setw(12) << v.z;

      *os << " " << std::setw(12) << q.x;
      *os << " " << std::setw(12) << q.y;
      *os << " " << std::setw(12) << q.z;
      *os << " " << std::setw(12) << q.w;

      *os << " " << std::setw(12) << p.x;
      *os << " " << std::setw(12) << p.y;
      *os << " " << std::setw(12) << p.z;
      *os << " " << std::setw(12) << p.w;

      *os << " " << std::setw(12) << m;

      *os << " " << std::setw(12) << in.x;
      *os << " " << std::setw(12) << in.y;
      *os << " " << std::setw(12) << in.z;

      *os << " " << std::setw(4)  << type;
      *os << " " << std::setw(4)  << id;
      *os << std::endl;
    }
    float3 t = tst[idx];
    *os << " " << std::setw(8) << t.x;
    *os << " " << std::setw(8) << t.y;
    *os << " " << std::setw(8) << 0.5f/t.z;
    *os << std::endl;
    float4 b = bst[idx];
    *os << " " << std::setw(8) << b.x;
    *os << " " << std::setw(8) << b.y;
    *os << " " << std::setw(8) << b.z;
    *os << " " << std::setw(8) << 0.5f/b.w;
    *os << std::endl;

    *os << " " << std::setw(8) << L[idx].x;
    *os << " " << std::setw(8) << L[idx].y;
    *os << " " << std::setw(8) << L[idx].z;
    *os << std::endl;

    delete os;
  }
}

/*
 WRITE(12) NSTEP,
          (
	  X(I),Y(I),Z(I),          // coordinate of gravity center
	  X1(I),Y1(I),Z1(I),
	  X2(I),Y2(I),Z2(I),
	  X3(I),Y3(I),Z3(I),
	  X4(I),Y4(I),Z4(I),
	  X5(I),Y5(I),Z5(I),
          A(I),B(I),C(I),D(I),     // angle of molecule
	  A1(I),B1(I),C1(I),D1(I),
	  A2(I),B2(I),C2(I),D2(I),
	  A3(I),B3(I),C3(I),D3(I),
	  A4(I),B4(I),C4(I),D4(I),
	  VA(I),VB(I),VC(I),
	  VA1(I),VB1(I),VC1(I),
	  VA2(I),VB2(I),VC2(I),
	  VA3(I),VB3(I),VC3(I),
	  VA4(I),VB4(I),VC4(I),
          XT(I),YT(I),ZT(I),
	  XO(I),YO(I),ZO(I),       // coordinate of gravity center
	  I=1,N),
          VL0,                     // length of nanotube
	  VL1,VL2,VL3,VL4,VL5,
	  WMC1,
	  SL0,SL1,SL2,SL3,SL4,SL5,
	  QC1,RF,
	  PLZ1                     // chiral vector of zigzag carbon nanotube
*/
/*
void REMDGPU::OutputTakaiwaContinue(const char *prefix,const int *index){
  GetBackAllVariablesFromGPU();

  for(unsigned int rep=0;rep<nreplica;rep++){
    std::stringstream strs;
    strs << prefix << std::setw(4) << std::setfill('0') << rep << ".dat";

    std::ostream *os;
    os = new std::ofstream(strs.str(),std::ios_base::out);
    if(os->fail()){
      std::cerr << "error: open " << strs.str() << " failed" << std::endl;
      return;
    }

    const float ImaxInv = 2.3283064365387e-10f;
    int idx = index[rep];

    *os << nmol << " " << natom << std::endl;
    for(unsigned int i=0;i<nmol;i++){
      float  s = tst[rep].x;
      float3 r = {((float)gcoor[i+nmol*idx].x*ImaxInv+0.5f)*L[idx].x,
		  ((float)gcoor[i+nmol*idx].y*ImaxInv+0.5f)*L[idx].y,
		  ((float)gcoor[i+nmol*idx].z*ImaxInv+0.5f)*L[idx].z};
      qtype4 q = { angle[i+nmol*idx].x,
		   angle[i+nmol*idx].y,
		   angle[i+nmol*idx].z,
		   angle[i+nmol*idx].w};

      *os << std::scientific;
      *os << " " << std::setw(12) << r.z;     //X(I),Y(I),Z(I)
      *os << " " << std::setw(12) << r.y;
      *os << " " << std::setw(12) << r.x;

      *os << " " << "0.0d0 0.0d0 0.0d0";      //X1(I),Y1(I),Z1(I)
      *os << " " << "0.0d0 0.0d0 0.0d0";      //X2(I),Y2(I),Z2(I)
      *os << " " << "0.0d0 0.0d0 0.0d0";      //X3(I),Y3(I),Z3(I)
      *os << " " << "0.0d0 0.0d0 0.0d0";      //X4(I),Y4(I),Z4(I)
      *os << " " << "0.0d0 0.0d0 0.0d0";      //X5(I),Y5(I),Z5(I)

      *os << " " << std::setw(12) << q.x;     //A(I),B(I),C(I),D(I)
      *os << " " << std::setw(12) << q.y;
      *os << " " << std::setw(12) << q.z;
      *os << " " << std::setw(12) << q.w;

      *os << " " << "0.0d0 0.0d0 0.0d0 0.0d0";//A1(I),B1(I),C1(I),D1(I)
      *os << " " << "0.0d0 0.0d0 0.0d0 0.0d0";//A2(I),B2(I),C2(I),D2(I)
      *os << " " << "0.0d0 0.0d0 0.0d0 0.0d0";//A3(I),B3(I),C3(I),D3(I)
      *os << " " << "0.0d0 0.0d0 0.0d0 0.0d0";//A4(I),B4(I),C4(I),D4(I)

      *os << " " << "0.0d0 0.0d0 0.0d0";      //VA(I),VB(I),VC(I)
      *os << " " << "0.0d0 0.0d0 0.0d0";      //VA1(I),VB1(I),VC1(I),
      *os << " " << "0.0d0 0.0d0 0.0d0";      //VA2(I),VB2(I),VC2(I),
      *os << " " << "0.0d0 0.0d0 0.0d0";      //VA3(I),VB3(I),VC3(I),
      *os << " " << "0.0d0 0.0d0 0.0d0";      //VA4(I),VB4(I),VC4(I),

      *os << " " << "0.0d0 0.0d0 0.0d0";      //XT(I),YT(I),ZT(I),

      *os << " " << std::setw(12) << r.z;     //XO(I),YO(I),ZO(I),
      *os << " " << std::setw(12) << r.y;
      *os << " " << std::setw(12) << r.x;
    }
    *os << " " << std::setw(8) << L[idx].z;   //VL0

    *os << " " << "0.0d0 0.0d0 0.0d0 0.0d0 0.0d0";       //VL1,VL2,VL3,VL4,VL5,
    *os << " " << "0.0d0";                               //WMC1,
    *os << " " << "0.0d0 0.0d0 0.0d0 0.0d0 0.0d0 0.0d0"; //SL0,SL1,SL2,SL3,SL4,SL5,
    *os << " " << "0.0d0 0.0d0";                         //QC1,RF,
    *os << " " << std::setw(8) << 16.0;                  //PLZ1

    delete os;
  }
}
//*/

void REMDGPU::PrintAll(std::ostream &s){
  const float ImaxInv = 2.3283064365387e-10f;
  for(unsigned int rep=0;rep<nreplica;rep++){
    const float3 ItoF = {L[rep].x*ImaxInv,L[rep].y*ImaxInv,L[rep].z*ImaxInv};
    s << "#Coordinate" << std::endl;
    for(unsigned int mol=0;mol<nmol;mol++){
      s << mol << " " << (float)gcoor[mol].x*ItoF.x+0.5*L[rep].x << " " << (float)gcoor[mol].y*ItoF.y+0.5*L[rep].y << " " << (float)gcoor[mol].z*ItoF.z+0.5*L[rep].z << std::endl;
    }
    s << "#Velocity (and mass)" << std::endl;
    for(unsigned int mol=0;mol<nmol;mol++){
      s << mol << " " << vel[mol].x << " " << vel[mol].y << " " << vel[mol].z << " " << vel[mol].w << std::endl;
    }
    s << "#Force" << std::endl;
    for(unsigned int mol=0;mol<nmol;mol++){
      s << mol << " " << force[mol].x << " " << force[mol].y << " " << force[mol].z << std::endl;
    }
    s << "#Angle" << std::endl;
    for(unsigned int mol=0;mol<nmol;mol++){
      qtype sq = angle[mol].x*angle[mol].x + angle[mol].y*angle[mol].y + angle[mol].z*angle[mol].z + angle[mol].w*angle[mol].w;
      s << mol << " " << angle[mol].x << " " << angle[mol].y << " " << angle[mol].z << " " << angle[mol].w << " " << sq << std::endl;
    }
    s << "#Angular Velocity" << std::endl;
    for(unsigned int mol=0;mol<nmol;mol++){
      s << mol << " " << angvel[mol].x << " " << angvel[mol].y << " " << angvel[mol].z << " " << angvel[mol].w << std::endl;
    }
    s << "#Torque" << std::endl;
    for(unsigned int mol=0;mol<nmol;mol++){
      s << mol << " " << torque[mol].x << " " << torque[mol].y << " " << torque[mol].z << " " << torque[mol].w << std::endl;
    }

    s << "#Coordinate of atoms" << std::endl;
    for(unsigned int mol=0;mol<natom;mol++){
      s << mol << " " << (float)acoor[mol].x*ItoF.x+0.5*L[rep].x << " " << (float)acoor[mol].y*ItoF.y+0.5*L[rep].y << " " << (float)acoor[mol].z*ItoF.z+0.5*L[rep].z << std::endl;
    }
    s << "#Force of atoms" << std::endl;
    for(unsigned int mol=0;mol<natom;mol++){
      s << mol << " " << aforce[mol].x << " " << aforce[mol].y << " " << aforce[mol].z << " " << aforce[mol].w << std::endl;
    }
    s << "#Thermostat" << std::endl;
    s << tst[rep].x << " " << tst[rep].y << " " << tst[rep].z << std::endl;
    s << "#Barostat" << std::endl;
    s << bst[rep].x << " " << bst[rep].y << " " << bst[rep].z << " " << bst[rep].w << std::endl;
    s << "#Cell length" << std::endl;
    s << L[rep].x << " " << L[rep].y << " " << L[rep].z << std::endl;
  }
}

void REMDGPU::Diff(Molecules **mlcls,std::ostream &s){
  for(int rep=0;rep<nreplica;rep++){
    Molecule *m = mlcls[rep]->mlcl;
    s << "#Coordinate" << std::endl;
    for(int i=0;i<nmol;i++){
      float3 r = {(float)gcoor[i+nmol*rep].x/(float)IntMax+0.5f,
		  (float)gcoor[i+nmol*rep].y/(float)IntMax+0.5f,
		  (float)gcoor[i+nmol*rep].z/(float)IntMax+0.5f};
      s << i << " " << (r.x - m[i].r[0])/m[i].r[0] << " " << (r.y - m[i].r[1])/m[i].r[1] << " " << (r.z - m[i].r[2])/m[i].r[2] << std::endl;
    }
    s << "#Velocity" << std::endl;
    for(int i=0;i<nmol;i++){
      float4 v = vel[i+nmol*rep];
      s << i << " " << (v.x - m[i].v[0])/m[i].v[0] << " " << (v.y - m[i].v[1])/m[i].v[1] << " " << (v.z - m[i].v[2])/m[i].v[2] << std::endl;
    }
    s << "#Force" << std::endl;
    for(int i=0;i<nmol;i++){
      float4 f = force[i+nmol*rep];
      s << i << " " << (f.x - m[i].f[0])/m[i].f[0] << " " << (f.y - m[i].f[1])/m[i].f[1] << " " << (f.z - m[i].f[2])/m[i].f[2] << std::endl;
    }
    s << "#Angle" << std::endl;
    for(int i=0;i<nmol;i++){
      qtype4 q = angle[i+nmol*rep];
      s << i << " " << (q.x - m[i].q[0])/m[i].q[0] << " " << (q.y - m[i].q[1])/m[i].q[1] << " " << (q.z - m[i].q[2])/m[i].q[2] << " " << (q.w - m[i].q[3])/m[i].q[3] << std::endl;
    }
    s << "#Angvel" << std::endl;
    for(int i=0;i<nmol;i++){
      ptype4 p = angvel[i+nmol*rep];
      s << i << " " << (p.x - m[i].p[0])/m[i].p[0] << " " << (p.y - m[i].p[1])/m[i].p[1] << " " << (p.z - m[i].p[2])/m[i].p[2] << " " << (p.w - m[i].p[3])/m[i].p[3] << std::endl;
    }
    s << "#Torque" << std::endl;
    for(int i=0;i<nmol;i++){
      float4 n = torque[i+nmol*rep];
      s << i << " " << (n.x - m[i].n[0])/m[i].n[0] << " " << (n.y - m[i].n[1])/m[i].n[1] << " " << (n.z - m[i].n[2])/m[i].n[2] << " " << (n.w - m[i].n[3])/m[i].n[3] << std::endl;
    }
    Atom *a = mlcls[rep]->atom;
    s << "#Coordinate of atoms" << std::endl;
    for(int i=0;i<natom;i++){
      float3 r = {((float)acoor[i+natom*rep].x/(float)IntMax+0.5f)*L[rep].x,
		  ((float)acoor[i+natom*rep].y/(float)IntMax+0.5f)*L[rep].y,
		  ((float)acoor[i+natom*rep].z/(float)IntMax+0.5f)*L[rep].z};
      int cpuind = (nsite*i)%natom+i/nmol;
      s << i << " " << (r.x - a[cpuind].r[0])/a[cpuind].r[0] << " " << (r.y - a[cpuind].r[1])/a[cpuind].r[1] << " " << (r.z - a[cpuind].r[2])/a[cpuind].r[2] << std::endl;
    }
    s << "#Force of atoms" << std::endl; // This never match because gpu force includes intra force
    for(int i=0;i<natom;i++){
      float4 f = aforce[i+nmol*rep];
      int cpuind = (nsite*i)%natom+i/nmol;
      s << i << " " << (f.x - a[cpuind].f[0])/a[cpuind].f[0] << " " << (f.y - a[cpuind].f[1])/a[cpuind].f[1] << " " << (f.z - a[cpuind].f[2])/a[cpuind].f[2] << std::endl;
    }
  }
}

#ifdef FASTMATH
  #define _EXP __expf
  #define _SIN __sinf
  #define _COS __cosf
  #define _SINCOS __sincosf
#else
  #define _EXP exp
  #define _SIN sin
  #define _COS cos
  #define _SINCOS sincosf
#endif

template <int MODE>
__device__
float4 CalcForcePotVirGPU
(int4               *coor,
 float4             *force,
 float2             *kvecr,
 int3               *kvc,
 const float3       L,
 const float        rcutSq,
 const float        alpha,
 const unsigned int nwave,
 const unsigned int natom,
 const unsigned int nsite,
 const float R
){
  const float ImaxInv = 2.3283064365387e-10f;
  const float3 ItoF = {L.x*ImaxInv,L.y*ImaxInv,L.z*ImaxInv};
  float4 potvir = {0.f,0.f,0.f,0.f};

  const float alpha2 = alpha*alpha;
  const float alphai = 1.f / alpha;
  const float coef = 2.f / sqrt(M_PI);

  int nmol = natom / nsite;
  __shared__ int4 coor_shared[NUMSHARED];
#ifdef HIGHACC
  VGGSI ptmp,vtmpx,vtmpy,vtmpz;
  INIT_NARUMI(ptmp);
  INIT_NARUMI(vtmpx);
  INIT_NARUMI(vtmpy);
  INIT_NARUMI(vtmpz);
  float ys,tmp1,tmp2;
#endif
  // short range force
  for(unsigned int i=0;i<(natom+bdm-1)/bdm;i++){
    int4 coor_i;
    if(i*bdm+tid<natom){
      coor_i = coor[i*bdm+tid];
    }else{
      coor_i.w = NMOLTYPE*nsite;
    }
    const float aq_i = alpha*chg[coor_i.w];
#ifdef HIGHACC
    VGGSI ftmpx,ftmpy,ftmpz,ftmpw;
    INIT_NARUMI(ftmpx);
    INIT_NARUMI(ftmpy);
    INIT_NARUMI(ftmpz);
    INIT_NARUMI(ftmpw);
#endif
    float4 force_i = {0.f,0.f,0.f,0.f};
    if(((MODE>>CSHIFT)&MASK)==1){
      if(i*bdm+tid<nmol){
      float2 r_wall = {coor_i.x*ItoF.x,coor_i.y*ItoF.y};
      float2 table  = tex1D(tex_wall,sqrtf(r_wall.x*r_wall.x+r_wall.y*r_wall.y)/R);
#ifdef HIGHACC
      ys = table.x*r_wall.x;
      ADD_NARUMI(ftmpx,ys,tmp1,tmp2);
      ys = table.x*r_wall.y;
      ADD_NARUMI(ftmpy,ys,tmp1,tmp2);
      ys = table.y;
      ADD_NARUMI(ptmp,ys,tmp1,tmp2);
#else
      force_i.x += table.x*r_wall.x;
      force_i.y += table.x*r_wall.y;
      potvir.w += table.y;
      }
#endif
    }
    if(((MODE>>CSHIFT)&MASK)==2){
      if(i*bdm+tid<nmol){
	float  r_wall = coor_i.z*ItoF.z;
	float2 table  = tex1D(tex_wall,fabs(r_wall/R));
	force_i.z += table.x*r_wall/fabs(r_wall);
	potvir.w  += table.y;
      }
    }
    for(unsigned int s=0;s<(natom+NUMSHARED-1)/NUMSHARED;s++){
      __syncthreads();
      if(tid<NUMSHARED){
	coor_shared[tid] = coor[s*NUMSHARED + tid];
      }
      __syncthreads();
      for(unsigned int j=0;j<NUMSHARED;++j){
	const int4 coor_j = coor_shared[j];
	//if(coor_i.w == NMOLTYPE*nsite) continue;
	if((s*NUMSHARED+j)%nmol==(i*bdm+tid)%nmol || (s*NUMSHARED+j)>=natom) continue;
	const float x = (float)(coor_i.x - coor_j.x)*ItoF.x;
	const float y = (float)(coor_i.y - coor_j.y)*ItoF.y;
	const float z = (float)(coor_i.z - coor_j.z)*ItoF.z;

	const float r02 = x*x + y*y + z*z;
	if(r02 > rcutSq) continue;
	const float r01   = sqrt(r02);
	const float r01i  = 1.f / r01;
	const float r02i  = r01i*r01i;
	//Ewald part
	const float ar  = alpha * r01;
	const float ari = alphai * r01i;
	const float aqq = aq_i * chg[coor_j.w];
	const float eari= erfc(ar)*ari;
	float f = alpha2*aqq*(eari + coef*exp(-ar*ar))*ari*ari;
	float e = 0.5*aqq*eari; // halven for summation
	//LJ part
	const float rs02i = sgm02[coor_i.w][coor_j.w]*r02i;
	const float rs06i = rs02i*rs02i*rs02i;

	f += 6.0f * eps[coor_i.w][coor_j.w] * rs06i * (2.f*rs06i - 1.f) * r02i;
	e += 0.5f * eps[coor_i.w][coor_j.w] * rs06i * (rs06i - 1.f); // halven for summation
#ifdef HIGHACC
	ys = f * x;
	ADD_NARUMI(ftmpx,ys,tmp1,tmp2);
	ys = f * y;
	ADD_NARUMI(ftmpy,ys,tmp1,tmp2);
	ys = f * z;
	ADD_NARUMI(ftmpz,ys,tmp1,tmp2);
	ys = e;
	ADD_NARUMI(ptmp, ys,tmp1,tmp2);
	ys = f * x * x;
	ADD_NARUMI(vtmpx,ys,tmp1,tmp2);
	ys = f * y * y;
	ADD_NARUMI(vtmpy,ys,tmp1,tmp2);
	ys = f * z * z;
	ADD_NARUMI(vtmpz,ys,tmp1,tmp2);

	if(j%16==0){
	  UPDATE_NARUMI(ftmpx);
	  UPDATE_NARUMI(ftmpy);
	  UPDATE_NARUMI(ftmpz);
	  UPDATE_NARUMI(ftmpw);

	  UPDATE_NARUMI(vtmpx);
	  UPDATE_NARUMI(vtmpy);
	  UPDATE_NARUMI(vtmpz);
	  UPDATE_NARUMI(ptmp);
	}
#else
	force_i.x += f*x;
	force_i.y += f*y;
	force_i.z += f*z;
	//force_i.w += e;

	potvir.w += e;
	potvir.x += f*x*x;
	potvir.y += f*y*y;
	potvir.z += f*z*z;
#endif
      }
    }
#ifdef HIGHACC
    force_i.x = COPY_NARUMI(ftmpx);
    force_i.y = COPY_NARUMI(ftmpy);
    force_i.z = COPY_NARUMI(ftmpz);
    //force_i.w = COPY_NARUMI(ftmpw);
#endif
    force[i*bdm+tid] = force_i;
    __syncthreads();
  }
  potvir.x *= 0.5f;
  potvir.y *= 0.5f;
  potvir.z *= 0.5f;

#if 1
  // long range force
  const float pa02 = (M_PI/alpha)*(M_PI/alpha);
  const float3 Li = {1.f/L.x, 1.f/L.y, 1.f/L.z};
  const float coef1 = 2.f * Li.x * Li.y * Li.z;
  const float coef2 = 1.f / (2.f * M_PI * L.x * L.y * L.z);

  float mr;

  //wave part
  for (int i = 0; i < (nwave+bdm-1)/bdm; i++){
    __syncthreads();
    float3 vec;
    if(i*bdm+tid<nwave){
      vec.x = (float)kvc[i*bdm+tid].x * Li.x;
      vec.y = (float)kvc[i*bdm+tid].y * Li.y;
      vec.z = (float)kvc[i*bdm+tid].z * Li.z;
    }else{
      vec.x = vec.y = vec.z = 0.f;
    }
#ifdef HIGHACC
    VGGSI recvecx,recvecy;
    INIT_NARUMI(recvecx);
    INIT_NARUMI(recvecy);
#else
    float2 recvec = {0.f, 0.f};
#endif
    //forward dft
    for(unsigned int s=0;s<(natom+NUMSHARED-1)/NUMSHARED;s++){
      __syncthreads();
      if((s*NUMSHARED+tid) < natom){
	coor_shared[tid] = coor[s*NUMSHARED + tid];
      }else{
	coor_shared[tid].x = 0;
	coor_shared[tid].y = 0;
	coor_shared[tid].z = 0;
	coor_shared[tid].w = NMOLTYPE*nsite;
      }
      __syncthreads();
      for (unsigned int j = 0; j < NUMSHARED; j++){
	if((s*NUMSHARED+j)>=natom) continue;
	//const int4 r = coor_shared[j];
	const int4 r = coor[s*NUMSHARED+j];
	mr  = (float)r.x * ItoF.x * vec.x;
	mr += (float)r.y * ItoF.y * vec.y;
	mr += (float)r.z * ItoF.z * vec.z;
	mr *= 6.2831853f; // 2*M_PI (M_PI = 3.14159265358973)
#ifdef HIGHACC
	ys = chg[r.w] * _SIN(mr);
	ADD_NARUMI(recvecx,ys,tmp1,tmp2);
	ys = chg[r.w] * _COS(mr);
	ADD_NARUMI(recvecy,ys,tmp1,tmp2);
	if(j%16==0){
	UPDATE_NARUMI(recvecx);
	UPDATE_NARUMI(recvecy);
	}
#else
	recvec.x += chg[r.w] * _SIN(mr);
	recvec.y += chg[r.w] * _COS(mr);
	//float sinmr,cosmr;
	//_SINCOS(mr,&sinmr,&cosmr);
	//recvec.x += chg[r.w] * sinmr;
	//recvec.y += chg[r.w] * cosmr;
#endif
      }
    }
    if(vec.x==0 && vec.y==0 && vec.z==0) continue; // you need continue here to avoid synchrads() in branch
    const float3 vecvec = {vec.x*vec.x, vec.y*vec.y, vec.z*vec.z};
    const float  vec02  = vecvec.x + vecvec.y + vecvec.z;
    const float  vec02i = 1.f/vec02;
    const float  coef3  = _EXP(-pa02*vec02) * vec02i;

#ifdef HIGHACC
    float2 recvec;
    recvec.x = COPY_NARUMI(recvecx);
    recvec.y = COPY_NARUMI(recvecy);
#endif

    float vtmp = recvec.x*recvec.x + recvec.y*recvec.y;
    recvec.x *= coef3;
    recvec.y *= coef3;
    kvecr[i*bdm+tid] = recvec;
    const float  coef4  = coef2 * coef3 * vtmp;
    vtmp = 2.f * (vec02i + pa02);

#ifdef HIGHACC
    ys = coef4*(1.f - vtmp * vec.x * vec.x);
    ADD_NARUMI(vtmpx,ys,tmp1,tmp2);
    ys = coef4*(1.f - vtmp * vec.y * vec.y);
    ADD_NARUMI(vtmpy,ys,tmp1,tmp2);
    ys = coef4*(1.f - vtmp * vec.z * vec.z);
    ADD_NARUMI(vtmpz,ys,tmp1,tmp2);
    if(i%16==0){
      UPDATE_NARUMI(vtmpx);
      UPDATE_NARUMI(vtmpy);
      UPDATE_NARUMI(vtmpz);
    }
#else
    potvir.x += coef4*(1.f - vtmp * vec.x * vec.x);
    potvir.y += coef4*(1.f - vtmp * vec.y * vec.y);
    potvir.z += coef4*(1.f - vtmp * vec.z * vec.z);
#endif
  }
  __syncthreads();

  //backward dft
  __shared__ float2 kvecr_shared[NUMSHARED];
  for(int j=0;j<(natom+bdm-1)/bdm;j++){
    int4   ir;
    if(j*bdm+tid<natom){
      ir = coor[j*bdm+tid];
    }else{
      ir.x = 0;
      ir.y = 0;
      ir.z = 0;
      ir.w = NMOLTYPE*nsite;
    }
    float3 r = {(float)ir.x*ItoF.x, (float)ir.y*ItoF.y, (float)ir.z*ItoF.z};
#ifdef HIGHACC
    VGGSI ftmpx,ftmpy,ftmpz,ftmpw;
    float ys,tmp1,tmp2;
    INIT_NARUMI(ftmpx);
    INIT_NARUMI(ftmpy);
    INIT_NARUMI(ftmpz);
    INIT_NARUMI(ftmpw);
#endif
    float4 ftmp = {0.f, 0.f, 0.f, 0.f};
    for(unsigned int s=0;s<(nwave+NUMSHARED-1)/NUMSHARED;s++){
      __syncthreads();
      if((s*NUMSHARED+tid)<nwave){
	kvecr_shared[tid] = kvecr[s*NUMSHARED + tid];
      }else{
	kvecr_shared[tid].x = 0.f;
	kvecr_shared[tid].y = 0.f;
      }
      __syncthreads();
      for (int i = 0; i < NUMSHARED; i++){
	float3 vec;
	float2 vecr;
	//const float2 vecr = kvecr_shared[i];
	if((NUMSHARED*s+i) < nwave){
	  vec.x = (float)kvc[NUMSHARED*s+i].x*Li.x;
	  vec.y = (float)kvc[NUMSHARED*s+i].y*Li.y;
	  vec.z = (float)kvc[NUMSHARED*s+i].z*Li.z;
	  vecr  = kvecr_shared[i];
        }else{
	  vec.x = 0.f;
	  vec.y = 0.f;
	  vec.z = 0.f;
	  vecr.x = 0.f;
	  vecr.y = 0.f;
        }
	mr  = r.x * vec.x;
	mr += r.y * vec.y;
	mr += r.z * vec.z;
	mr *= 6.2831853f;
	const float sinmr = _SIN(mr);
	const float cosmr = _COS(mr);
	//float sinmr,cosmr;
	//_SINCOS(mr,&sinmr,&cosmr);
	const float f = chg[ir.w] * coef1 * (vecr.y * sinmr - vecr.x * cosmr);
#ifdef HIGHACC
	ys = f * vec.x;
	ADD_NARUMI(ftmpx,ys,tmp1,tmp2);
	ys = f * vec.y;
	ADD_NARUMI(ftmpy,ys,tmp1,tmp2);
	ys = f * vec.z;
	ADD_NARUMI(ftmpz,ys,tmp1,tmp2);
	ys = chg[ir.w] * coef2 * (vecr.y * cosmr + vecr.x * sinmr);
	ADD_NARUMI(ftmpw,ys,tmp1,tmp2);
	ADD_NARUMI(ptmp,ys,tmp1,tmp2);
	if(i%16==0){
	UPDATE_NARUMI(ftmpx);
	UPDATE_NARUMI(ftmpy);
	UPDATE_NARUMI(ftmpz);
	UPDATE_NARUMI(ftmpw);
	UPDATE_NARUMI(ptmp);
        }
#else
	ftmp.x += f * vec.x;
	ftmp.y += f * vec.y;
	ftmp.z += f * vec.z;
	//ftmp.w += chg[ir.w] * coef2 * (vecr.y * cosmr + vecr.x * sinmr);
	potvir.w += chg[ir.w] * coef2 * (vecr.y * cosmr + vecr.x * sinmr);
#endif
      }
    }
#ifdef HIGHACC
    ftmp.x = COPY_NARUMI(ftmpx);
    ftmp.y = COPY_NARUMI(ftmpy);
    ftmp.z = COPY_NARUMI(ftmpz);
    ftmp.w = COPY_NARUMI(ftmpw);
#endif
    force[j*bdm+tid].x += ftmp.x;
    force[j*bdm+tid].y += ftmp.y;
    force[j*bdm+tid].z += ftmp.z;
    force[j*bdm+tid].w += ftmp.w;
    __syncthreads();
  }
  __syncthreads();

#ifdef HIGHACC
  potvir.x = COPY_NARUMI(vtmpx);
  potvir.y = COPY_NARUMI(vtmpx);
  potvir.z = COPY_NARUMI(vtmpx);
  potvir.w = COPY_NARUMI(ptmp);
#endif
#endif // long range force
  potvir = block_all_reduce(potvir);
  return potvir;

#undef EXP
#undef SIN
#undef COS
}


template <int MODE>
__device__
float4 CalcForcePotVirSwitchingGPU2
(int4         *gcoor,
 float4       *gforce,
 f3i1         *acoor,
 float4       *aforce,
 Result       &rslt,
 const float3 L,
 const float  rcutSq,
 const float  rswitch,
 const float  R,
 const int nmol,
 const int nsite
){
  const float ImaxInv = 2.3283064365387e-10f;
  const float3 ItoF = {L.x*ImaxInv,L.y*ImaxInv,L.z*ImaxInv};
  __shared__ int4 gcoor_shared[NTHREADS];
  __shared__ f3i1 acoor_shared[MAXNATOM*NTHREADS];

  // switching function variables
  const float rc  = sqrtf(rcutSq);
  const float rl  = rswitch;
  const float rlc = rl - rc;
  const float rcli = 1.f/(rc - rl);
  const float coef = 1.f / (rlc*rlc*rlc*rlc*rlc);
  float sw=1.f,dsw=0.f;

  float3 vir = make_float3(0.f);
  float lj_sum=0.f,clmb_sum=0.f,wall_sum=0.f;

  for(int i=0;i<(nmol+bdm-1)/bdm;i++){
    int ind_i = i*bdm + tid;
    int4 gcoor_i;
    f3i1 ai[MAXNATOM];
    float4 fi[MAXNATOM];

    if(ind_i < nmol){
      gcoor_i = gcoor[ind_i];
      gforce[ind_i] = make_float4(0.f);
      for(int ii=0;ii<nsite;ii++){
	ai[ii] = acoor[ind_i + ii*nmol];
	fi[ii] = make_float4(0.0);
      }
    }

    // register variables ?
    float x,y,z,r01,r02,r01i,r02i,rs06i;
    float4 force_i_sum = make_float4(0.f);

    for(int s=0;s<(nmol+NUMSHARED-1)/NUMSHARED;s++){
      // store variables to shared memoery
      __syncthreads();
      if(tid<NUMSHARED && tid+s*NUMSHARED<nmol){
	gcoor_shared[tid] = gcoor[tid+s*NUMSHARED];
	for(int k=0;k<nsite;k++){
	  acoor_shared[tid + k*NUMSHARED] = acoor[tid + s*NUMSHARED + k*nmol];
	}
      }
      __syncthreads();

      if(ind_i < nmol)
#pragma unroll
      for(int j=0;j<NUMSHARED;j++){
	if(j+s*NUMSHARED >= nmol) break;     // skip j >= nmol
	if(ind_i == j+s*NUMSHARED) continue; // skip i == j
	float4 force_i = make_float4(0.f);

	//cut-off using distance between gravity center
	//const int4 gcoor_j = gcoor[j];
	const int4 gcoor_j = gcoor_shared[j];
	const float gz = (float)(gcoor_i.z - gcoor_j.z)*ItoF.z;
	const float gx = (float)(gcoor_i.x - gcoor_j.x)*ItoF.x;
	const float gy = (float)(gcoor_i.y - gcoor_j.y)*ItoF.y;
	r02 = gx*gx + gy*gy + gz*gz;
	if(r02 > rcutSq) continue; // skip molecule outside the cut-off radious
	// switching function
	r01 = sqrtf(r02);
#if 0
	if(r01 > rl){
	  const float r01c = r01 - rc;
	  const float r01l = r01 - rl;
	  sw  = coef*r01c*r01c*r01c*(10.f*r01l*r01l - 5.f*r01l*r01c + r01c*r01c);
	  dsw = coef*30.f*r01c*r01c*r01l*r01l;
	}else{
	  sw  = 1.f;
	  dsw = 0.f;
	}
#else
	const float2 switching = tex1D(tex_switch,(r01-rl)*rcli);
	sw  = switching.x;
	dsw = switching.y;
#endif
	r01i = 1.f/r01;
	const float dswx = dsw*gx*r01i;
	const float dswy = dsw*gy*r01i;
	const float dswz = dsw*gz*r01i;

	float lje=0.f,clmbe=0.f,etmp,e,ftmp;
	// LJ part
	f3i1 aj = acoor_shared[j];
	x = gx + ai[0].x - aj.x;
	y = gy + ai[0].y - aj.y;
	z = gz + ai[0].z - aj.z;
	r02   = x*x + y*y + z*z;
	r02i  = 1.f/r02;
	rs06i = sgm02[ai[0].w][aj.w]*r02i;
	rs06i = rs06i*rs06i*rs06i;

	ftmp = 12.f * eps[ai[0].w][aj.w] * rs06i * (rs06i - 0.5f) * r02i * sw;
	fi[0].x += ftmp*x;
	fi[0].y += ftmp*y;
	fi[0].z += ftmp*z;
	force_i.x += ftmp*x;
	force_i.y += ftmp*y;
	force_i.z += ftmp*z;

	etmp = eps[ai[0].w][aj.w] * rs06i * (rs06i - 1.f);
	lje += etmp;
	lj_sum += etmp * sw;
	//Coulomb part
#pragma unroll
	for(int ii=1;ii<nsite;ii++){
#pragma unroll
	  for(int jj=1;jj<nsite;jj++){
	    aj = acoor_shared[j+jj*NUMSHARED];
	    const float qq = chg[ai[ii].w] * chg[aj.w];
	    x = gx + ai[ii].x - aj.x;
	    y = gy + ai[ii].y - aj.y;
	    z = gz + ai[ii].z - aj.z;
	    r02   = x*x + y*y + z*z;
	    r02i  = 1.f/r02;
	    r01i = sqrtf(r02i);

	    ftmp = qq*r01i*r02i*sw;
	    fi[ii].x += ftmp*x;
	    fi[ii].y += ftmp*y;
	    fi[ii].z += ftmp*z;
	    force_i.x += ftmp*x;
	    force_i.y += ftmp*y;
	    force_i.z += ftmp*z;

	    etmp = qq*r01i;
	    clmbe += etmp;
	    clmb_sum += etmp * sw;
	  }// jj roop
	}// ii roop
	e = lje + clmbe;
	force_i.x -= e * dswx;
	force_i.y -= e * dswy;
	force_i.z -= e * dswz;
	force_i_sum += force_i;
	vir.x += force_i.x * gx;
	vir.y += force_i.y * gy;
	vir.z += force_i.z * gz;
	//if(tid==0) printf("%d %f %f %f\n",tid,force_i.x,force_i.y,force_i.z);
      }// j roop
    }// s roop
    // confined force and potential
    if(ind_i < nmol)
    if(((MODE>>CSHIFT)&MASK)==1){
      x = gcoor_i.x*ItoF.x + ai[0].x;
      y = gcoor_i.y*ItoF.y + ai[0].y;

      r02 = x*x + y*y;
      r01 = sqrtf(r02);

      float2 table  = tex1D(tex_wall,r01/R);
      fi[0].x       += table.x * x;
      force_i_sum.x += table.x * x;
      fi[0].y       += table.x * y;
      force_i_sum.y += table.x * y;

      //potvir.w  += table.y;
      wall_sum  += table.y;
    }
    if(((MODE>>CSHIFT)&MASK)==2){
      if(ind_i<nmol){
	const float  z   = gcoor_i.z*ItoF.z + ai[0].z;
	const float  za  = fabs(z);
	const float  dir = (za!=0.f) ? z / za : 1.f;
	const float2 table  = tex1D(tex_wall,2.0*za/R);

	fi[0].z       += table.x * dir;
	force_i_sum.z += table.x * dir;
	wall_sum      += table.y;
	//printf("%e %e %e %e\n",za,table.x,table.y,dir);
      }
    }
    if(ind_i<nmol){
      //printf("%d %f %f %f\n",ind_i,force_i_sum.x,force_i_sum.y,force_i_sum.z);
      gforce[ind_i] = force_i_sum;
      for(int ii=0;ii<nsite;ii++){
	aforce[ind_i+ii*nmol] = fi[ii];
      }
    }
  }// i roop
  __syncthreads();

  lj_sum   = 0.5f * block_all_reduce(lj_sum);
  clmb_sum = 0.5f * block_all_reduce(clmb_sum);
  wall_sum = block_all_reduce(wall_sum);
  if(tid==0){
    rslt.lj   = lj_sum;
    rslt.clmb = clmb_sum;
    rslt.wall = wall_sum;
  }

  vir   = 0.5f * block_all_reduce(vir);
  float4 potvir;
  potvir.x = vir.x;
  potvir.y = vir.y;
  potvir.z = vir.z;
  potvir.w = lj_sum + clmb_sum + wall_sum;

  return potvir;
}

template <int MODE>
__device__
float4 CalcForcePotVirDirectGPU
(int4               *coor,
 float4             *force,
 float2             *kvecr,
 int3               *kvc,
 const float3       L,
 const float        rcutSq,
 const float        alpha,
 const unsigned int nwave,
 const unsigned int natom,
 const unsigned int nsite,
 const float R
){
  const float ImaxInv = 2.3283064365387e-10f;
  const float3 ItoF = {L.x*ImaxInv,L.y*ImaxInv,L.z*ImaxInv};
  int nmol = natom / nsite;
  float4 potvir = {0.f,0.f,0.f,0.f};
  float4 potvir2= {0.f,0.f,0.f,0.f};
  __shared__ int4 coor_shared[NUMSHARED];
  int n = 0;
  do{
  potvir2 = potvir;
  potvir = make_float4(0.f);
  for(unsigned int i=0;i<(natom+bdm-1)/bdm;i++){
    int4 coor_i;
    if(i*bdm+tid<natom){
      coor_i = coor[i*bdm+tid];
    }else{
      coor_i.w = NMOLTYPE*nsite;
    }
    const float q_i = chg[coor_i.w];
    float4 force_i = {0.f,0.f,0.f,0.f};
    if(((MODE>>CSHIFT)&MASK)==1){
      if(i*bdm+tid<nmol){
      float2 r_wall = {coor_i.x*ItoF.x,coor_i.y*ItoF.y};
      float2 table  = tex1D(tex_wall,sqrtf(r_wall.x*r_wall.x+r_wall.y*r_wall.y)/R);
      force_i.x += table.x*r_wall.x;
      force_i.y += table.x*r_wall.y;
      potvir.w += table.y;
      }
    }

    const float xi = (float)coor_i.x*ItoF.x;
    const float yi = (float)coor_i.y*ItoF.y;
    const float zi = (float)coor_i.z*ItoF.z;
    for(unsigned int s=0;s<(natom+NUMSHARED-1)/NUMSHARED;s++){
      __syncthreads();
      if(tid<NUMSHARED){
	coor_shared[tid] = coor[s*NUMSHARED + tid];
      }
      __syncthreads();
      for(unsigned int j=0;j<NUMSHARED;++j){
	const int4 coor_j = coor_shared[j];
	if((s*NUMSHARED+j)>=natom) continue;
	for(int nz=-n;nz<=n;nz++){
	  //if(((MODE>>CSHIFT)&MASK)==2 && nz!=0) continue;
	  const float z = zi - (float)coor_j.z*ItoF.z + nz*L.z;
	  for(int ny=-n;ny<=n;ny++){
	    //if(((MODE>>CSHIFT)&MASK)==1 && ny!=0) continue;
	    const float y = yi - (float)coor_j.y*ItoF.y + ny*L.y;
	    for(int nx=-n;nx<=n;nx++){
	      //if(((MODE>>CSHIFT)&MASK)==1 && nx!=0) continue;
	      const float x = xi - (float)coor_j.x*ItoF.x + nx*L.x;
	      //if(j==0 && tid==0) printf("%d %d %d: %d %d %d\n",nx,ny,nz,abs((int)(x/L.x))-abs(nx),abs((int)(y/L.y))-abs(ny),abs((int)(z/L.z))-abs(nz));
	      if(nx==0 && ny==0 && nz==0 && (s*NUMSHARED+j)%nmol==(i*bdm+tid)%nmol)continue;

	      const float r02  = x*x + y*y + z*z;
	      const float r01  = sqrt(r02);
	      //if(tid==0 && j==0) printf("%d %d %d: %lf\n",nx,ny,nz,r01);
	      const float r01i = 1.f / r01;
	      const float r02i = r01i*r01i;
	      //coulomb part
	      const float qq = q_i * chg[coor_j.w];
	      float f = qq*r02i;
	      float e = qq*r01i;
	      //LJ part
	      if(r02 <= rcutSq){
		const float rs02i = sgm02[coor_i.w][coor_j.w]*r02i;
		const float rs06i = rs02i*rs02i*rs02i;
		f += 6.f * eps[coor_i.w][coor_j.w] * rs06i * (2.f*rs06i - 1.f) * r02i;
		e += eps[coor_i.w][coor_j.w] * rs06i * (rs06i - 1.f);
	      }
	      if(nx==0 && ny==0 && nz==0) e *= 0.5f; // halven for newton 3rd law

	      force_i.x += f*x;
	      force_i.y += f*y;
	      force_i.z += f*z;
	      //force_i.w += e;

	      potvir.w += e;
	      potvir.x += f*x*x;
	      potvir.y += f*y*y;
	      potvir.z += f*z*z;
	    }// nx roop
	  }// ny roop
	}// nz roop
      }// j roop
    }// s roop
    force[i*bdm+tid] = force_i;
    __syncthreads();
  }// i roop
  potvir = block_all_reduce(potvir);
  if(tid==0){
    //printf("%d %f\n",n,potvir.w);
    //printf("iteration %d: pot = %f\n",n,potvir.w);
    //printf("force[0]: %f %f %f\n",force[0].x,force[0].y,force[0].z);
  }
  n++;
  __syncthreads();
  }while(n<=15);

  return potvir;
}
/*
template<int MODE>
__global__
void calc_direct_potvir
(
 int4   *r_dev,               // coordinate of gravity center (x,y,z) and molecule type (w)
 float4 *v_dev,               // translational velocity (x,y,z) and mass (w)
 float4 *f_dev,               // force (x,y,z) and 
 qtype4 *q_dev,               // angle
 ptype4 *p_dev,               // angular velocity
 float4 *n_dev,               // torque
 int4   *ar_dev,              // coordinate of atoms
 float4 *af_dev,              // force of atoms
 float2 *kr_dev,              // coordinate in wave space
 int3   *kvc_dev,             // wave vector
 float4 *potvir_dev,
 float3 *tst_dev,             // thermostat property
 float4 *bst_dev,             // barostat property
 float3 *L_dev,               // cell edge length (x,y,z)
 float2 *cond0_dev,           // temperature and pressure of each replica
 float2 *cond1_dev,           // temperature and pressure of each replica
 float4 *rslt_dev,            // results (potential, volume, Ptmp, Ttmp)
 float  *H0_dev,              // Hamiltonian of system
 const unsigned int nmol,     // # of molecules
 //const unsigned int natom,    // # of atoms
 const float dt,              // delta timer
 const float rc02,            // square of cut off length
 const float alpha,           // alpha for ewald sum
 const unsigned int  nwave,   // # of wave vector
 const unsigned long interval, // # of steps for a interval
 const float wall_length
 ){
  int4   *r_rep = r_dev + nmol*bid;
  float4 *v_rep = v_dev + nmol*bid;
  float4 *f_rep = f_dev + nmol*bid;
  qtype4 *q_rep = q_dev + nmol*bid;
  ptype4 *p_rep = p_dev + nmol*bid;
  float4 *n_rep = n_dev + nmol*bid;

  int4   *ar_rep = ar_dev + nsite*nmol*bid;
  float4 *af_rep = af_dev + nsite*nmol*bid;

  float2 *kr_rep  = kr_dev  + nwave*bid;
  int3   *kvc_rep = kvc_dev + nwave*bid;

  float sigma_c = 1.423; // sigma of carbon

  float3 L      = L_dev[bid];
  double V;
  if(((MODE>>CSHIFT)&MASK) == 0) V = L.x*L.y*L.z;
  if(((MODE>>CSHIFT)&MASK) == 1) V = 3.1415926f*(wall_length-0.5f*sigma_c)*(wall_length-0.5f*sigma_c)*L.z;
  if(((MODE>>CSHIFT)&MASK) == 2) V = L.x*L.y*wall_length;

  float  dthalf = 0.5f*dt;
  float3 t      = tst_dev[bid]; // .x = s, .y = Ps, .z = 1.f/Q
  float4 b      = bst_dev[bid]; // .x = Pv[0], .y = Pv[1], .x = Pv[2]

  float4 potvir = potvir_dev[bid];
  float  H0  = H0_dev[bid];
  float  gkT = 6.f*(float)nmol*cond0_dev[bid].x;
  float  P   = cond0_dev[bid].y;

    // convert molecules to atom
    for(unsigned int k=0;k<nsite;k++){
      for(unsigned int j=tid;j<nmol;j+=bdm){
	const int4   r = r_rep[j];
	const qtype4 q = q_rep[j];
	const int  ind = nsite*r.w + k;
	const qtype4 sq = {q.x*q.x, q.y*q.y, q.z*q.z, q.w*q.w};
	const qtype  xy = 2.f*q.x*q.y, xz = 2.f*q.x*q.z, xw = 2.f*q.x*q.w;
	const qtype  yz = 2.f*q.y*q.z, yw = 2.f*q.y*q.w, zw = 2.f*q.z*q.w;
	qtype4 ar;
	ar.x  = (sq.x + sq.y - sq.z - sq.w)*(qtype)atm[ind].x;
	ar.x += (yz - xw)*(qtype)atm[ind].y;
	ar.x += (yw + xz)*(qtype)atm[ind].z;
	ar.x /= (qtype)L.x;

	ar.y  = (yz + xw)*(qtype)atm[ind].x;
	ar.y += (sq.x - sq.y + sq.z - sq.w)*(qtype)atm[ind].y;
	ar.y += (zw - xy)*(qtype)atm[ind].z;
	ar.y /= (qtype)L.y;

	ar.z  = (yw - xz)*(qtype)atm[ind].x;
	ar.z += (zw + xy)*(qtype)atm[ind].y;
	ar.z += (sq.x - sq.y - sq.z + sq.w)*(qtype)atm[ind].z;
	ar.z /= (qtype)L.z;

	ar_rep[k*nmol+j].x = (int)(ar.x*(qtype)IntMax)+r.x;
	ar_rep[k*nmol+j].y = (int)(ar.y*(qtype)IntMax)+r.y;
	ar_rep[k*nmol+j].z = (int)(ar.z*(qtype)IntMax)+r.z;
	ar_rep[k*nmol+j].w = ind;
      }
    }
    for(unsigned int k=0;k<nsite;k++){
      for(unsigned int j=tid;j<nmol;j+=bdm){
	af_rep[k*nmol+j] = make_float4(0.f);
	//av_rep[k*nsite+j] = zero;
      }
    }
    __syncthreads();
    // force calculation
    potvir =   CalcForcePotVirDirectGPU<MODE>
      (ar_rep,af_rep,kr_rep,kvc_rep,
       L,rc02,alpha,nwave,nsite*nmol,wall_length);
    __syncthreads();

    //convert atoms to molecules
    float3 vir = make_float3(0.f);
    for(unsigned int j=tid;j<nmol;j+=bdm){
      float4 f = {0.f,0.f,0.f,0.f};
      float4 n = {0.f,0.f,0.f,0.f};
      int4   r_i = r_rep[j];
      for(unsigned int k=0;k<nsite;k++){
	const unsigned int ind = nsite*r_rep[j].w + k;
	const int4   ar = ar_rep[nmol*k+j];
	const float4 af = af_rep[nmol*k+j];
	const qtype4 q = q_rep[j];
	f.x += af.x;
	f.y += af.y;
	f.z += af.z;
	f.w += af.w;
	const qtype4 sq = {q.x*q.x, q.y*q.y, q.z*q.z, q.w*q.w};
	const qtype xy = 2.f*q.x*q.y, xz = 2.f*q.x*q.z, xw = 2.f*q.x*q.w;
	const qtype yz = 2.f*q.y*q.z, yw = 2.f*q.y*q.w, zw = 2.f*q.z*q.w;
	float3 fb; // convert space force to body force
	fb.x  = (float)(sq.x + sq.y - sq.z - sq.w) * af.x;
	fb.x += (float)(yz + xw) * af.y;
	fb.x += (float)(yw - xz) * af.z;

	fb.y  = (float)(yz - xw) * af.x;
	fb.y += (float)(sq.x - sq.y + sq.z - sq.w) * af.y;
	fb.y += (float)(zw + xy) * af.z;

	fb.z  = (float)(yw + xz) * af.x;
	fb.z += (float)(zw - xy) * af.y;
	fb.z += (float)(sq.x - sq.y - sq.z + sq.w) * af.z;

	n.x += atm[ind].x * fb.x + atm[ind].y * fb.y + atm[ind].z * fb.z;
	n.y += atm[ind].y * fb.z - atm[ind].z * fb.y;
	n.z += atm[ind].z * fb.x - atm[ind].x * fb.z;
	n.w += atm[ind].x * fb.y - atm[ind].y * fb.x;

	// difference between atomic pressure and molecular pressure
	float3 r_ai = make_float3((float)(ar.x - r_i.x)/(float)IntMax*L.x,
				  (float)(ar.y - r_i.y)/(float)IntMax*L.y,
				  (float)(ar.z - r_i.z)/(float)IntMax*L.z);
	vir.x += r_ai.x * af.x;
	vir.y += r_ai.y * af.y;
	vir.z += r_ai.z * af.z;
      }
      f_rep[j] = f;
      n_rep[j] = n;
    }

    if(tid==0) rslt_dev[bid].x = potvir.w;
}
//*/

__inline__ __device__
float4 CalcKin(int4 *r_rep,float4 *v_rep,qtype4 *q_rep,ptype4 *p_rep,float3 L,float3 t,unsigned int nmol){
  float4 sum = {0.f,0.f,0.f,0.f};
  float tmp;
  for(int i=tid;i<nmol;i+=bdm){
    const int4   r = r_rep[i];
    const float4 v = v_rep[i];
    const qtype4 q = q_rep[i];
    const ptype4 p = p_rep[i];

    sum.x += 0.5f * v.w * v.x*v.x/(L.x*L.x * t.x*t.x);
    sum.y += 0.5f * v.w * v.y*v.y/(L.y*L.y * t.x*t.x);
    sum.z += 0.5f * v.w * v.z*v.z/(L.z*L.z * t.x*t.x);

    tmp = (-p.x*q.y + p.y*q.x + p.z*q.w - p.w*q.z) * 0.25f / (inertia[r.w].x * t.x);
    sum.w += 2.f * tmp * tmp * inertia[r.w].x;

    tmp = (-p.x*q.z - p.y*q.w + p.z*q.x + p.w*q.y) * 0.25f / (inertia[r.w].y * t.x);
    sum.w += 2.f * tmp * tmp * inertia[r.w].y;

    tmp = (-p.x*q.w + p.y*q.z - p.z*q.y + p.w*q.x) * 0.25f / (inertia[r.w].z * t.x);
    sum.w += 2.f * tmp * tmp * inertia[r.w].z;
  }
  return block_all_reduce(sum);
}

__inline__ __device__
float3 CalcMom(float4 *v_rep,float3 L,float3 t,unsigned int nmol){
  float3 sum = {0.f,0.f,0.f};
  for(int i=tid;i<nmol;i+=bdm){
    const float4 v = v_rep[i];
    sum.x += v.w * v.x / (L.x*t.x);
    sum.y += v.w * v.y / (L.y*t.x);
    sum.z += v.w * v.z / (L.z*t.x);
  }
  return sum;
}

// 2.0f * dthalf = dt
 __inline__ __device__
void D6(float3 &t,const float dthalf){
  const float tmp = 1.f + (t.y * t.z * dthalf);
  t.x *= tmp * tmp;
  t.y /= tmp;
}

__inline__ __device__
void D5(float4 &v, const float4 f, ptype4 &p, const qtype4 q, const float4 n, const float3 t, const float3 L, const float dt, const float dthalf){
   v.x += f.x * L.x * t.x * dthalf / v.w;
   v.y += f.y * L.y * t.x * dthalf / v.w;
   v.z += f.z * L.z * t.x * dthalf / v.w;
   p.x += (q.x*(ptype)n.x - q.y*(ptype)n.y - q.z*(ptype)n.z - q.w*(ptype)n.w) * (ptype)t.x * (ptype)dt;
   p.y += (q.y*(ptype)n.x + q.x*(ptype)n.y - q.w*(ptype)n.z + q.z*(ptype)n.w) * (ptype)t.x * (ptype)dt;
   p.z += (q.z*(ptype)n.x + q.w*(ptype)n.y + q.x*(ptype)n.z - q.y*(ptype)n.w) * (ptype)t.x * (ptype)dt;
   p.w += (q.w*(ptype)n.x - q.z*(ptype)n.y + q.y*(ptype)n.z + q.x*(ptype)n.w) * (ptype)t.x * (ptype)dt;
 }

const float max_aspect_ratio = 4.0f / 1.0f;
 // Pv/W * dt/2 = Pv/2W * dt
template <const int MODE>
__inline__ __device__
void D4(float3 &t, const double4 b,
        float3 &L, const float wall_length,
        double &V, double &A,
        const float dt, const float dthalf)
{
  if(((MODE>>PSHIFT)&MASK) == 1){
    if(((MODE>>TSHIFT)&MASK) == 1)
      t.y -= b.x * b.x * b.w * dthalf;
    V += (double)(t.x*b.x*b.w*dt);
    L.x = L.y = L.z = (float)pow(V,1.0/3.0);
  }
  if(((MODE>>PSHIFT)&MASK) == 2){
    if(((MODE>>TSHIFT)&MASK) == 1)
      t.y -= b.z * b.z * b.w * dthalf;
    L.z += t.x * b.z * b.w * dt;
    V = A * L.z;
  }
  if(((MODE>>PSHIFT)&MASK) == 3){
    if(((MODE>>TSHIFT)&MASK) == 1)
      t.y -= (b.x*b.x + b.y*b.y) * b.w * dthalf;
    //if(max_aspect_ratio > L.x / L.y) L.x += t.x * b.x * b.w * dt;
    //if(max_aspect_ratio > L.y / L.x) L.y += t.x * b.y * b.w * dt;
    L.x += t.x * b.x * b.w * dt;
    L.y += t.x * b.y * b.w * dt;
    A = L.x * L.y;
    if(((MODE>>CSHIFT)&MASK) == 0) V = A * L.z;
    if(((MODE>>CSHIFT)&MASK) == 2) V = A * wall_length;
  }
}
/*
  float4 P3q,P3p;
  P3q.x = -q.w; P3p.x = -p.w;
  P3q.y =  q.z; P3p.y =  p.z;
  P3q.z = -q.y; P3p.z = -p.y;
  P3q.w =  q.x; P3p.w =  p.x;
//*/
//*
 template <const int MODE>
__inline__ __device__
void D3(const int4 r,ptype4 &p,qtype4 &q, float4 &kin,const float3 inertia[NMOLTYPE], const float3 t, const float dthalf){
  const qtype xi = (-p.x*q.w + p.y*q.z - p.z*q.y + p.w*q.x) / (qtype)(4.f*inertia[r.w].z * t.x);
  if(((MODE>>TSHIFT)&MASK) == 1){
    kin.w += 2.f * (float)xi * (float)xi * inertia[r.w].z * dthalf;
  }
  const qtype cosxidt = cos(xi*(qtype)dthalf);
  const qtype sinxidt = sin(xi*(qtype)dthalf);
  const qtype4 Pq = {-q.w, q.z,-q.y, q.x};
  const qtype4 Pp = {-p.w, p.z,-p.y, p.x};
  q.x = q.x * cosxidt + Pq.x * sinxidt;
  q.y = q.y * cosxidt + Pq.y * sinxidt;
  q.z = q.z * cosxidt + Pq.z * sinxidt;
  q.w = q.w * cosxidt + Pq.w * sinxidt;
  p.x = p.x * cosxidt + Pp.x * sinxidt;
  p.y = p.y * cosxidt + Pp.y * sinxidt;
  p.z = p.z * cosxidt + Pp.z * sinxidt;
  p.w = p.w * cosxidt + Pp.w * sinxidt;
}
//*/
template <const int MODE>
__inline__ __device__
void D2(const int4 r,ptype4 &p,qtype4 &q, float4 &kin,const float3 inertia[NMOLTYPE], const float3 t, const float dthalf){
  const qtype xi = (-p.x*q.z - p.y*q.w + p.z*q.x + p.w*q.y) / (qtype)(4.f * inertia[r.w].y * t.x);
  if(((MODE>>TSHIFT)&MASK) == 1){
    kin.w += 2.f * xi * xi * inertia[r.w].y * dthalf;
  }
  const qtype cosxidt = cos(xi*(qtype)dthalf);
  const qtype sinxidt = sin(xi*(qtype)dthalf);
  const qtype4 Pq = {-q.z,-q.w, q.x, q.y};
  const qtype4 Pp = {-p.z,-p.w, p.x, p.y};
  q.x = q.x * cosxidt + Pq.x * sinxidt;
  q.y = q.y * cosxidt + Pq.y * sinxidt;
  q.z = q.z * cosxidt + Pq.z * sinxidt;
  q.w = q.w * cosxidt + Pq.w * sinxidt;
  p.x = p.x * cosxidt + Pp.x * sinxidt;
  p.y = p.y * cosxidt + Pp.y * sinxidt;
  p.z = p.z * cosxidt + Pp.z * sinxidt;
  p.w = p.w * cosxidt + Pp.w * sinxidt;
}

/*
  Pq.x = -q.y;Pp.x = -p.y;
  Pq.y =  q.x;Pp.y =  p.x;
  Pq.z =  q.w;Pp.z =  p.w;
  Pq.w = -q.z;Pp.w = -p.z;
//*/
template <const int MODE>
__inline__ __device__
  void D1(int4 &r,float4 &v,
	  ptype4 &p,qtype4 &q,
	  float4 &kin,
	  const float3 inertia[NMOLTYPE], const float3 t, const float3 L,
	  const float dt)
{
  r.x += (int)( (float)IntMax * (v.x/(L.x*L.x*t.x)*dt) );
  r.y += (int)( (float)IntMax * (v.y/(L.y*L.y*t.x)*dt) );
  r.z += (int)( (float)IntMax * (v.z/(L.z*L.z*t.x)*dt) );
  const qtype xi = (-p.x*q.y + p.y*q.x + p.z*q.w - p.w*q.z) / (qtype)(4.f * inertia[r.w].x * t.x);
  const qtype cosxidt = cos(xi*(qtype)dt);
  const qtype sinxidt = sin(xi*(qtype)dt);
  const qtype4 Pq = {-q.y, q.x, q.w,-q.z};
  const qtype4 Pp = {-p.y, p.x, p.w,-p.z};
  q.x  = q.x * cosxidt + Pq.x * sinxidt;
  q.y  = q.y * cosxidt + Pq.y * sinxidt;
  q.z  = q.z * cosxidt + Pq.z * sinxidt;
  q.w  = q.w * cosxidt + Pq.w * sinxidt;
  p.x  = p.x * cosxidt + Pp.x * sinxidt;
  p.y  = p.y * cosxidt + Pp.y * sinxidt;
  p.z  = p.z * cosxidt + Pp.z * sinxidt;
  p.w  = p.w * cosxidt + Pp.w * sinxidt;
  if(((MODE>>TSHIFT)&MASK) == 1 || ((MODE>>PSHIFT)&MASK) > 0){
    kin.x += 0.5f*v.w * v.x*v.x/(L.x*L.x * t.x*t.x) * dt;
    kin.y += 0.5f*v.w * v.y*v.y/(L.y*L.y * t.x*t.x) * dt;
    kin.z += 0.5f*v.w * v.z*v.z/(L.z*L.z * t.x*t.x) * dt;
    kin.w += 2.f * xi * xi * inertia[r.w].x * dt;
  }
}
template<int MODE>
__global__
void vel_scale
(int4   *r_dev,
 float4 *v_dev,
 qtype4 *q_dev,
 ptype4 *p_dev,
 float3 *tst_dev,
 float3 *L_dev,
 float2 *cond0_dev,
 const int nmol)
{
  int4   *r_rep = r_dev + nmol*bid;
  float4 *v_rep = v_dev + nmol*bid;
  qtype4 *q_rep = q_dev + nmol*bid;
  ptype4 *p_rep = p_dev + nmol*bid;
  float3 L      = L_dev[bid];
  float3 t      = tst_dev[bid]; // .x = s, .y = Ps, .z = 1.f/Q

  float T = cond0_dev[bid].x;

  float4 kin = CalcKin(r_rep,v_rep,q_rep,p_rep,L,t,nmol);

  float vscale = sqrtf( 1.5f * (float)nmol * T / (kin.x+kin.y+kin.z) );
  float pscale = sqrtf( 1.5f * (float)nmol * T / kin.w );

  float4 mom = make_float4(0.f);

  for(int i=tid;i<nmol;i+=bdm){
    v_rep[i].x *= vscale;
    v_rep[i].y *= vscale;
    v_rep[i].z *= vscale;

    p_rep[i].x *= pscale;
    p_rep[i].y *= pscale;
    p_rep[i].z *= pscale;
    p_rep[i].w *= pscale;

    mom.x += v_rep[i].w * v_rep[i].x;
    mom.y += v_rep[i].w * v_rep[i].y;
    mom.z += v_rep[i].w * v_rep[i].z;
    mom.w += v_rep[i].w;
  }
  mom = block_all_reduce(mom);
  for(int i=tid;i<nmol;i+=bdm){
    if(((MODE>>CSHIFT)&MASK)!=1) v_rep[i].x -= mom.x/mom.w;
    if(((MODE>>CSHIFT)&MASK)!=1) v_rep[i].y -= mom.y/mom.w;
    if(((MODE>>CSHIFT)&MASK)!=2) v_rep[i].z -= mom.z/mom.w;
  }
  // need bug fix
  /*
  for(int i=tid;i<nmol;i+=bdm){
    ptype4 p = p_rep[i];
    const float scale = 1.f/sqrt(p.x*p.x + p.y*p.y + p.z*p.z + p.w*p.w);
    p.x *= scale;
    p.y *= scale;
    p.z *= scale;
    p.w *= scale;
    p_rep[i] = p;
  }
  //*/
}

__device__
void remove_angular_momentum_in_tube
(
 int4   *r_rep,
 float4 *v_rep,
 qtype4 *q_rep,
 ptype4 *p_rep,
 float3 t,
 float3 L,
 const int nmol
 )
{
  const float ImaxInv = 2.3283064365387e-10f;
  const float3 ItoF = {L.x*ImaxInv,L.y*ImaxInv,L.z*ImaxInv};

  float2 mom = make_float2(0.f);
  for(int i=tid;i<nmol;i+=bdm){
    float3 r = {r_rep[i].x*ItoF.x, r_rep[i].y*ItoF.y, r_rep[i].z*ItoF.z};
    float4 v = {v_rep[i].x/(L.x*t.x),v_rep[i].y/(L.y*t.x),v_rep[i].z/(L.z*t.x),v_rep[i].w};
    mom.x += v.w*(r.x*v.y - r.y*v.x);
    mom.y += v.w;
  }

  mom = block_all_reduce(mom);

  for(int i=tid;i<nmol;i+=bdm){
    const float3 r    = {r_rep[i].x*ItoF.x, r_rep[i].y*ItoF.y, r_rep[i].z*ItoF.z};
    const float  r02  = r.x*r.x + r.y*r.y;
    if(r02 > 1.f){ //if r = (0, 0, 0), v becomes infinity
      v_rep[i].x += mom.x * r.y / (mom.y * r02) * L.x * t.x;
      v_rep[i].y -= mom.x * r.x / (mom.y * r02) * L.y * t.x;
    }
  }
}


template<int MODE>
__global__
void remd_kernel
(
 int4   *r_dev,               // coordinate of gravity center (x,y,z) and molecule type (w)
 float4 *v_dev,               // translational velocity (x,y,z) and mass (w)
 float4 *f_dev,               // force (x,y,z) and 
 qtype4 *q_dev,               // angle
 ptype4 *p_dev,               // angular velocity
 float4 *n_dev,               // torque
 int4   *ar_dev,              // coordinate of atoms
 float4 *af_dev,              // force of atoms
 float2 *kr_dev,              // coordinate in wave space
 int3   *kvc_dev,             // wave vector
 float4 *potvir_dev,
 float3 *tst_dev,             // thermostat property
 float4 *bst_dev,             // barostat property
 float3 *L_dev,               // cell edge length (x,y,z)
 float2 *cond0_dev,           // temperature and pressure of each replica
 float2 *cond1_dev,           // temperature and pressure of each replica
 Result *rslt_dev,            // results (potential, volume, Ptmp, Ttmp)
 float  *H0_dev,              // Hamiltonian of system
 const unsigned int nmol,     // # of molecules
 //const unsigned int natom,  // # of atoms
 const unsigned int nsite,    // # of atoms in a water
 const float dt,              // delta time
 const float rc02,            // square of cut off length
 const float alpha,           // alpha for ewald sum
 const unsigned int  nwave,   // # of wave vector
 const float rswitch,         // distance where switching start
 const unsigned long interval,// # of steps for a interval
 const float wall_length
 ){
  int4   *r_rep = r_dev + nmol*bid;
  float4 *v_rep = v_dev + nmol*bid;
  float4 *f_rep = f_dev + nmol*bid;
  qtype4 *q_rep = q_dev + nmol*bid;
  ptype4 *p_rep = p_dev + nmol*bid;
  float4 *n_rep = n_dev + nmol*bid;

  int4   *ar_rep = ar_dev + nsite*nmol*bid;
  float4 *af_rep = af_dev + nsite*nmol*bid;
#ifndef SWITCHING
  float2 *kr_rep  = kr_dev  + nwave*bid;
  int3   *kvc_rep = kvc_dev + nwave*bid;
#endif
  float3 L      = L_dev[bid];
  double V,A;
  const float sigma_c = 3.4f; // sigma of carbon
  if(((MODE>>CSHIFT)&MASK) == 0){
    A = L.x*L.y;
    V = A * L.z;
  }
  if(((MODE>>CSHIFT)&MASK) == 1){
    A = 3.1415926f*(wall_length - 0.5f*sigma_c)*(wall_length - 0.5f*sigma_c);
    //A = 3.1415926f * wall_length * wall_length;
    V = A * L.z;
  }
  if(((MODE>>CSHIFT)&MASK) == 2){
    A = L.x * L.y;
    V = A * wall_length;
  }

  float  dthalf = 0.5f*dt;
  float3 t      = tst_dev[bid]; // .x = s, .y = Ps, .z = 1.f/Q
  //float4 b      = bst_dev[bid]; // .x = Pv[0], .y = Pv[1], .x = Pv[2]
  double4 b     = make_double4(bst_dev[bid]); // .x = Pv[0], .y = Pv[1], .x = Pv[2]

  float4 potvir = potvir_dev[bid];
  float  H0     = H0_dev[bid];
  float  gkT    = 6.f*(float)nmol*cond0_dev[bid].x;
  float  P      = cond0_dev[bid].y;

  float4 kin; // kinetic energy {x,y,z} = {tra.x , tra.y, tra.z, rot}
  float3 mom; // momentum of the system

  //float P_LRC = 16.f/3.f*3.1415926f * (float)(nmol*nmol) * 0.650194 * powf(3.16557f,3) * (2.f * powf(3.16557f/sqrtf(rc02),9) - 3.f * powf(3.16557f/sqrtf(rc02),3));

#ifdef HIGHACC
  VGGSI tmpqp;
  float ys,tmp1,tmp2;
#endif
  // velocity scaling after temperature exchange
  if(cond0_dev[bid].x != cond1_dev[bid].x){
    if(((MODE>>TSHIFT)&MASK) == 1){
      //scaling velocities
      const float scale = sqrtf(cond0_dev[bid].x/cond1_dev[bid].x);
      float mass = 0.f;
      mom = make_float3(0.f);
      for(unsigned int j=tid;j<nmol;j+=bdm){
	float4 v = v_rep[j];
	v.x  *= scale/t.x;
	v.y  *= scale/t.x;
	v.z  *= scale/t.x;
	mom  += v*v.w;
	mass += v.w;
	v_rep[j] = v;

	ptype4 p = p_rep[j];
	p.x *= scale/t.x;
	p.y *= scale/t.x;
	p.z *= scale/t.x;
	p.w *= scale/t.x;
	p_rep[j] = p;
      }
      t.y *= scale;
      if(((MODE>>PSHIFT)&MASK) == 1) b.x *= scale*t.x;
      if(((MODE>>PSHIFT)&MASK) == 2) b.z *= scale*t.x;
      if(((MODE>>PSHIFT)&MASK) == 3){b.x *= scale*t.x; b.y *= scale*t.x;}
      //reset s
      //without reset, s becomes smaller/larger for low/high temperature, and
      //s sometimes becomes smaller/larger and smaller/larger with replica exchanges
      t.x = 1.f;
      //remove momentum
      mom  = block_all_reduce(mom);
      mass = block_all_reduce(mass);
      mom /= mass;
      for(unsigned int j=tid;j<nmol;j+=bdm){
	float4 v = v_rep[j];
	if(((MODE>>CSHIFT)&MASK) == 0) v   -= mom;
	if(((MODE>>CSHIFT)&MASK) == 1) v.z -= mom.z;
	if(((MODE>>CSHIFT)&MASK) == 2){v.x -= mom.x; v.y -= mom.y;}
	v_rep[j] = v;
      }
    }
    if(((MODE>>CSHIFT)&MASK) == 1){
      remove_angular_momentum_in_tube(r_rep,v_rep,q_rep,p_rep,t,L,nmol);
    }

    //calc new hamiltonian
    kin = CalcKin(r_rep,v_rep,q_rep,p_rep,L,t,nmol);
    H0  = potvir.w + sum(kin);
    if(((MODE>>TSHIFT)&MASK) == 1){
      H0 += gkT*log(t.x) + t.y*t.y*t.z;
    }

    if(((MODE>>PSHIFT)&MASK) == 1)
      H0 += (b.x*b.x + b.y*b.y + b.z*b.z)*b.w;
    if(((MODE>>PSHIFT)&MASK) == 2)
      H0 += b.z*b.z*b.w;
    if(((MODE>>PSHIFT)&MASK) == 3)
      H0 += (b.x*b.x + b.y*b.y)*b.w;
    if(((MODE>>PSHIFT)&MASK) > 0) H0 += P*V;
  }

  if(cond0_dev[bid].y != cond1_dev[bid].y){
    H0  = potvir.w + rslt_dev[bid].kin;
    if(((MODE>>TSHIFT)&MASK) == 1){
      H0 += gkT*log(t.x) + t.y*t.y*t.z;
    }
    if(((MODE>>PSHIFT)&MASK) == 1)
      H0 += (b.x*b.x + b.y*b.y + b.z*b.z)*b.w;
    if(((MODE>>PSHIFT)&MASK) == 2)
      H0 += b.z*b.z*b.w;
    if(((MODE>>PSHIFT)&MASK) == 3)
      H0 += (b.x*b.x + b.y*b.y)*b.w;
    if(((MODE>>PSHIFT)&MASK) > 0) H0 += P*V;
  }
  // main roop
  for(unsigned long s=0;s<interval;s++){
#ifdef STEEPEST
    const float dr   = dt;
    const qtype dphi = 2.f / t.z;
    const float dv   = 2.f / b.w / P;
    for(unsigned int i=tid;i<nmol;i+=bdm){
      r_rep[i].x += (int)( (float)IntMax * f_rep[i].x * dr);
      r_rep[i].y += (int)( (float)IntMax * f_rep[i].y * dr);
      r_rep[i].z += (int)( (float)IntMax * f_rep[i].z * dr);
      qtype4 q = q_rep[i];
      float4 n = n_rep[i];
      q_rep[i].x += 2.*(q.x*n.x - q.y*n.y - q.z*n.z - q.w*n.w)* dphi;
      q_rep[i].y += 2.*(q.y*n.x + q.x*n.y - q.w*n.z + q.z*n.w)* dphi;
      q_rep[i].z += 2.*(q.z*n.x + q.w*n.y + q.x*n.z - q.y*n.w)* dphi;
      q_rep[i].w += 2.*(q.w*n.x - q.z*n.y + q.y*n.z + q.x*n.w)* dphi;
      q_rep[i]   /= norm(q_rep[i]);
    }

    if(((MODE>>PSHIFT)&MASK) == 1){
      V += (double)((potvir.x+potvir.y+potvir.z)/(3.f*V) - P)*dv;
      L.x = L.y = L.z = (float)pow(V,1.0/3.0);
    }
    if(((MODE>>PSHIFT)&MASK) == 2){
      L.z += (potvir.z / L.z - P * A) * dv;
      V = A * L.z;
    }
    if(((MODE>>PSHIFT)&MASK) == 3){
      if(((MODE>>CSHIFT)&MASK) == 2){
	L.x += (potvir.x/L.x - P * L.y*wall_length)*dv;
	L.y += (potvir.y/L.y - P * L.x*wall_length)*dv;
      }else{
	L.x += (potvir.x/L.x - P * L.y*L.z)*dv;
	L.y += (potvir.y/L.y - P * L.x*L.z)*dv;
      }
      A = L.x * L.y;
      if(((MODE>>CSHIFT)&MASK) == 0) V = A * L.z;
      if(((MODE>>CSHIFT)&MASK) == 2) V = A * wall_length;
    }
    //*/
#else //STEEPEST

    if(((MODE>>TSHIFT)&MASK) == 1){
      D6(t,dthalf);
      t.y -= potvir.w*dthalf; // D5 integration
      if(((MODE>>PSHIFT)&MASK) > 0){
	t.y -= P * V * dthalf; // D5 integration
      }
    }
    if(((MODE>>PSHIFT)&MASK) == 1){
      b.x += t.x*((potvir.x+potvir.y+potvir.z)/(3.f*V) - P)*dthalf;
    }
    if(((MODE>>PSHIFT)&MASK) == 2){
      b.z += t.x * potvir.z / L.z * dthalf;
      b.z -= t.x * P * A * dthalf;
    }
    if(((MODE>>PSHIFT)&MASK) == 3){
      if(((MODE>>CSHIFT)&MASK) == 2){
	b.x += t.x*(potvir.x/L.x - P * L.y*wall_length)*dthalf;
	b.y += t.x*(potvir.y/L.y - P * L.x*wall_length)*dthalf;
      }else{
	b.x += t.x*(potvir.x/L.x - P * L.y*L.z)*dthalf;
	b.y += t.x*(potvir.y/L.y - P * L.x*L.z)*dthalf;
      }
    }
    if(((MODE>>TSHIFT)&MASK) == 1 || ((MODE>>PSHIFT)&MASK)>0){
      kin.x = kin.y = kin.z = kin.w = 0.f;
    }
    for(unsigned int j=tid;j<nmol;j+=bdm){
      // load to register
      int4   r = r_rep[j];
      float4 v = v_rep[j];
      qtype4 q = q_rep[j];
      ptype4 p = p_rep[j];
      const float4 f = f_rep[j];
      const float4 n = n_rep[j];
      D5(v,f,p,q,n,t,L,dt,dthalf);
      if(j==tid){ D4<MODE>(t,b,L,wall_length,V,A,dt,dthalf);}
      D3<MODE>(r,p,q,kin,inertia,t,dthalf);
      D2<MODE>(r,p,q,kin,inertia,t,dthalf);
      D1<MODE>(r,v,p,q,kin,inertia,t,L,dt);
      D2<MODE>(r,p,q,kin,inertia,t,dthalf);
      D3<MODE>(r,p,q,kin,inertia,t,dthalf);

      // store to grobal memory
      r_rep[j] = r;
      v_rep[j] = v;
      q_rep[j] = q;
      p_rep[j] = p;
    }
    //D1 integration for t and b
    if(((MODE>>TSHIFT)&MASK) == 1 || ((MODE>>PSHIFT)&MASK) > 0){
      kin = block_all_reduce(kin);
    }
    if(((MODE>>TSHIFT)&MASK) == 1){
      t.y += kin.x + kin.y + kin.z + kin.w - ( gkT*(1.f+log(t.x)) - H0 )*dt; // D1, D2, D3, D5 integration for t.y
    }
    //D1 integration for P constant (dt is included in kin)
    if(((MODE>>PSHIFT)&MASK) == 1){
      b.x += t.x * 2.f * (kin.x+kin.y+kin.z) / (3.f*V);
    }
    if(((MODE>>PSHIFT)&MASK) == 2){
      b.z += t.x * 2.f * kin.z / L.z;
    }
    if(((MODE>>PSHIFT)&MASK) == 3){
      b.x += t.x * 2.f * kin.x / L.x;
      b.y += t.x * 2.f * kin.y / L.y;
    }
    D4<MODE>(t,b,L,wall_length,V,A,dt,dthalf);
#endif //else of STEEPEST
    __syncthreads();

    /*-- convert molecules to atom --*/
    for(unsigned int k=0;k<nsite;k++){
      for(unsigned int j=tid;j<nmol;j+=bdm){
	const int4   r = r_rep[j];
	const qtype4 q = q_rep[j];
	const int  ind = nsite*r.w + k;

	const qtype4 sq = {q.x*q.x, q.y*q.y, q.z*q.z, q.w*q.w};
	const qtype  xy = 2.f*q.x*q.y, xz = 2.f*q.x*q.z, xw = 2.f*q.x*q.w;
	const qtype  yz = 2.f*q.y*q.z, yw = 2.f*q.y*q.w, zw = 2.f*q.z*q.w;
	qtype4 ar;
	ar.x  = (sq.x + sq.y - sq.z - sq.w)*(qtype)atm[ind].x;
	ar.x += (yz - xw)*(qtype)atm[ind].y;
	ar.x += (yw + xz)*(qtype)atm[ind].z;

	ar.y  = (yz + xw)*(qtype)atm[ind].x;
	ar.y += (sq.x - sq.y + sq.z - sq.w)*(qtype)atm[ind].y;
	ar.y += (zw - xy)*(qtype)atm[ind].z;

	ar.z  = (yw - xz)*(qtype)atm[ind].x;
	ar.z += (zw + xy)*(qtype)atm[ind].y;
	ar.z += (sq.x - sq.y - sq.z + sq.w)*(qtype)atm[ind].z;

#ifdef SWITCHING
	f3i1 *ar_repf = (f3i1*)ar_rep;
	ar_repf[k*nmol+j].x = (float)ar.x;
	ar_repf[k*nmol+j].y = (float)ar.y;
	ar_repf[k*nmol+j].z = (float)ar.z;
#else  //SWITCHING
	ar.x /= (qtype)L.x;
	ar.y /= (qtype)L.y;
	ar.z /= (qtype)L.z;

	ar_rep[k*nmol+j].x = (int)(ar.x*(qtype)IntMax)+r.x;
	ar_rep[k*nmol+j].y = (int)(ar.y*(qtype)IntMax)+r.y;
	ar_rep[k*nmol+j].z = (int)(ar.z*(qtype)IntMax)+r.z;
#endif //SWITCHING 
	ar_rep[k*nmol+j].w = ind;
      }
    }
    for(unsigned int k=0;k<nsite;k++){
      for(unsigned int j=tid;j<nmol;j+=bdm){
	af_rep[k*nmol+j] = make_float4(0.f);
      }
    }
    __syncthreads();
    /*-- force calculation --*/
#ifdef SWITCHING

#ifndef TUNING
    potvir = CalcForcePotVirSwitchingGPU2<MODE>
      (r_rep,f_rep,(f3i1*)ar_rep,af_rep,rslt_dev[bid],L,rc02,rswitch,wall_length,nmol,nsite);
#else
    potvir = CalcForcePotVirSwitchingTIP4P<MODE>(r_rep,f_rep,(f3i1*)ar_rep,af_rep,rslt_dev[bid],L,rc02,rswitch,wall_length,nmol,nsite);
#endif

#else //SWITCHING

    potvir = CalcForcePotVirGPU<MODE>
      (ar_rep,af_rep,kr_rep,kvc_rep,
       L,rc02,alpha,nwave,nsite*nmol,wall_length);

#endif //SWITCHING
    __syncthreads();

    /*-- convert atoms to molecules --*/
    float3 vir = make_float3(0.f);
    for(unsigned int j=tid;j<nmol;j+=bdm){
      float4 f   = make_float4(0.f);
      float4 n   = make_float4(0.f);
      int4   r_i = r_rep[j];
      for(unsigned int k=0;k<nsite;k++){
	const unsigned int ind = nsite*r_rep[j].w + k;
	const float4 af = af_rep[nmol*k+j];
	const qtype4 q  = q_rep[j];
	f.x += af.x;
	f.y += af.y;
	f.z += af.z;
	const qtype4 sq = {q.x*q.x, q.y*q.y, q.z*q.z, q.w*q.w};
	const qtype xy = 2.f*q.x*q.y, xz = 2.f*q.x*q.z, xw = 2.f*q.x*q.w;
	const qtype yz = 2.f*q.y*q.z, yw = 2.f*q.y*q.w, zw = 2.f*q.z*q.w;
	//float3 fb; // convert space force to body force
	double3 fb; // convert space force to body force
	fb.x  = (sq.x + sq.y - sq.z - sq.w) * af.x;
	fb.x += (yz + xw) * af.y;
	fb.x += (yw - xz) * af.z;

	fb.y  = (yz - xw) * af.x;
	fb.y += (sq.x - sq.y + sq.z - sq.w) * af.y;
	fb.y += (zw + xy) * af.z;

	fb.z  = (yw + xz) * af.x;
	fb.z += (zw - xy) * af.y;
	fb.z += (sq.x - sq.y - sq.z + sq.w) * af.z;

	n.x += atm[ind].x * fb.x + atm[ind].y * fb.y + atm[ind].z * fb.z;
	n.y += atm[ind].y * fb.z - atm[ind].z * fb.y;
	n.z += atm[ind].z * fb.x - atm[ind].x * fb.z;
	n.w += atm[ind].x * fb.y - atm[ind].y * fb.x;
#ifndef SWITCHING
	// difference between atomic pressure and molecular pressure
	const int4   ar = ar_rep[nmol*k+j];
	float3 r_ai = make_float3((float)(ar.x - r_i.x)/(float)IntMax*L.x,
				  (float)(ar.y - r_i.y)/(float)IntMax*L.y,
				  (float)(ar.z - r_i.z)/(float)IntMax*L.z);
	vir.x += r_ai.x * af.x;
	vir.y += r_ai.y * af.y;
	vir.z += r_ai.z * af.z;
#endif  //SWITCHING
      }
#ifndef SWITCHING
      f_rep[j] = f;
#endif  //SWITCHING
      n_rep[j] = n;
    }
    // reduce virial correction
    vir = block_all_reduce(vir);
    potvir.x -= vir.x;
    potvir.y -= vir.y;
    potvir.z -= vir.z;

#ifdef STEEPEST
#else
    // integration after force calculation
    for(unsigned int j=tid;j<nmol;j+=bdm){
      // load to register
      float4 v = v_rep[j];
      ptype4 p = p_rep[j];
      const float4 f = f_rep[j];
      const qtype4 q = q_rep[j];
      const float4 n = n_rep[j];
      D5(v,f,p,q,n,t,L,dt,dthalf);

      // store to grobal memory
      v_rep[j] = v;
      p_rep[j] = p;
    }

    if(((MODE>>TSHIFT)&MASK) == 1){
      t.y -= potvir.w*dthalf; // D5 integration
      if(((MODE>>PSHIFT)&MASK) > 0){
	t.y -= P * V * dthalf;
      }
    }
    if(((MODE>>PSHIFT)&MASK) == 1){
      b.x += t.x*((potvir.x+potvir.y+potvir.z)/(3.f*V) - P)*dthalf;
    }
    if(((MODE>>PSHIFT)&MASK) == 2){
      b.z += t.x * potvir.z / L.z * dthalf;
      b.z -= t.x * P * A * dthalf;
    }
    if(((MODE>>PSHIFT)&MASK) == 3){
      if(((MODE>>CSHIFT)&MASK) == 2){
	b.x += t.x*(potvir.x/L.x - P * L.y*wall_length)*dthalf;
	b.y += t.x*(potvir.y/L.y - P * L.x*wall_length)*dthalf;
      }else{
	b.x += t.x*(potvir.x/L.x - P * L.y*L.z)*dthalf;
	b.y += t.x*(potvir.y/L.y - P * L.x*L.z)*dthalf;
      }
    }
    if(((MODE>>TSHIFT)&MASK) == 1){
      D6(t,dthalf);
    }
#endif //STEEPEST
    __syncthreads();
  }
#ifdef STEEPEST
  kin = make_float4(0.f,0.f,0.f,0.f);
#else
  kin = CalcKin(r_rep,v_rep,q_rep,p_rep,L,t,nmol);
#endif
  float tot = potvir.w + sum(kin);
  //thermostat energy
  float2 etst = make_float2(0.f);
  if(((MODE>>TSHIFT)&MASK) == 1){
    etst.x = gkT*log(t.x);
    etst.y = t.y*t.y*t.z;
  }
  tot += sum(etst);

  //barostat energy
  float2 ebst = make_float2(0.f);
  if(((MODE>>PSHIFT)&MASK) == 1)
    ebst.x += (b.x*b.x + b.y*b.y + b.z*b.z)*b.w;
  if(((MODE>>PSHIFT)&MASK) == 2)
    ebst.x += b.z*b.z*b.w;
  if(((MODE>>PSHIFT)&MASK) == 3)
    ebst.x += (b.x*b.x + b.y*b.y)*b.w;
  if(((MODE>>PSHIFT)&MASK) > 0) ebst.y += P*V;
  tot += sum(ebst);

  if(tid==0){
    // store results
    rslt_dev[bid].pot = potvir.w;
    rslt_dev[bid].vol = (float)V;

    rslt_dev[bid].vir.x = potvir.x;
    rslt_dev[bid].vir.y = potvir.y;
    rslt_dev[bid].vir.z = potvir.z;

    rslt_dev[bid].pres.x = (2.f*kin.x + potvir.x) / V;
    rslt_dev[bid].pres.y = (2.f*kin.y + potvir.y) / V;
    rslt_dev[bid].pres.z = (2.f*kin.z + potvir.z) / V;

    rslt_dev[bid].tra = kin.w;
    rslt_dev[bid].rot = kin.x + kin.y + kin.z;
    rslt_dev[bid].kin = sum(kin);

    rslt_dev[bid].tst = etst;
    rslt_dev[bid].bst = ebst;
    rslt_dev[bid].tot = tot;
    rslt_dev[bid].drf = tot - H0;

    rslt_dev[bid].temp = 2.f*sum(kin) / (6.f*(float)nmol);

    rslt_dev[bid].t.x = t.x;
    rslt_dev[bid].t.y = t.y;

    rslt_dev[bid].b.x = b.x;
    rslt_dev[bid].b.y = b.y;
    rslt_dev[bid].b.z = b.z;
    // store vairables to continue simulation
    tst_dev[bid] = t;
    //bst_dev[bid] = b;
    bst_dev[bid] = make_float4(b);

    potvir_dev[bid] = potvir;
    L_dev[bid] = L;
    H0_dev[bid] = H0;
    cond1_dev[bid] = cond0_dev[bid];
  }
}

//#define TEST
#ifdef TEST
int main(int argc,char** argv){
  Parameter param;
  param.read(std::cin);
  param.print(std::cout);

  const unsigned int nreplica = param.nreplica;

  std::ostream *ocpu;
  ocpu = new std::ofstream("cpu.txt",std::ios_base::out);

  Molecules **mlcls;
  IOManager **io;
  mlcls = new Molecules*[nreplica];
  io    = new IOManager*[nreplica];
  for(int i=0;i<nreplica;i++){
    mlcls[i] = new Molecules(param);
    io[i]    = new IOManager(param);
    io[i]->ReadInputs(mlcls[i]);
    mlcls[i]->ExecuteSteps();
    mlcls[i]->PrintAll(*ocpu);
  }

  double tmp[nreplica];
  for(unsigned int rep=0;rep<nreplica;rep++){
    tmp[rep] = 1.0;
  }
  REMDGPU *remd;
  remd = new REMDGPU(mlcls,nreplica);
  remd->ExecuteSteps(tmp,tmp);

  delete ocpu;
}
#endif
