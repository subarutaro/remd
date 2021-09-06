#ifndef CUH_4SITE
#define CUH_4SITE

#include "remdgpu.h"

#ifndef WarpSize
#define WarpSize 32
#endif

/*
const float3 r_OW  = { 6.556274e-02f, 0.0000000e+00f, 0.0000000e+00f};
const float3 r_HW1 = {-5.203197e-01f, 7.5695022e-01f, 0.0000000e+00f};
const float3 r_HW2 = {-5.203197e-01f,-7.5695022e-01f, 0.0000000e+00f};
const float3 r_MW  = {-8.443726e-02f, 0.0000000e+00f, 0.0000000e+00f};

const float4 mass = {1.126000e+00f, 7.219740e-02f, 3.841210e-02f, 1.106090e-01f};
const float sigma = 3.15365e+00f;
const float epsilon = 6.485200e-01f;
const float4 lj = {4.f*epsilon*powf(sigma,12), 4.f*epsilon*powf(sigma,6), 48.f*epsilon*powf(sigma,12), 24.f*epsilon*powf(sigma,6)};

const float q_HW = 5.2000e-01f;
const float q_MW = -1.0400e+00f;
const float3 q_all = {q_HW*q_HW, q_HW*q_MW, q_MW*q_MW};
//*/

#ifdef TIP3P
//4.f * 0.6501940 * 3.16557f^12 and 4.f * 3.16557f * 0.6501940^6
__device__ const float lje12  = 2633455.63681f;
__device__ const float lje06  = 2617.06480953f;

__device__ const float q_HO = -0.8476f*0.4238f*1389.35080129f;
__device__ const float q_HH =  0.4238f*0.4238f*1389.35080129f;
__device__ const float q_OO =  0.8476f*0.8476f*1389.35080129f;
#endif

// 1 / (4 Pi e0 e* l*) = 1389.35080129
#ifdef TIP4P
__device__ const float lje12 = 2510413.58406f;
__device__ const float lje06 = 2551.90393043f;

__device__ const float q_HO = -0.5408f*1389.35080129f;
__device__ const float q_HH =  0.2704f*1389.35080129f;
__device__ const float q_OO =  1.0816f*1389.35080129f;
#endif

#ifdef TIP4PICE
//4.f * 0.8821644f * 3.16680f^12 and 4.f * 0.8821644f * 3.16680f^6
__device__ const float lje12 = 3589691.72791f;
__device__ const float lje06 = 3559.0438319f;

__device__ const float q_HO = -0.5897f*1.1794f*1389.35080129f;
__device__ const float q_HH =  0.5897f*0.5897f*1389.35080129f;
__device__ const float q_OO =  1.1794f*1.1794f*1389.35080129f;
#endif

__device__ const float ljf12  = 12.f * lje12;
__device__ const float ljf06  =  6.f * lje06;

__device__ __forceinline__
f3i1 __shfl(const f3i1 a,const int index){
  f3i1 val;
  val.x = __shfl(a.x,index);
  val.y = __shfl(a.y,index);
  val.z = __shfl(a.z,index);
  //val.w = __shfl(a.w,index);
  return val;
}
__device__ __forceinline__
f3i1 __shfl_xor(const f3i1 a,const int mask){
  f3i1 val;
  val.x = __shfl_xor(a.x,mask);
  val.y = __shfl_xor(a.y,mask);
  val.z = __shfl_xor(a.z,mask);
  //val.w = __shfl_xor(a.w,mask);
  return val;
}

__device__ __forceinline__
void CalcLJ(const float &gx,const float &gy,const float &gz,const f3i1 &ai,const f3i1 &aj,float4  &fi,float4  &gfi,float &e,const float &sw){
  const float x = gx + ai.x - aj.x;
  const float y = gy + ai.y - aj.y;
  const float z = gz + ai.z - aj.z;
  const float r02  = x*x + y*y + z*z;
  const float r02i = 1.f/r02;
  const float r06i = r02i*r02i*r02i;
  const float ftmp = r06i * (ljf12*r06i - ljf06) * r02i * sw;
  fi.x += ftmp*x;
  fi.y += ftmp*y;
  fi.z += ftmp*z;
  gfi.x += ftmp*x;
  gfi.y += ftmp*y;
  gfi.z += ftmp*z;

  const float etmp = r06i * (lje12*r06i - lje06);
  e += etmp;
}

__device__ __forceinline__
void CalcCoulomb(const float &gx,const float &gy,const float &gz,const f3i1 &ai,const f3i1 &aj,const float &qq,float4 &fi,float4 &gfi,float &e,const float &sw){
  const float x = gx + ai.x - aj.x;
  const float y = gy + ai.y - aj.y;
  const float z = gz + ai.z - aj.z;
  const float r02   = x*x + y*y + z*z;
  const float r02i  = 1.f/r02;
  const float r01i = sqrtf(r02i);

  const float ftmp = qq*r01i*r02i*sw;
  fi.x  += ftmp*x;
  fi.y  += ftmp*y;
  fi.z  += ftmp*z;
  gfi.x += ftmp*x;
  gfi.y += ftmp*y;
  gfi.z += ftmp*z;

  const float etmp = qq*r01i;
  e += etmp;
}

#define CALC_FORCE_NEWTON3(width)					\
  float4 force_i = make_float4(0.f);					\
  const float gx = (float)(gcoor_i.x - __shfl(gcoor_i.x,(laneid+j)%width))*ItoF.x; \
  const float gy = (float)(gcoor_i.y - __shfl(gcoor_i.y,(laneid+j)%width))*ItoF.y; \
  const float gz = (float)(gcoor_i.z - __shfl(gcoor_i.z,(laneid+j)%width))*ItoF.z; \
  const float r02 = gx*gx + gy*gy + gz*gz;				\
  const float r01 = sqrtf(r02);						\
  const float2 switching = tex1D(tex_switch,(r01-rl)*rcli);		\
  const float sw  = switching.x;					\
  const float dsw = switching.y;					\
  const float r01i = 1.f/r01;						\
  const float dswx = dsw*gx*r01i;					\
  const float dswy = dsw*gy*r01i;					\
  const float dswz = dsw*gz*r01i;					\
  float lje = 0.f,clmbe = 0.f;						\
  CalcLJ(gx,gy,gz,ai[0],__shfl(ai[0],(laneid+j)%width),fi[0],force_i,lje,sw); \
  CalcCoulomb(gx,gy,gz,ai[1],__shfl(ai[1],(laneid+j)%width),q_HH,fi[1],force_i,clmbe,sw); \
  CalcCoulomb(gx,gy,gz,ai[1],__shfl(ai[2],(laneid+j)%width),q_HH,fi[1],force_i,clmbe,sw); \
  CalcCoulomb(gx,gy,gz,ai[1],__shfl(ai[3],(laneid+j)%width),q_HO,fi[1],force_i,clmbe,sw); \
  CalcCoulomb(gx,gy,gz,ai[2],__shfl(ai[1],(laneid+j)%width),q_HH,fi[2],force_i,clmbe,sw); \
  CalcCoulomb(gx,gy,gz,ai[2],__shfl(ai[2],(laneid+j)%width),q_HH,fi[2],force_i,clmbe,sw); \
  CalcCoulomb(gx,gy,gz,ai[2],__shfl(ai[3],(laneid+j)%width),q_HO,fi[2],force_i,clmbe,sw); \
  CalcCoulomb(gx,gy,gz,ai[3],__shfl(ai[1],(laneid+j)%width),q_HO,fi[3],force_i,clmbe,sw); \
  CalcCoulomb(gx,gy,gz,ai[3],__shfl(ai[2],(laneid+j)%width),q_HO,fi[3],force_i,clmbe,sw); \
  CalcCoulomb(gx,gy,gz,ai[3],__shfl(ai[3],(laneid+j)%width),q_OO,fi[3],force_i,clmbe,sw); \
  const float e = lje + clmbe;						\
  force_i.x -= e * dswx;						\
  force_i.y -= e * dswy;						\
  force_i.z -= e * dswz;						\
  force_i_sum += force_i;

#define CALC_FORCE()							\
  float4 force_i = make_float4(0.f);					\
  const float gx = (float)(gcoor_i.x - __shfl(gcoor_j.x,j))*ItoF.x;	\
  const float gy = (float)(gcoor_i.y - __shfl(gcoor_j.y,j))*ItoF.y;	\
  const float gz = (float)(gcoor_i.z - __shfl(gcoor_j.z,j))*ItoF.z;	\
  const float r02 = gx*gx + gy*gy + gz*gz;				\
  const float r01 = sqrtf(r02);						\
  const float2 switching = tex1D(tex_switch,(r01-rl)*rcli);		\
  const float sw  = switching.x;					\
  const float dsw = switching.y;					\
  const float r01i = 1.f/r01;						\
  const float dswx = dsw*gx*r01i;					\
  const float dswy = dsw*gy*r01i;					\
  const float dswz = dsw*gz*r01i;					\
  float lje = 0.f,clmbe = 0.f;						\
  CalcLJ(gx,gy,gz,ai[0],__shfl(oxy_j,j),fi[0],force_i,lje,sw);		\
  CalcCoulomb(gx,gy,gz,ai[1],__shfl(hyd0_j,j),q_HH,fi[1],force_i,clmbe,sw); \
  CalcCoulomb(gx,gy,gz,ai[1],__shfl(hyd1_j,j),q_HH,fi[1],force_i,clmbe,sw); \
  CalcCoulomb(gx,gy,gz,ai[1],__shfl(dum_j, j),q_HO,fi[1],force_i,clmbe,sw); \
  CalcCoulomb(gx,gy,gz,ai[2],__shfl(hyd0_j,j),q_HH,fi[2],force_i,clmbe,sw); \
  CalcCoulomb(gx,gy,gz,ai[2],__shfl(hyd1_j,j),q_HH,fi[2],force_i,clmbe,sw); \
  CalcCoulomb(gx,gy,gz,ai[2],__shfl(dum_j, j),q_HO,fi[2],force_i,clmbe,sw); \
  CalcCoulomb(gx,gy,gz,ai[3],__shfl(hyd0_j,j),q_HO,fi[3],force_i,clmbe,sw); \
  CalcCoulomb(gx,gy,gz,ai[3],__shfl(hyd1_j,j),q_HO,fi[3],force_i,clmbe,sw); \
  CalcCoulomb(gx,gy,gz,ai[3],__shfl(dum_j, j),q_OO,fi[3],force_i,clmbe,sw); \
  const float e = lje + clmbe;						\
  force_i.x -= e * dswx;						\
  force_i.y -= e * dswy;						\
  force_i.z -= e * dswz;						\
  force_i_sum += force_i;

#define CALC_FORCE_CONFINED()			\
  const float x = gcoor_i.x*ItoF.x + ai[0].x;	\
  const float y = gcoor_i.y*ItoF.y + ai[0].y;	\
  const float r02 = x*x + y*y;			\
  const float r01 = sqrtf(r02);			\
  float2 table  = tex1D(tex_wall,r01/R);	\
  fi[0].x       += table.x * x;			\
  fi[0].y       += table.x * y;			\
  force_i_sum.x += table.x * x;			\
  force_i_sum.y += table.x * y;

#define divup(a,b) ((a + b -1) / b)

template <int MODE>
__device__
float4 CalcForcePotVirSwitchingTuned
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
  const int warpid = threadIdx.x / WarpSize;
  const int nwarp  = blockDim.x / WarpSize;
  const int laneid = threadIdx.x % WarpSize;
  const int nblock = divup(nmol,WarpSize);
  const int nblock_j = 1;

  const float ImaxInv = 2.3283064365387e-10f;
  const float3 ItoF = {L.x*ImaxInv,L.y*ImaxInv,L.z*ImaxInv};

  // switching function variables
  const float rc  = sqrtf(rcutSq);
  const float rl  = rswitch;
  const float rlc = rl - rc;
  const float rcli = 1.f/(rc - rl);

  float3 vir = make_float3(0.f);
  float lj_sum=0.f,clmb_sum=0.f,wall_sum=0.f;

  const int nshared = NTHREADS+2*nblock_j*WarpSize;
  __shared__ int4 gcoor_shared[nshared];
  __shared__ f3i1 acoor_shared[4*nshared];

  int blockid = warpid;
  for(int i=0;i<divup(nmol,bdm);i++){
    const int ind_i = blockid*WarpSize + laneid;
    int ind_i_round = ind_i;
    if(ind_i_round >= nmol) ind_i_round -= nmol;

    const int4 gcoor_i = gcoor[ind_i_round];
    f3i1 ai[4];
    float4 fi[4];

    ai[0] = acoor[ind_i_round];
    ai[1] = acoor[ind_i_round + nmol];
    ai[2] = acoor[ind_i_round + 2*nmol];
    ai[3] = acoor[ind_i_round + 3*nmol];

    fi[0] = make_float4(0.f);
    fi[1] = make_float4(0.f);
    fi[2] = make_float4(0.f);
    fi[3] = make_float4(0.f);

    gcoor_shared[tid + nblock_j*WarpSize] = gcoor_i;
    acoor_shared[tid + nblock_j*WarpSize + nshared*0] = ai[0];
    acoor_shared[tid + nblock_j*WarpSize + nshared*1] = ai[1];
    acoor_shared[tid + nblock_j*WarpSize + nshared*2] = ai[2];
    acoor_shared[tid + nblock_j*WarpSize + nshared*3] = ai[3];

    if(tid < nblock_j*WarpSize){
      int index = ind_i - nblock_j*WarpSize;
      if(index<0) index += nmol;
      gcoor_shared[tid] = gcoor[index];
      acoor_shared[tid + nshared*0] = acoor[index + nmol*0];
      acoor_shared[tid + nshared*1] = acoor[index + nmol*1];
      acoor_shared[tid + nshared*2] = acoor[index + nmol*2];
      acoor_shared[tid + nshared*3] = acoor[index + nmol*3];

      index = ind_i + bdm;
      if(index>=nmol) index -= nmol;
      int sindex = nblock_j*WarpSize + bdm + tid;
      gcoor_shared[sindex] = gcoor[index];
      acoor_shared[sindex + nshared*0] = acoor[index + nmol*0];
      acoor_shared[sindex + nshared*1] = acoor[index + nmol*1];
      acoor_shared[sindex + nshared*2] = acoor[index + nmol*2];
      acoor_shared[sindex + nshared*3] = acoor[index + nmol*3];
    }
    __syncthreads();
    if(blockid >= nblock) continue;

    float4 force_i_sum = make_float4(0.f);
    for(int w=-nblock_j;w<=nblock_j;w++){
      int wid_j = blockid + w;
      if(wid_j <  0)      wid_j += nblock;
      if(wid_j >= nblock) wid_j -= nblock;

      if(wid_j == blockid){
	#pragma unroll 31
	for(int j=1;j<WarpSize;j++){
	  CALC_FORCE_NEWTON3(WarpSize);
	  if(ind_i < nmol){
	    lj_sum += lje*sw;
	    clmb_sum += clmbe*sw;
	    vir.x += force_i.x * gx;
	    vir.y += force_i.y * gy;
	    vir.z += force_i.z * gz;
	  }
	}
      }else{
#if 0
	const int  ind_j   = WarpSize*wid_j + laneid;
	const int4 gcoor_j = gcoor[ind_j];
	const f3i1 oxy_j   = acoor[ind_j + 0*nmol];
	const f3i1 hyd0_j  = acoor[ind_j + 1*nmol];
	const f3i1 hyd1_j  = acoor[ind_j + 2*nmol];
	const f3i1 dum_j   = acoor[ind_j + 3*nmol];
#else
	const int  ind_j   = (nblock_j + w) * WarpSize + tid;
	int ind_j_old = WarpSize*wid_j + laneid;
	if(ind_j_old >= nmol) ind_j_old -= nmol;
	const int4 gcoor_j = gcoor_shared[ind_j];
	const f3i1 oxy_j   = acoor_shared[ind_j + 0*nshared];
	const f3i1 hyd0_j  = acoor_shared[ind_j + 1*nshared];
	const f3i1 hyd1_j  = acoor_shared[ind_j + 2*nshared];
	const f3i1 dum_j   = acoor_shared[ind_j + 3*nshared];
#endif
	#pragma unroll 32
	for(int j=0;j<WarpSize;j++){
	  CALC_FORCE();
	  if(ind_i < nmol){
	    lj_sum += lje*sw;
	    clmb_sum += clmbe*sw;
	    vir.x += force_i.x * gx;
	    vir.y += force_i.y * gy;
	    vir.z += force_i.z * gz;
	  }
	}// j roop
      }// if w == 0
    }// w roop
    // confined force and potential
    CALC_FORCE_CONFINED();
    if(ind_i < nmol){
      wall_sum  += table.y;
      gforce[ind_i] = force_i_sum;
      aforce[ind_i]        = fi[0];
      aforce[ind_i+nmol]   = fi[1];
      aforce[ind_i+2*nmol] = fi[2];
      aforce[ind_i+3*nmol] = fi[3];
    }
    blockid += nwarp;
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
#endif
