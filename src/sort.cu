#include "quicksort.cuh"
#include "bitonic_sort.cuh"
#include "remdgpu.h"

#define QSORT

__device__ int val_dev[102400],key_dev[102400];

__global__
void sort(int4   *r_dev,float4 *v_dev,float4 *f_dev,
	  qtype4 *q_dev,ptype4 *p_dev,float4 *n_dev,int nelem){
  int4   *r_rep = r_dev + blockDim.x*blockIdx.x;
  float4 *v_rep = v_dev + blockDim.x*blockIdx.x;
  float4 *f_rep = f_dev + blockDim.x*blockIdx.x;
  qtype4 *q_rep = q_dev + blockDim.x*blockIdx.x;
  ptype4 *p_rep = p_dev + blockDim.x*blockIdx.x;
  float4 *n_rep = n_dev + blockDim.x*blockIdx.x;

  int *key = key_dev + blockDim.x*blockIdx.x;
  int *val = val_dev + blockDim.x*blockIdx.x;

  key[threadIdx.x] = r_rep[threadIdx.x].z;
  val[threadIdx.x] = threadIdx.x;
#ifndef QSORT
  if(nelem > blockDim.x){
    key[blockDim.x+threadIdx.x] =  2147483647;
    val[blockDim.x+threadIdx.x] = blockDim.x+threadIdx.x;
  }
#endif
  __syncthreads();
  if(threadIdx.x == 0){
#ifdef QSORT
    quicksort_key<<<1,blockDim.x>>>(key,val,0,blockDim.x-1);
#else
    int nthread = 2;
    while(nthread < nelem) nthread <<= 1;
    switch(nthread){
    case(2):
      bitonic_sort_key<   2,int><<<1,   2>>>(key,val);
      break;
    case(4):
      bitonic_sort_key<   4,int><<<1,   4>>>(key,val);
      break;
    case(8):
      bitonic_sort_key<   8,int><<<1,   8>>>(key,val);
      break;
    case(16):
      bitonic_sort_key<  16,int><<<1,  16>>>(key,val);
      break;
    case(32):
      bitonic_sort_key<  32,int><<<1,  32>>>(key,val);
      break;
    case(64):
      bitonic_sort_key<  64,int><<<1,  64>>>(key,val);
      break;
    case(128):
      bitonic_sort_key< 128,int><<<1, 128>>>(key,val);
      break;
    case(256):
      bitonic_sort_key< 256,int><<<1, 256>>>(key,val);
      break;
    case(512):
      bitonic_sort_key< 512,int><<<1, 512>>>(key,val);
      break;
    case(1024):
      bitonic_sort_key<1024,int><<<1,1024>>>(key,val);
      break;
      /*
	case(2048):
	bitonic_sort_key<2048,int><<<1,2048>>>(key,val);
	break;
      //*/
    }
#endif
    cudaDeviceSynchronize();
  }
  __syncthreads();
  const int index = val[threadIdx.x];
  //if(index>=nelem) printf("error: dummy is included! (%d %d)\n",threadIdx.x,index);
  //printf("index of %d is %d\n",threadIdx.x,index);
  const int4   rtmp = r_rep[index];
  const float4 vtmp = v_rep[index];
  const float4 ftmp = f_rep[index];
  const qtype4 qtmp = q_rep[index];
  const ptype4 ptmp = p_rep[index];
  const float4 ntmp = n_rep[index];
  __syncthreads();
  r_rep[threadIdx.x] = rtmp;
  v_rep[threadIdx.x] = vtmp;
  f_rep[threadIdx.x] = ftmp;
  q_rep[threadIdx.x] = qtmp;
  p_rep[threadIdx.x] = ptmp;
  n_rep[threadIdx.x] = ntmp;
  __syncthreads();
}
