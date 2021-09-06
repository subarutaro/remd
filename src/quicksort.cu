#include "scan.cuh"
#include "reduce.cuh"

#include "quicksort.cuh"

template <class X>
__device__ __forceinline__
void swap(X &a,X &b){
  X tmp = a; a = b; b = tmp;
}

// reduce and scan are limited as maxixmum # of threads = 1024
template <class X>
__global__
void quicksort(X *gdata,int left,int right,int depth){
  if(left == right) return;
  const int tid = threadIdx.x;

  X *data = gdata + blockDim.x*blockIdx.x + left;
  const int ndata = right - left + 1;
  X input = data[tid];
  X pivot = data[ndata/2];
  __syncthreads();

  int lmask = input <  pivot ? 1 : 0;
  int lsum  = block_inclusive_scan(lmask) - 1;
  int nleft = block_all_reduce(lmask);

  int emask = input == pivot ? 1 : 0;
  int esum  = block_inclusive_scan(emask) - 1;
  int neven = block_all_reduce(emask);

  int rmask = input >  pivot ? 1 : 0;
  int rsum  = block_inclusive_scan(rmask) - 1;

  if(lmask == 1) data[lsum] = input;
  else if(emask == 1) data[nleft+esum] = input;
  else data[nleft+neven+rsum] = input;
  __syncthreads();
  //*
  printf("%2d %4d %4d %4d %d %4d %d %4d %4d %4d %4d\n",
	 depth,left,right,tid,lmask,lsum,rmask,rsum,nleft,pivot,input);
  //*/
  int lright = left+nleft-1;
  int rleft  = left+nleft+neven;
  if(tid == 0){
    if(lright > left){
      cudaStream_t s1;
      cudaStreamCreateWithFlags(&s1, cudaStreamNonBlocking);
      quicksort<<<1,nleft,32,s1>>>(gdata, left, lright,depth+1);
    }
    if(right > rleft){
      cudaStream_t s2;
      cudaStreamCreateWithFlags(&s2, cudaStreamNonBlocking);
      quicksort<<<1,ndata-(nleft+neven),32,s2>>>(gdata, rleft, right,depth+1);
      cudaDeviceSynchronize();
    }
  }
  __syncthreads();
}

template __global__ void quicksort<int>(int*,int,int,int);

template < class X0, class X1>
__global__
void quicksort_key(X0 *gkeys, X1 *gvalues,int left,int right,int depth){
  if(left == right) return;

  const int tid = blockDim.x * blockIdx.x + threadIdx.x;
  X0 *keys   = gkeys + left;
  X0 *values = gvalues + left;

  const int ndata = right - left + 1;
  X0 key   = keys[tid];
  X0 pivot = keys[ndata/2];
  X1 value = values[tid];
  __syncthreads();

  int lmask = key <  pivot ? 1 : 0;
  int lsum  = block_inclusive_scan(lmask) - 1;
  int nleft = block_all_reduce(lmask);

  int emask = key == pivot ? 1 : 0;
  int esum  = block_inclusive_scan(emask) - 1;
  int neven = block_all_reduce(emask);

  int rmask = key >  pivot ? 1 : 0;
  int rsum  = block_inclusive_scan(rmask) - 1;

  if(lmask == 1){
    keys[lsum]   = key;
    values[lsum] = value;
  }else if(emask == 1){
    keys[nleft+esum]   = key;
    values[nleft+esum] = value;
  }else{
    keys[nleft+neven+rsum]   = key;
    values[nleft+neven+rsum] = value;
  }
  __syncthreads();
  /*
  printf("%2d %4d %4d %4d %d %4d %d %4d %4d %4d %4d %4d\n",
	 depth,left,right,tid,lmask,lsum,rmask,rsum,nleft,pivot,key,value);
  //*/
  int lright = left+nleft-1;
  int rleft  = left+nleft+neven;
  if(tid == 0){
    if(lright > left){
      cudaStream_t s1;
      cudaStreamCreateWithFlags(&s1, cudaStreamNonBlocking);
      quicksort_key<<<1,nleft,32,s1>>>(gkeys, gvalues, left, lright,depth+1);
      cudaDeviceSynchronize();
    }
    if(right > rleft){
      cudaStream_t s2;
      cudaStreamCreateWithFlags(&s2, cudaStreamNonBlocking);
      quicksort_key<<<1,ndata-(nleft+neven),32,s2>>>(gkeys, gvalues, rleft, right,depth+1);
      cudaDeviceSynchronize();
    }
  }
  __syncthreads();
  /*
  if(threadIdx.x == 0)
    for(int i=0;i<blockDim.x;i++) printf("%d %d %d %d\n",depth,i,keys[i],values[i]);
  //*/
}

template __global__ void quicksort_key<int,int>(int*,int*,int,int,int);
