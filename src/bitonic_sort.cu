
#include "bitonic_sort.cuh"

template <class X>
__device__ __forceinline__
void swap(X &a,X&b){
  const X tmp = a; a = b; b = tmp;
}

template <class X>
__device__ __forceinline__
void swap(X &key,X &val,int mask,int dir){
  const X y = __shfl_xor(key,mask);
  if(key < y == dir){
    key = y;
    val = __shfl_xor(val,mask);
  }
}

//extract k th bit of i
__device__ inline
 int bfe(int i, int k){
  return (i>>k) & 0x01;
}

template <int N, class X>
__global__
void bitonic_sort_key(X *gkeys,X *gvals){
  X *keys = gkeys + blockDim.x*blockIdx.x;
  X *vals = gvals + blockDim.x*blockIdx.x;
  int i = threadIdx.x;
#if 0
  for(int k=2;k <= N;k<<=1){
    for(int j=k>>1;j>0;j>>=1){
      int ixj = i^j;
      if ((ixj)>i) {
	if ((i&k)==0) {
	  /* Sort ascending */
	  if (keys[i]>keys[ixj]) {
	    swap(keys[i],keys[ixj]);
	    swap(vals[i],vals[ixj]);
	  }
	}
	if ((i&k)!=0) {
	  /* Sort descending */
	  if (keys[i]<keys[ixj]) {
	    swap(keys[i],keys[ixj]);
	    swap(vals[i],vals[ixj]);
	  }
	}
      }
      __syncthreads();
    }
  }
#else
  const int laneid = threadIdx.x % WarpSize;
  //const int warpid = threadIdx.x / WarpSize;
  X key = keys[i];
  X val = vals[i];

  if(N>=2){
    swap(key, val, 0x01, bfe(laneid, 1) ^ bfe(laneid, 0)); //  2
  }
  if(N>=4){
    swap(key, val, 0x02, bfe(laneid, 2) ^ bfe(laneid, 1)); //  4
    swap(key, val, 0x01, bfe(laneid, 2) ^ bfe(laneid, 0));
  }
  if(N>=8){
    swap(key, val, 0x04, bfe(laneid, 3) ^ bfe(laneid, 2)); //  8
    swap(key, val, 0x02, bfe(laneid, 3) ^ bfe(laneid, 1));
    swap(key, val, 0x01, bfe(laneid, 3) ^ bfe(laneid, 0));
  }
  if(N>=16){
    swap(key, val, 0x08, bfe(laneid, 4) ^ bfe(laneid, 3)); // 16
    swap(key, val, 0x04, bfe(laneid, 4) ^ bfe(laneid, 2));
    swap(key, val, 0x02, bfe(laneid, 4) ^ bfe(laneid, 1));
    swap(key, val, 0x01, bfe(laneid, 4) ^ bfe(laneid, 0));
  }
  if(N>=32){
    swap(key, val, 0x10, bfe(laneid, 4)); // 32
    swap(key, val, 0x08, bfe(laneid, 3));
    swap(key, val, 0x04, bfe(laneid, 2));
    swap(key, val, 0x02, bfe(laneid, 1));
    swap(key, val, 0x01, bfe(laneid, 0));
  }
  __shared__ X skey[N];
  __shared__ X sval[N];
  skey[i] = key;
  sval[i] = val;
  __syncthreads();
  for(int k=32;k <= N;k<<=1){
    for(int j=k>>1;j>0;j>>=1){
      int ixj = i^j;
      if ((ixj)>i) {
	if ((i&k)==0) {
	  /* Sort ascending */
	  if (skey[i]>skey[ixj]) {
	    swap(skey[i],skey[ixj]);
	    swap(sval[i],sval[ixj]);
	  }
	}
	if ((i&k)!=0) {
	  /* Sort descending */
	  if (skey[i]<skey[ixj]) {
	    swap(skey[i],skey[ixj]);
	    swap(sval[i],sval[ixj]);
	  }
	}
      }
      __syncthreads();
    }
  }
  keys[i] = skey[i];
  vals[i] = sval[i];
#endif
}

template __global__ void bitonic_sort_key<   2,int>(int*,int*);
template __global__ void bitonic_sort_key<   4,int>(int*,int*);
template __global__ void bitonic_sort_key<   8,int>(int*,int*);
template __global__ void bitonic_sort_key<  16,int>(int*,int*);
template __global__ void bitonic_sort_key<  32,int>(int*,int*);
template __global__ void bitonic_sort_key<  64,int>(int*,int*);
template __global__ void bitonic_sort_key< 128,int>(int*,int*);
template __global__ void bitonic_sort_key< 256,int>(int*,int*);
template __global__ void bitonic_sort_key< 512,int>(int*,int*);
template __global__ void bitonic_sort_key<1024,int>(int*,int*);
//template __global__ void bitonic_sort_key<2048,int>(int*,int*); // there is a bug for 2048 element
