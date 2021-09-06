#ifndef H_SCAN
#define H_SCAN

#include <cstdio>

#ifndef WarpSize
#define WarpSize 32
#endif

#ifndef LogWarpSize
#define LogWarpSize 5
#endif

template <class X>
__device__ __forceinline__
X warp_inclusive_scan(X val){
  const int laneid = threadIdx.x % WarpSize;
#if 0
  for(int offstet = 1; offset < WarpSize; offset <<= 1){
    X tmp = __shfl_up(val,offset);
    if(laneid >= offset) val += tmp;
  }
#else
  X tmp;
  tmp = __shfl_up(val,1);
  if(laneid >=  1) val += tmp;
  tmp = __shfl_up(val,2);
  if(laneid >=  2) val += tmp;
  tmp = __shfl_up(val,4);
  if(laneid >=  4) val += tmp;
  tmp = __shfl_up(val,8);
  if(laneid >=  8) val += tmp;
  tmp = __shfl_up(val,16);
  if(laneid >= 16) val += tmp;
#endif
  return val;
}

template <class X>
__device__ __forceinline__
X warp_exclusive_scan(X val){
  const int laneid = threadIdx.x % WarpSize;
  val = __shfl_up(val,1);
  if(laneid==0) val = (X)0;
#if 0
  for(int offstet = 1; offset < WarpSize; offset <<= 1){
    X tmp = __shfl_up(val,offset);
    if(laneid >= offset) val += tmp;
  }
#else
  X tmp;
  tmp = __shfl_up(val,1);
  if(laneid >=  1) val += tmp;
  tmp = __shfl_up(val,2);
  if(laneid >=  2) val += tmp;
  tmp = __shfl_up(val,4);
  if(laneid >=  4) val += tmp;
  tmp = __shfl_up(val,8);
  if(laneid >=  8) val += tmp;
  tmp = __shfl_up(val,16);
  if(laneid >= 16) val += tmp;
#endif
  return val;
}

template <class X>
__device__ __forceinline__
X block_inclusive_scan(X val){
  const int warpid = threadIdx.x / WarpSize;
  const int laneid = threadIdx.x % WarpSize;
  __shared__ X shared[WarpSize];

  val = warp_inclusive_scan(val);

  if(laneid == WarpSize-1) shared[warpid] = val;
  __syncthreads();

  X tmp;
  if(warpid == 0){
    tmp = shared[laneid];
    tmp = warp_exclusive_scan(tmp);
    shared[laneid] = tmp;
  }
  __syncthreads();

  tmp = shared[warpid];
  val += tmp;

  return val;
}

template <class X>
__device__ __forceinline__
X block_exclusive_scan(X val){
  const int warpid = threadIdx.x / WarpSize;
  const int laneid = threadIdx.x % WarpSize;

  X tmp = block_inclusive_scan(val);
  __shared__ X shared[WarpSize];
  if(laneid == WarpSize - 1) shared[warpid] = tmp;
  __syncthreads();
  val = __shfl_up(tmp,1);

  if(laneid == 0){
    if(warpid > 0) val = shared[warpid-1];
    else           val = (X)0;
  }
  return val;
}
#endif

//#define SCAN_TEST
#ifdef SCAN_TEST

__global__
void warp_inclusive_scan_test(int *input,int *output){
  int val = input[threadIdx.x];
  val = warp_inclusive_scan(val);
  output[threadIdx.x] = val;
}
__global__
void warp_exclusive_scan_test(int *input,int *output){
  int val = input[threadIdx.x];
  val = warp_exclusive_scan(val);
  output[threadIdx.x] = val;
}
__global__
void block_inclusive_scan_test(int *input,int *output){
  int val = input[threadIdx.x];
  val = block_inclusive_scan(val);
  output[threadIdx.x] = val;
}
__global__
void block_exclusive_scan_test(int *input,int *output){
  int val = input[threadIdx.x];
  val = block_exclusive_scan(val);
  output[threadIdx.x] = val;
}

int main(int argc,char** argv){
  int *input, *w_in, *w_ex,*b_in,*b_ex;
  int *input_dev, *output_dev;

  const int size = atoi(argv[1]);

  input = (int*)malloc(sizeof(int)*size);
  w_in  = (int*)malloc(sizeof(int)*size);
  w_ex  = (int*)malloc(sizeof(int)*size);
  b_in  = (int*)malloc(sizeof(int)*size);
  b_ex  = (int*)malloc(sizeof(int)*size);
  cudaMalloc((void**) &input_dev,sizeof(int)*size);
  cudaMalloc((void**)&output_dev,sizeof(int)*size);

  for(int i=0;i<size;i++) input[i] = 1;
  cudaMemcpy(input_dev,input,sizeof(int)*size,cudaMemcpyHostToDevice);

  warp_inclusive_scan_test<<<1,size>>>(input_dev,output_dev);
  cudaMemcpy(w_in,output_dev,sizeof(int)*size,cudaMemcpyDeviceToHost);
  warp_exclusive_scan_test<<<1,size>>>(input_dev,output_dev);
  cudaMemcpy(w_ex,output_dev,sizeof(int)*size,cudaMemcpyDeviceToHost);
  block_inclusive_scan_test<<<1,size>>>(input_dev,output_dev);
  cudaMemcpy(b_in,output_dev,sizeof(int)*size,cudaMemcpyDeviceToHost);
  block_exclusive_scan_test<<<1,size>>>(input_dev,output_dev);
  cudaMemcpy(b_ex,output_dev,sizeof(int)*size,cudaMemcpyDeviceToHost);
  for(int i=0;i<size;i++)
    printf("%5d %5d %5d %5d %5d\n",i,w_in[i],w_ex[i],b_in[i],b_ex[i]);
}
#endif
