#ifndef CUH_REDUCE
#define CUH_REDUCE

#ifndef WarpSize
#define WarpSize 32
#endif

#ifndef LogWarpSize
#define LogWarpSize 5
#endif

#include "vector.cuh"

__device__ __forceinline__
int nwarp(){
  return (blockDim.x + WarpSize-1)/WarpSize;
}

template <class X>
__forceinline__ __device__
X warp_reduce(X val){
#if 0
  for(int offset = WarpSize/2; offset > 0; offset /= 2){
    val += __shfl_down(val,offset);
  }
#else
  val += __shfl_down(val,16);
  val += __shfl_down(val, 8);
  val += __shfl_down(val, 4);
  val += __shfl_down(val, 2);
  val += __shfl_down(val, 1);
#endif
  return val;
}

template <class X>
__forceinline__ __device__
X warp_all_reduce(X val){
#if 0
  for(int mask = WarpSize/2; mask > 0; mask /= 2){
    val += __shfl_xor(val,mask);
  }
#else
  val += __shfl_xor(val,16);
  val += __shfl_xor(val, 8);
  val += __shfl_xor(val, 4);
  val += __shfl_xor(val, 2);
  val += __shfl_xor(val, 1);
#endif
  return val;
}

// broadcast first thread value to other threads in warp
template <class X>
__forceinline__ __device__
X warp_broadcast(X val){
#if 0
  for(int offset = 1; offset <= WarpSize/2; offset *= 2){
    val = __shfl_up(val,offset);
  }
#else
  val = __shfl_up(val, 1);
  val = __shfl_up(val, 2);
  val = __shfl_up(val, 4);
  val = __shfl_up(val, 8);
  val = __shfl_up(val,16);
#endif
  return val;
}

// this function returns the reduced value on the first thread only
template <class X>
__forceinline__ __device__
X block_reduce(X val){
  static __shared__ X shared[WarpSize];
  int lane = threadIdx.x % WarpSize;
  int wid  = threadIdx.x / WarpSize;

  val = warp_reduce(val);
  if(lane==0) shared[wid] = val;
  __syncthreads();
  if(wid==0){
    val = (lane < nwarp()) ? shared[lane] : (X)0;
    val = warp_reduce(val);
  }
  return val;
}

// this function returns the reduced value on all the threads
template <class X>
__inline__ __device__
X block_all_reduce(X val){
  __shared__ X shared[WarpSize];
  int lane = threadIdx.x % WarpSize;
  int wid  = threadIdx.x / WarpSize;

  val = warp_reduce(val);
  if(lane==0) shared[wid] = val;
  __syncthreads();
  if(wid==0){
    val = (lane < nwarp()) ? shared[lane] : (X)0;
    val = warp_all_reduce(val);
    shared[lane] = val;
  }
  __syncthreads();
  if(lane==0) val = shared[wid];
  val = warp_broadcast(val);
  return val;
}

// for float4
template <>
__inline__ __device__
float4 block_reduce<float4>(float4 val){
  static __shared__ float4 shared[WarpSize];
  int lane = threadIdx.x % WarpSize;
  int wid  = threadIdx.x / WarpSize;

  val.x = warp_reduce(val.x);
  val.y = warp_reduce(val.y);
  val.z = warp_reduce(val.z);
  val.w = warp_reduce(val.w);
  if(lane==0) shared[wid] = val;
  __syncthreads();
  if(wid==0){
    val = (lane < nwarp()) ? shared[lane] : make_float4(0.f);
    val.x = warp_reduce(val.x);
    val.y = warp_reduce(val.y);
    val.z = warp_reduce(val.z);
    val.w = warp_reduce(val.w);
  }
  return val;
}

template <>
__inline__ __device__
float4 block_all_reduce(float4 val){
  static __shared__ float4 shared[WarpSize];
  int lane = threadIdx.x % WarpSize;
  int wid  = threadIdx.x / WarpSize;

  val.x = warp_reduce(val.x);
  val.y = warp_reduce(val.y);
  val.z = warp_reduce(val.z);
  val.w = warp_reduce(val.w);
  if(lane==0) shared[wid] = val;
  __syncthreads();
  if(wid==0){
    val = (lane < nwarp()) ? shared[lane] : make_float4(0.f);
    val.x = warp_all_reduce(val.x);
    val.y = warp_all_reduce(val.y);
    val.z = warp_all_reduce(val.z);
    val.w = warp_all_reduce(val.w);
    shared[lane] = val;
  }
  __syncthreads();
  if(lane==0) val = shared[wid];
  val.x = warp_broadcast(val.x);
  val.y = warp_broadcast(val.y);
  val.z = warp_broadcast(val.z);
  val.w = warp_broadcast(val.w);
  return val;
}

// for float3
template <>
__inline__ __device__
float3 block_reduce<float3>(float3 val){
  static __shared__ float3 shared[WarpSize];
  int lane = threadIdx.x % WarpSize;
  int wid  = threadIdx.x / WarpSize;

  val.x = warp_reduce(val.x);
  val.y = warp_reduce(val.y);
  val.z = warp_reduce(val.z);
  if(lane==0) shared[wid] = val;
  __syncthreads();
  if(wid==0){
    val = (lane < nwarp()) ? shared[lane] : make_float3(0.f);
    val.x = warp_reduce(val.x);
    val.y = warp_reduce(val.y);
    val.z = warp_reduce(val.z);
  }
  return val;
}

template <>
__inline__ __device__
float3 block_all_reduce(float3 val){
  static __shared__ float3 shared[WarpSize];
  int lane = threadIdx.x % WarpSize;
  int wid  = threadIdx.x / WarpSize;

  val.x = warp_reduce(val.x);
  val.y = warp_reduce(val.y);
  val.z = warp_reduce(val.z);
  if(lane==0) shared[wid] = val;
  __syncthreads();
  if(wid==0){
    val = (lane < nwarp()) ? shared[lane] : make_float3(0.f);
    val.x = warp_all_reduce(val.x);
    val.y = warp_all_reduce(val.y);
    val.z = warp_all_reduce(val.z);
    shared[lane] = val;
  }
  __syncthreads();
  if(lane==0) val = shared[wid];
  val.x = warp_broadcast(val.x);
  val.y = warp_broadcast(val.y);
  val.z = warp_broadcast(val.z);
  return val;
}

// for float2
template <>
__inline__ __device__
float2 block_reduce(float2 val){
  static __shared__ float2 shared[WarpSize];
  int lane = threadIdx.x % WarpSize;
  int wid  = threadIdx.x / WarpSize;

  val.x = warp_reduce(val.x);
  val.y = warp_reduce(val.y);
  if(lane==0) shared[wid] = val;
  __syncthreads();
  if(wid==0){
    val = (lane < nwarp()) ? shared[lane] : make_float2(0.0f);
    val.x = warp_reduce(val.x);
    val.y = warp_reduce(val.y);
  }
  return val;
}

template <>
__inline__ __device__
float2 block_all_reduce(float2 val){
  static __shared__ float2 shared[WarpSize];
  int lane = threadIdx.x % WarpSize;
  int wid  = threadIdx.x / WarpSize;

  val.x = warp_reduce(val.x);
  val.y = warp_reduce(val.y);
  if(lane==0) shared[wid] = val;
  __syncthreads();
  if(wid==0){
    val = (lane < nwarp()) ? shared[lane] : make_float2(0.0f);
    val.x = warp_all_reduce(val.x);
    val.y = warp_all_reduce(val.y);
    shared[lane] = val;
  }
  __syncthreads();
  if(lane==0) val = shared[wid];
  val.x = warp_broadcast(val.x);
  val.y = warp_broadcast(val.y);
  return val;
}
#endif // end of CUH_REDUCE
