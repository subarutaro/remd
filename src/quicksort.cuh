#ifndef H_QUICKSORT
#define H_QUICKSORT
template <class X>
__global__ void quicksort(X *gdata,int left,int right,int depth = 0);

template <class X0,class X1>
__global__ void quicksort_key(X0 *keys,X1 *values,int left,int right,int depth = 0);
#endif
