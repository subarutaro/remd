#ifndef BITONICSORT_CUH
#define BITONICSORT_CUH

#ifndef WarpSize
#define WarpSize 32
#endif


template <int N,class X>
__global__
void bitonic_sort_key(X*,X*);


#endif
