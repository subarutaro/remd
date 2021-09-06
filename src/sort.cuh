#ifndef CUH_SORT
#define CUH_SORT
__global__
void sort(int4   *r_dev,float4 *v_dev,float4 *f_dev,
	  qtype4 *q_dev,ptype4 *p_dev,float4 *n_dev,int nelem);
#endif
