//#include <cutil.h>
#include <cufft.h>

typedef struct __align__(16){
  float x;
  float y;
  float z;
  unsigned int w;
} VGGVDWATM;

typedef struct __align__(16){
  float x;
  float y;
  float z;
  unsigned int w;
} VGGGBATM1;

typedef struct __align__(16){
  float ex;
  float ey;
  float ez;
  unsigned int w;
} VGGGBATM2;

typedef struct __align__(32){
  float x;
  float y;
  float z;
  float ex;
  float ey;
  float ez;
  unsigned int w;
  unsigned int ww;
} VGGGBATM;

typedef struct __align__(16){
  int x;
  int y;
  int z;
  unsigned int w;
} VGGVDWATMINT;

typedef struct __align__(16){
  int x;
  int y;
  int z;
  float w;
} VGGATMINT;

typedef union {
  struct {
    int i0;
    int i1;
  } i2;
  double d;
} VGGDI2;

typedef struct __align__(8){
  float hs;
  int   li;
} VGGSI;

typedef struct{
  double* coor_i;
  double* coor_j;
  int*    kvec_j;
  double* charge_i;
  double* charge_j;
  int   num_atm_i;
  int   num_atm_j;
  int   num_kvecj;
  double* force;
  double* pot;
  double* stress;
  double  alpha;
  int     idevice;
  double* recip;
} VGGTHREAD;

#include <vg_common.h>
#include <vg_util.h>
#include <vgg_util.h>

#define VGG_FFT_FORWARD      1.e0
#define VGG_FFT_BACKWARD     -1.e0
#define FLG_VGGNOTCHARMM     0.e0

#ifdef TIMER
#define FLG_VGGTIMER         1
#else // TIMER
#define FLG_VGGTIMER         0
#endif // TIMER
#ifdef CUDA2
#define FLG_CUDA2            1
#else // CUDA2
#define FLG_CUDA2            0
#endif // CUDA2
#define ERRORCHECKER         0
#if ERRORCHECKER && FLG_CUDA2
//#define CHECK_ERROR(x) CUDAError(x,1);
#define VGGERR(x) printf("%d %s\n",x,cudaGetErrorString(cudaGetLastError()));
#else // ERRORCHECKER
#define VGGERR(x)
#endif // ERRORCHECKER

#if FLG_CUDA2
#define GPUMUL(x,y) (__fmul_rn((x),(y)))
#else // FLG_CUDA2
#define GPUMUL(x,y) ((x)*(y))
#endif // FLG_CUDA2

#if FLG_VGGTIMER
#define VGG_TIMER_START(X1,X2)          get_cputime(&X1,&X2);
#define VGG_TIMER_STOP_HOST(X1,X2,X3)   get_cputime(&X1,&X2); X3 += X2;
#define VGG_TIMER_STOP_DEVICE(X1,X2,X3) cudaThreadSynchronize(); get_cputime(&X1,&X2); X3 += X2;
#else // FLG_VGGTIMER
#define VGG_TIMER_START(X1,X2)
#define VGG_TIMER_STOP_HOST(X1,X2,X3)
#define VGG_TIMER_STOP_DEVICE(X1,X2,X3)
#endif // FLG_VGGTIMER

#define SET_FLOAT4_ZERO(f4)				\
  f4.x = 0.f; f4.y = 0.f; f4.z = 0.f; f4.w = 0.f;
#define SET_FLOAT2_ZERO(f2)				\
  f2.x = 0.f; f2.y = 0.f;

#define LARGE_SHIFT    15
#define LARGE          (float)(3 << (LARGE_SHIFT-1))
#define LOWER_SHIFT    7
#define LOWER_LOOP     (1 << LOWER_SHIFT)
#define LOWER_FACTOR   (float)(1LL << (23-LARGE_SHIFT+32-LOWER_SHIFT))
#define LOWER_FACTOR_1 (1.f / LOWER_FACTOR)
#define MASK(n)        ((0x1<<(n)) -1)

#define INIT_NARUMI(si)   si.hs = LARGE; si.li = 0;
#define COPY_NARUMI(si)   (si.hs - LARGE + si.li * (double)LOWER_FACTOR_1)
#define UPDATE_NARUMI(si)						\
  si.hs += GPUMUL((float)(si.li & (MASK(LOWER_SHIFT)<<(32-LOWER_SHIFT))),LOWER_FACTOR_1); \
  si.li = si.li & MASK(32-LOWER_SHIFT);
#define ADD_NARUMI(si,ys,tmp1,tmp2)		\
  tmp1 = si.hs + (ys);				\
  tmp2 = tmp1 - si.hs;				\
  si.li += (int)(GPUMUL(((ys)-tmp2),LOWER_FACTOR));	\
  si.hs = tmp1;
#define ADDSI_NARUMI(out,in,ys,tmp1,tmp2)	\
  ys = in.hs - LARGE;				\
  ADD_NARUMI(out,ys,tmp1,tmp2);			\
  ys = in.li * (float)LOWER_FACTOR_1;		\
  ADD_NARUMI(out,ys,tmp1,tmp2);

#if 0
#define FREE_HOST(ptr) {cudaFreeHost(ptr); ptr = NULL;}
#define ALLOC_HOST(ptr,size,type) {cudaMallocHost((void**)&ptr,size);}
#else
#define FREE_HOST(ptr) {free(ptr); ptr = NULL;}
#define ALLOC_HOST(ptr,size,type) {ptr = (type*)malloc(size);}
#endif

// for only PME
#if 0
#define FREE_PME_HOST(ptr) {cudaFreeHost(ptr); ptr = NULL;}
#define ALLOC_PME_HOST(PTR,TYPE,SIZE) {cudaMallocHost((void**)&PTR,(SIZE));}
#else
#define FREE_PME_HOST(ptr) {free(ptr); ptr = NULL;}
#define ALLOC_PME_HOST(PTR,TYPE,SIZE) {PTR = (TYPE*)malloc((SIZE));}
#endif

