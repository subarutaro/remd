#define MAX_THREAD  8
#ifdef PTHREADS
#define FLG_PTHREAD     1
#else // PTHREADS
#define FLG_PTHREAD     0
#endif // PTHREADS
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>
#if FLG_PTHREAD
#include <pthread.h>
#endif // FLG_PTHREAD

#define EWALDCOEF1 1.1283792f
// CERFC1~6 erfc(x)
#define CERFC1 0.0705230784f
#define CERFC2 0.0422820123f
#define CERFC3 0.0092705272f
#define CERFC4 0.0001520143f
#define CERFC5 0.0002765672f
#define CERFC6 0.0000430638f
// CEXP1~7 for 2/sqrt(pi)*exp(-x)
#define CEXP1 1.6211389f
// 1.621138938277405f
#define CEXP2 0.40528260f
// 0.405282601474736f
#define CEXP3 0.050672885f
// 0.05067288524197f
#define CEXP4 0.0042009728f
// 0.004200972755851f
#define CEXP5 0.00027812584f
// 0.0002781258385287f
#define CEXP6 0.0000088031087f
// 0.000008803108662634f
#define CEXP7 0.0000011195586f
// 0.000001119558550774f

#define VG_PME_WITHOUT_DIFF 1
#define VG_PME_WITH_DIFF    2
#define VG_FFT_FORWARD      1.e0
#define VG_FFT_BACKWARD     -1.e0
#define VG_COUNT_DOUBLE     1
#define VG_COUNT_SINGLE     0
#define VG_EXCLUSION_REAL   0
#define VG_EXCLUSION_BOTH   1
#define VG_RANDSEED_MAN     0
#define VG_RANDSEED_AUTO    1

#define STOPWITHMSG(MSG)						\
  { printf("VGLIBRARY ERROR : %s (FILE: %s, LINE %d)\n",MSG,__FILE__,__LINE__);	\
    fprintf(stderr,"VGLIBRARY ERROR : %s (FILE: %s, LINE %d)\n",MSG,__FILE__,__LINE__); \
    exit(EXIT_FAILURE); }

#define MEMFREE(ptr)            { free(ptr); ptr = NULL; }
#define MEMALLOC(ptr,size,type) { ptr = (type*)malloc(size); }

#if 1
#define MEMALIGN(err,ptr,bnd,size,type) { ptr = (type*)malloc(size); }
#else // 
#define MEMALIGN(err,ptr,bnd,size,type)	{ err = posix_memalign(ptr,bnd,size);}
#endif // 

#define SHIFTCELLLISTX 0
#define SHIFTCELLLISTY 4
#define SHIFTCELLLISTZ 8
#define MASKCELLLISTX 0xf
#define MASKCELLLISTY 0xf0
#define MASKCELLLISTZ 0xf00

#define SETNEIGHBORCELLINDEX(ii,jj,kk,SUBDIV) \
  ((((ii)+SUBDIV)<<SHIFTCELLLISTX)+(((jj)+SUBDIV)<<SHIFTCELLLISTY)+(((kk)+SUBDIV)<<SHIFTCELLLISTZ))
#define GETNEIGHBORCELLX(INDEX) (INDEX&MASKCELLLISTX)
#define GETNEIGHBORCELLY(INDEX) ((INDEX&MASKCELLLISTY)>>SHIFTCELLLISTY)
#define GETNEIGHBORCELLZ(INDEX) ((INDEX&MASKCELLLISTZ)>>SHIFTCELLLISTZ)
