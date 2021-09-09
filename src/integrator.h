#ifndef H_INTEGRATOR
#define H_INTEGRATOR

#include "vector.h"

#define E_CONSTANT 0
#define T_CONSTANT 1
#define T_SCALE 2

#define V_CONSTANT 0
#define P_CONSTANT 1
#define PxAyz_CONSTANT 2
#define PxyLz_CONSTANT 3

const int NBIT      = 2;
const int MASK      = ((1<<NBIT) - 1);
const int TSHIFT    = 0;
const int PSHIFT    = (NBIT);
const int CSHIFT    = (2*NBIT);

const int NVE            = ((0<<CSHIFT) + (0<<PSHIFT) + (0<<TSHIFT));
const int NVT            = ((0<<CSHIFT) + (0<<PSHIFT) + (1<<TSHIFT));
const int NVTSCALE       = ((0<<CSHIFT) + (0<<PSHIFT) + (2<<TSHIFT));
const int NPH            = ((0<<CSHIFT) + (1<<PSHIFT) + (0<<TSHIFT));
const int NPT            = ((0<<CSHIFT) + (1<<PSHIFT) + (1<<TSHIFT));
const int NPTSCALE       = ((0<<CSHIFT) + (1<<PSHIFT) + (2<<TSHIFT));
const int NAxyPzH        = ((0<<CSHIFT) + (2<<PSHIFT) + (0<<TSHIFT));
const int NAxyPzT        = ((0<<CSHIFT) + (2<<PSHIFT) + (1<<TSHIFT));
const int NAxyPzTSCALE   = ((0<<CSHIFT) + (2<<PSHIFT) + (2<<TSHIFT));
const int NPxyLzH        = ((0<<CSHIFT) + (3<<PSHIFT) + (0<<TSHIFT));
const int NPxyLzT        = ((0<<CSHIFT) + (3<<PSHIFT) + (1<<TSHIFT));
const int NPxyLzTSCALE   = ((0<<CSHIFT) + (3<<PSHIFT) + (2<<TSHIFT));
const int NVE1D          = ((1<<CSHIFT) + (0<<PSHIFT) + (0<<TSHIFT));
const int NVT1D          = ((1<<CSHIFT) + (0<<PSHIFT) + (1<<TSHIFT));
const int NVTSCALE1D     = ((1<<CSHIFT) + (0<<PSHIFT) + (2<<TSHIFT));
const int NPH1D          = ((1<<CSHIFT) + (1<<PSHIFT) + (0<<TSHIFT));
const int NPT1D          = ((1<<CSHIFT) + (1<<PSHIFT) + (1<<TSHIFT));
const int NPTSCALE1D     = ((1<<CSHIFT) + (1<<PSHIFT) + (2<<TSHIFT));
const int NAxyPzH1D      = ((1<<CSHIFT) + (2<<PSHIFT) + (0<<TSHIFT));
const int NAxyPzT1D      = ((1<<CSHIFT) + (2<<PSHIFT) + (1<<TSHIFT));
const int NAxyPzTSCALE1D = ((1<<CSHIFT) + (2<<PSHIFT) + (2<<TSHIFT));
const int NPxyLzH1D      = ((1<<CSHIFT) + (3<<PSHIFT) + (0<<TSHIFT));
const int NPxyLzT1D      = ((1<<CSHIFT) + (3<<PSHIFT) + (1<<TSHIFT));
const int NPxyLzTSCALE1D = ((1<<CSHIFT) + (3<<PSHIFT) + (2<<TSHIFT));
const int NVE2D          = ((2<<CSHIFT) + (0<<PSHIFT) + (0<<TSHIFT));
const int NVT2D          = ((2<<CSHIFT) + (0<<PSHIFT) + (1<<TSHIFT));
const int NVTSCALE2D     = ((2<<CSHIFT) + (0<<PSHIFT) + (2<<TSHIFT));
const int NPH2D          = ((2<<CSHIFT) + (1<<PSHIFT) + (0<<TSHIFT));
const int NPT2D          = ((2<<CSHIFT) + (1<<PSHIFT) + (1<<TSHIFT));
const int NPTSCALE2D     = ((2<<CSHIFT) + (1<<PSHIFT) + (2<<TSHIFT));
const int NAxyPzH2D      = ((2<<CSHIFT) + (2<<PSHIFT) + (0<<TSHIFT));
const int NAxyPzT2D      = ((2<<CSHIFT) + (2<<PSHIFT) + (1<<TSHIFT));
const int NAxyPzTSCALE2D = ((2<<CSHIFT) + (2<<PSHIFT) + (2<<TSHIFT));
const int NPxyLzH2D      = ((2<<CSHIFT) + (3<<PSHIFT) + (0<<TSHIFT));
const int NPxyLzT2D      = ((2<<CSHIFT) + (3<<PSHIFT) + (1<<TSHIFT));
const int NPxyLzTSCALE2D = ((2<<CSHIFT) + (3<<PSHIFT) + (2<<TSHIFT));

#define TMODE ((MODE>>TSHIFT)&MASK)
#define PMODE ((MODE>>PSHIFT)&MASK)
#define CMODE ((MODE>>CSHIFT)&MASK)

template <int N,int M>
static inline
Vector< M,Vector<N,double> > transposed
(Vector< N,Vector<M,double> > v){
  Vector<M,Vector<N,double> > tmp;
  for(int i=0;i<M;i++)
  for(int j=0;j<N;j++)
    tmp[i][j] = v[j][i];
  return tmp;
}

static inline
dvec33 inverse(const dvec33 v){
  const double det
    = v[0][0]*v[1][1]*v[2][2] + v[1][0]*v[2][1]*v[0][2] + v[2][0]*v[0][1]*v[1][2]
    - v[0][0]*v[2][1]*v[1][2] - v[2][0]*v[1][1]*v[0][2] - v[1][0]*v[0][1]*v[2][2];
  if(det<1e-10) assert(det<1e-10);
  dvec33 tmp;
  tmp[0][0] = v[1][1]*v[2][2] - v[1][2]*v[2][1];
  tmp[0][1] = v[0][2]*v[2][1] - v[0][1]*v[2][2];
  tmp[0][2] = v[0][1]*v[1][2] - v[0][2]*v[1][1];

  tmp[1][0] = v[1][2]*v[2][0] - v[1][0]*v[2][2];
  tmp[1][1] = v[0][0]*v[2][2] - v[0][2]*v[2][0];
  tmp[1][2] = v[0][2]*v[1][0] - v[0][0]*v[1][2];

  tmp[2][0] = v[1][0]*v[2][1] - v[1][1]*v[2][0];
  tmp[2][1] = v[0][1]*v[2][0] - v[0][0]*v[2][1];
  tmp[2][2] = v[0][0]*v[1][1] - v[0][1]*v[1][0];
  return tmp/det;
}

static inline dvec44 S(const dvec4 q){
  dvec44 tmp;
  tmp[0][0] = q[0];tmp[0][1] =-q[1];tmp[0][2] =-q[2];tmp[0][3] =-q[3];
  tmp[1][0] = q[1];tmp[1][1] = q[0];tmp[1][2] =-q[3];tmp[1][3] = q[2];
  tmp[2][0] = q[2];tmp[2][1] = q[3];tmp[2][2] = q[0];tmp[2][3] =-q[1];
  tmp[3][0] = q[3];tmp[3][1] =-q[2];tmp[3][2] = q[1];tmp[3][3] = q[0];
  return tmp;
}

static inline dvec44 S_T(const dvec4 q){
  dvec44 tmp;
  tmp[0][0] = q[0];tmp[0][1] = q[1];tmp[0][2] = q[2];tmp[0][3] = q[3];
  tmp[1][0] =-q[1];tmp[1][1] = q[0];tmp[1][2] = q[3];tmp[1][3] =-q[2];
  tmp[2][0] =-q[2];tmp[2][1] =-q[3];tmp[2][2] = q[0];tmp[2][3] = q[1];
  tmp[3][0] =-q[3];tmp[3][1] = q[2];tmp[3][2] =-q[1];tmp[3][3] = q[0];
  return tmp;
}

static inline dvec4 P1(const dvec4 q){
  dvec4 tmp;
  tmp[0] =-q[1];
  tmp[1] = q[0];
  tmp[2] = q[3];
  tmp[3] =-q[2];
  return tmp;
}

static inline dvec4 P2(const dvec4 q){
  dvec4 tmp;
  tmp[0] =-q[2];
  tmp[1] =-q[3];
  tmp[2] = q[0];
  tmp[3] = q[1];
  return tmp;
}

static inline dvec4 P3(const dvec4 q){
  dvec4 tmp;
  tmp[0] =-q[3];
  tmp[1] = q[2];
  tmp[2] =-q[1];
  tmp[3] = q[0];
  return tmp;
}

#endif
