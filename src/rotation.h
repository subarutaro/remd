#ifndef H_ROTATION
#define H_ROTATION

#include "vector.h"
#include "timer.h"

static inline
dvec3 space_to_body(const dvec4 a,const dvec3 v){
  const double a0sq = a[0]*a[0],a1sq = a[1]*a[1],a2sq = a[2]*a[2],a3sq = a[3]*a[3];
  const double a01 = a[0]*a[1],a02 = a[0]*a[2],a03 = a[0]*a[3];
  const double a12 = a[1]*a[2],a13 = a[1]*a[3];
  const double a23 = a[2]*a[3];
  dvec3 tmp;
  tmp[0] = (a0sq + a1sq - a2sq - a3sq)*v[0];
  tmp[0] +=  2.*(a12 + a03)*v[1];
  tmp[0] +=  2.*(a13 - a02)*v[2];

  tmp[1] =  2.*(a12 - a03)*v[0];
  tmp[1] += (a0sq - a1sq + a2sq - a3sq)*v[1];
  tmp[1] +=  2.*(a23 + a01)*v[2];

  tmp[2] =  2.*(a13 + a02)*v[0];
  tmp[2] += 2.*(a23 - a01)*v[1];
  tmp[2] += (a0sq - a1sq - a2sq + a3sq)*v[2];

  return tmp;
}

static inline
dvec3 body_to_space(const dvec4 a,const dvec3 v){
  const double a0sq = a[0]*a[0],a1sq = a[1]*a[1],a2sq = a[2]*a[2],a3sq = a[3]*a[3];
  const double a01 = a[0]*a[1],a02 = a[0]*a[2],a03 = a[0]*a[3];
  const double a12 = a[1]*a[2],a13 = a[1]*a[3];
  const double a23 = a[2]*a[3];
  dvec3 tmp;
  tmp[0] = (a0sq + a1sq - a2sq - a3sq)*v[0];
  tmp[0] +=  2.*(a12 - a03)*v[1];
  tmp[0] +=  2.*(a13 + a02)*v[2];

  tmp[1] =  2.*(a12 + a03)*v[0];
  tmp[1] += (a0sq - a1sq + a2sq - a3sq)*v[1];
  tmp[1] +=  2.*(a23 - a01)*v[2];

  tmp[2] =  2.*(a13 - a02)*v[0];
  tmp[2] += 2.*(a23 + a01)*v[1];
  tmp[2] += (a0sq - a1sq - a2sq + a3sq)*v[2];
  return tmp;
}
#endif
