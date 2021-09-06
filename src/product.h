#include <cmath>
#include "vector.h"

//===scalar prod=======
template <int N>
static inline
double scalar_prod
(const Vector<N,double> v1,const Vector<N,double> v2){
  double tmp = 0.0;
  for(int d=0;d<N;d++) tmp += v1[d]*v2[d];
  return tmp;
}

template <int N>
static inline
Vector<N,double>  scalar_prod
(const Vector< N,Vector<N,double> > v1,const Vector<N,double> v2){
  Vector<N,double> tmp;
  for(int d=0;d<N;d++){
    tmp[d] = scalar_prod(v1[d],v2);
  }
  return tmp;
}

template <int N>
static inline
Vector< N,Vector<N,double> > scalar_prod
(const Vector< N,Vector<N,double> > v1,const Vector< N,Vector<N,double> > v2){
  const Vector< N,Vector<N,double> > v2i = transposed(v2);
  Vector< N,Vector<N,double> > tmp;
  for(int i=0;i<N;i++){
  for(int j=0;j<N;j++){
    tmp[i][j] = scalar_prod(v1[i],v2i[j]);
  }}
  return tmp;
}
//==========

//===vector_prod=======
static inline 
dvec3 vector_prod(const dvec3 v1,const dvec3 v2){
  dvec3 tmp;
  tmp[0] = v1[1]*v2[2] - v1[2]*v2[1];
  tmp[1] = v1[2]*v2[0] - v1[0]*v2[2];
  tmp[2] = v1[0]*v2[1] - v1[1]*v2[0];
  return tmp;
}
//==========

//===dyrac prod=======
static inline
dvec33 dyadic_prod(const dvec3 v1,const dvec3 v2){
  dvec33 tmp;
  for(int i=0;i<3;i++){
  for(int j=0;j<3;j++){
    tmp[i][j] = v1[i]*v2[j];
  }}
  return tmp;
}
//==========

static inline
dvec3 trace(const dvec33 v){
  dvec3 tmp;
  for(int d=0;d<3;d++) tmp[d] = v[d][d];
  return tmp;
}
//==========
