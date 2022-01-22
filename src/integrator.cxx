#include "molecules.h"
#include "integrator.h"
#ifdef ENABLE_AOS_TO_SOA_CONVERSION
void Molecules::AoStoSoA(){
#ifdef _OPENMP
  const int thread = omp_get_thread_num();
#else
  const int thread = 0;
#endif

#pragma omp barrier

  assert(isAoS == true);
  //#pragma omp parallel for
  for(int i=is[thread];i<ie[thread];i++){
    const Molecule& m = mlcl[i];
    mi[i] = 1.0 / m.m;
    rx[i] = m.r[0];
    ry[i] = m.r[1];
    rz[i] = m.r[2];

    vx[i] = m.v[0];
    vy[i] = m.v[1];
    vz[i] = m.v[2];

    fx[i] = m.f[0];
    fy[i] = m.f[1];
    fz[i] = m.f[2];

    qx[i] = m.q[0];
    qy[i] = m.q[1];
    qz[i] = m.q[2];
    qw[i] = m.q[3];

    px[i] = m.p[0];
    py[i] = m.p[1];
    pz[i] = m.p[2];
    pw[i] = m.p[3];

    nx[i] = m.n[0];
    ny[i] = m.n[1];
    nz[i] = m.n[2];
    nw[i] = m.n[3];

    ix[i] = m.i[0];
    iy[i] = m.i[1];
    iz[i] = m.i[2];
  }
  #pragma omp barrier
  #pragma omp single
  {
    isAoS = false;
  }
}
void Molecules::SoAtoAoS(){
#ifdef _OPENMP
  const int thread = omp_get_thread_num();
#else
  const int thread = 0;
#endif

#pragma omp barrier

  assert(isAoS == false);
  //#pragma omp parallel for
  for(int i=is[thread];i<ie[thread];i++){
    Molecule& m = mlcl[i];
#if 0
    if(m.r[0] != rx[i] || m.r[1] != ry[i] || m.r[2] != rz[i]){
      std::cout << "AoS: " << i << " " << m.r[0] << " " << m.r[1] << " " << m.r[2] << std::endl;
      std::cout << "SoA: " << i << " " << rx[i]  << " " <<  ry[i] << " " <<  rz[i] << std::endl;
    }
    assert(m.r[0] == rx[i]);
    assert(m.r[1] == ry[i]);
    assert(m.r[2] == rz[i]);
#endif
    m.m = 1.0 / mi[i];

    m.r[0] = rx[i];
    m.r[1] = ry[i];
    m.r[2] = rz[i];

    m.q[0] = qx[i];
    m.q[1] = qy[i];
    m.q[2] = qz[i];
    m.q[3] = qw[i];

    m.p[0] = px[i];
    m.p[1] = py[i];
    m.p[2] = pz[i];
    m.p[3] = pw[i];

    m.v[0] = vx[i];
    m.v[1] = vy[i];
    m.v[2] = vz[i];

    m.f[0] = fx[i];
    m.f[1] = fy[i];
    m.f[2] = fz[i];

    m.n[0] = nx[i];
    m.n[1] = ny[i];
    m.n[2] = nz[i];
    m.n[3] = nw[i];

    m.i[0] = ix[i];
    m.i[1] = iy[i];
    m.i[2] = iz[i];
  }
#pragma omp barrier
#pragma omp single
  {
    isAoS = true;
  }
}
#endif

template <int MODE>
inline void Molecules::init_D3D2D1D2D3_PsPv(){
#ifdef _OPENMP
  const int thread = omp_get_thread_num();
#else
  const int thread = 0;
#endif
  for(int i=0;i<nlane;i++){
    sum_Pv_omp[nlane*thread+i] = 0.0;
    sum_Ps_omp[nlane*thread+i] = 0.0;
  }
}
template <int MODE>
inline void Molecules::D3D2D1D2D3_PsPv(){
#ifdef _OPENMP
  const int thread = omp_get_thread_num();
#else
  const int thread = 0;
#endif
  for(int i=1;i<nlane;i++){
    sum_Ps_omp[nlane*thread] += sum_Ps_omp[nlane*thread+i];
    sum_Pv_omp[nlane*thread] += sum_Pv_omp[nlane*thread+i];
  }
  #pragma omp barrier
  #pragma omp single
  {
    dvec3  sum_Pv  = 0.0;
    if(((MODE>>PSHIFT)&MASK)>0 || ((MODE>>TSHIFT)&MASK)==T_CONSTANT){
      for(int i=0;i<nthreads;i++) sum_Pv  += sum_Pv_omp[nlane*i];
      sum_Pv[0] /= L[0]*L[0] * tst->s*tst->s;
      sum_Pv[1] /= L[1]*L[1] * tst->s*tst->s;
      sum_Pv[2] /= L[2]*L[2] * tst->s*tst->s;
    }
    if(((MODE>>PSHIFT)&MASK)==1){
      bst->Pv[0] += tst->s * sum(sum_Pv)/(Molecules::GetVolume()*3.0) * dthalf;
    }
    if(((MODE>>PSHIFT)&MASK)==2){
      bst->Pv[2] += tst->s * sum_Pv[2] / L[2] * dthalf;
      //printf("(kin)= %lff\n",sum_Pv[2] / L[2]);
    }
    if(((MODE>>PSHIFT)&MASK)==3){
      bst->Pv[0] += tst->s * sum_Pv[0] / L[0] * dthalf;
      bst->Pv[1] += tst->s * sum_Pv[1] / L[1] * dthalf;
    }

    if(((MODE>>TSHIFT)&MASK)==T_CONSTANT){
      //std::cout << "D2 sum_Ps: " << 2.0*sum_Ps << " s: " << tst->s << ", Ps: " << tst->Ps << std::endl;
      double sum_Ps = 0.0;
      for(int i=0;i<nthreads;i++) sum_Ps += sum_Ps_omp[nlane*i];
      tst->Ps += ( 0.25*sum(sum_Pv) + 2.0*sum_Ps - ( 1.0+log(tst->s) )*prop.gkT + prop.H0 ) * dt;
    }
  }
}

#ifdef ENABLE_AOS_TO_SOA_CONVERSION
template <int MODE>
inline void Molecules::D1(){
#ifdef _OPENMP
  const int thread = omp_get_thread_num();
#else
  const int thread = 0;
#endif

  const double si   = 1.0/tst->s;
  const dvec3  Li   = L.inv();
  const dvec3  coef_r  = Li*Li*si*dt;
  const double coef_xi = 0.25*si;
  //const double coef_tra = si*si*0.5;
  //double sum_tra = 0.0;
  const int ivs = is[thread];
  const int ive = is[thread] + ((ie[thread] - is[thread])/nlane)*nlane;
  dvec3*  sum_Pv = &sum_Pv_omp[nlane*thread];
  double* sum_Ps = &sum_Ps_omp[nlane*thread];
  for(int iv=ivs;iv<ive;iv+=nlane){
    #pragma omp simd
    for(int ii=0;ii<nlane;ii++){
      const int i = iv + ii;
      //update coordinate
      //update coordinate
      rx[i] += vx[i]*coef_r[0];
      ry[i] += vy[i]*coef_r[1];
      rz[i] += vz[i]*coef_r[2];
      if(rx[i] <  0.0) rx[i] += 1.0;
      if(rx[i] >= 1.0) rx[i] -= 1.0;
      if(ry[i] <  0.0) ry[i] += 1.0;
      if(ry[i] >= 1.0) ry[i] -= 1.0;
      if(rz[i] <  0.0) rz[i] += 1.0;
      if(rz[i] >= 1.0) rz[i] -= 1.0;
      dvec4  Pq; Pq[0] = -qy[i]; Pq[1] = qx[i]; Pq[2] = qw[i]; Pq[3] = -qz[i];
      dvec4  Pp; Pp[0] = -py[i]; Pp[1] = px[i]; Pp[2] = pw[i]; Pp[3] = -pz[i];
      const double xi = (px[i]*Pq[0] + py[i]*Pq[1] + pz[i]*Pq[2] + pw[i]*Pq[3])  * coef_xi/ix[i];
      const double xidt = xi*dt;
      const double cosxidt = cos(xidt);
      const double sinxidt = sin(xidt);
      qx[i] = qx[i]*cosxidt + Pq[0]*sinxidt;
      qy[i] = qy[i]*cosxidt + Pq[1]*sinxidt;
      qz[i] = qz[i]*cosxidt + Pq[2]*sinxidt;
      qw[i] = qw[i]*cosxidt + Pq[3]*sinxidt;
      px[i] = px[i]*cosxidt + Pp[0]*sinxidt;
      py[i] = py[i]*cosxidt + Pp[1]*sinxidt;
      pz[i] = pz[i]*cosxidt + Pp[2]*sinxidt;
      pw[i] = pw[i]*cosxidt + Pp[3]*sinxidt;

      if constexpr (((MODE>>PSHIFT)&MASK)>0 || ((MODE>>TSHIFT)&MASK)==T_CONSTANT ){
	//sum_Pv += m.v*m.v*Li*Li*m.m;
	const double m = 1.0 / mi[i];
	sum_Pv[ii][0] += 2.0 * m * vx[i]*vx[i];
	sum_Pv[ii][1] += 2.0 * m * vy[i]*vy[i];
	sum_Pv[ii][2] += 2.0 * m * vz[i]*vz[i];
      }
      if constexpr (((MODE>>TSHIFT)&MASK)==T_CONSTANT){
	sum_Ps[ii] += xi*xi*ix[i];
      }
    }
  }
  for(int i=ive;i<ie[thread];i++){
    //update coordinate
    //update coordinate
    rx[i] += vx[i]*coef_r[0];
    ry[i] += vy[i]*coef_r[1];
    rz[i] += vz[i]*coef_r[2];
    if(rx[i] <  0.0) rx[i] += 1.0;
    if(rx[i] >= 1.0) rx[i] -= 1.0;
    if(ry[i] <  0.0) ry[i] += 1.0;
    if(ry[i] >= 1.0) ry[i] -= 1.0;
    if(rz[i] <  0.0) rz[i] += 1.0;
    if(rz[i] >= 1.0) rz[i] -= 1.0;
    dvec4  Pq; Pq[0] = -qy[i]; Pq[1] = qx[i]; Pq[2] = qw[i]; Pq[3] = -qz[i];
    dvec4  Pp; Pp[0] = -py[i]; Pp[1] = px[i]; Pp[2] = pw[i]; Pp[3] = -pz[i];
    const double xi = (px[i]*Pq[0] + py[i]*Pq[1] + pz[i]*Pq[2] + pw[i]*Pq[3])  * coef_xi/ix[i];
    const double xidt = xi*dt;
    const double cosxidt = cos(xidt);
    const double sinxidt = sin(xidt);
    qx[i] = qx[i]*cosxidt + Pq[0]*sinxidt;
    qy[i] = qy[i]*cosxidt + Pq[1]*sinxidt;
    qz[i] = qz[i]*cosxidt + Pq[2]*sinxidt;
    qw[i] = qw[i]*cosxidt + Pq[3]*sinxidt;
    px[i] = px[i]*cosxidt + Pp[0]*sinxidt;
    py[i] = py[i]*cosxidt + Pp[1]*sinxidt;
    pz[i] = pz[i]*cosxidt + Pp[2]*sinxidt;
    pw[i] = pw[i]*cosxidt + Pp[3]*sinxidt;

    if constexpr (((MODE>>PSHIFT)&MASK)>0 || ((MODE>>TSHIFT)&MASK)==T_CONSTANT ){
      //sum_Pv += m.v*m.v*Li*Li*m.m;
      const double m = 1.0 / mi[i];
      const int ii = i%nlane;
      sum_Pv_omp[nlane*thread+ii][0] += 2.0 * m * vx[i]*vx[i];
      sum_Pv_omp[nlane*thread+ii][1] += 2.0 * m * vy[i]*vy[i];
      sum_Pv_omp[nlane*thread+ii][2] += 2.0 * m * vz[i]*vz[i];
    }
    if constexpr (((MODE>>TSHIFT)&MASK)==T_CONSTANT){
      const int ii = i%nlane;
      sum_Ps_omp[nlane*thread+ii] += xi*xi*ix[i];
    }
  }
}

template <int MODE>
inline void Molecules::D2(){
#ifdef _OPENMP
  const int thread = omp_get_thread_num();
#else
  const int thread = 0;
#endif

  const double coef_xi = 0.25 / tst->s;

  const int ivs = is[thread];
  const int ive = is[thread] + ((ie[thread] - is[thread])/nlane)*nlane;
  double* sum_Ps = &sum_Ps_omp[nlane*thread];
  for(int iv=ivs;iv<ive;iv+=nlane){
    #pragma omp simd
    for(int ii=0;ii<nlane;ii++){
      const int i = iv + ii;
      dvec4  Pq; Pq[0] = -qz[i]; Pq[1] = -qw[i]; Pq[2] = qx[i]; Pq[3] = qy[i];
      dvec4  Pp; Pp[0] = -pz[i]; Pp[1] = -pw[i]; Pp[2] = px[i]; Pp[3] = py[i];
      const double xi = (px[i]*Pq[0] + py[i]*Pq[1] + pz[i]*Pq[2] + pw[i]*Pq[3])*coef_xi/iy[i];
      const double xidt = xi*dthalf;
      const double cosxidt = cos(xidt);
      const double sinxidt = sin(xidt);
      qx[i] = qx[i]*cosxidt + Pq[0]*sinxidt;
      qy[i] = qy[i]*cosxidt + Pq[1]*sinxidt;
      qz[i] = qz[i]*cosxidt + Pq[2]*sinxidt;
      qw[i] = qw[i]*cosxidt + Pq[3]*sinxidt;

      px[i] = px[i]*cosxidt + Pp[0]*sinxidt;
      py[i] = py[i]*cosxidt + Pp[1]*sinxidt;
      pz[i] = pz[i]*cosxidt + Pp[2]*sinxidt;
      pw[i] = pw[i]*cosxidt + Pp[3]*sinxidt;

      if(((MODE>>TSHIFT)&MASK)==T_CONSTANT){
	sum_Ps_omp[nlane*thread+ii] += 0.5*xi*xi*iy[i];
      }
    }
  }
  for(int i=ive;i<ie[thread];i++){
    dvec4  Pq; Pq[0] = -qz[i]; Pq[1] = -qw[i]; Pq[2] = qx[i]; Pq[3] = qy[i];
    dvec4  Pp; Pp[0] = -pz[i]; Pp[1] = -pw[i]; Pp[2] = px[i]; Pp[3] = py[i];
    const double xi = (px[i]*Pq[0] + py[i]*Pq[1] + pz[i]*Pq[2] + pw[i]*Pq[3])*coef_xi/iy[i];
    const double xidt = xi*dthalf;
    const double cosxidt = cos(xidt);
    const double sinxidt = sin(xidt);
    qx[i] = qx[i]*cosxidt + Pq[0]*sinxidt;
    qy[i] = qy[i]*cosxidt + Pq[1]*sinxidt;
    qz[i] = qz[i]*cosxidt + Pq[2]*sinxidt;
    qw[i] = qw[i]*cosxidt + Pq[3]*sinxidt;

    px[i] = px[i]*cosxidt + Pp[0]*sinxidt;
    py[i] = py[i]*cosxidt + Pp[1]*sinxidt;
    pz[i] = pz[i]*cosxidt + Pp[2]*sinxidt;
    pw[i] = pw[i]*cosxidt + Pp[3]*sinxidt;

    if(((MODE>>TSHIFT)&MASK)==T_CONSTANT){
      const int ii = i%nlane;
      sum_Ps_omp[nlane*thread+ii] += 0.5*xi*xi*iy[i];
    }
  }
}

template <int MODE>
inline void Molecules::D3(){
#ifdef _OPENMP
  const int thread = omp_get_thread_num();
#else
  const int thread = 0;
#endif
  const double coef_xi = 0.25 / tst->s;
  //#pragma omp simd // make slower
  const int ivs = is[thread];
  const int ive = is[thread] + ((ie[thread] - is[thread])/nlane)*nlane;
  double* sum_Ps = &sum_Ps_omp[nlane*thread];
  for(int iv=ivs;iv<ive;iv+=nlane){
    #pragma omp simd
    for(int ii=0;ii<nlane;ii++){
      const int i = iv + ii;

      dvec4  Pq; Pq[0] = -qw[i]; Pq[1] = qz[i]; Pq[2] = -qy[i]; Pq[3] = qx[i];
      dvec4  Pp; Pp[0] = -pw[i]; Pp[1] = pz[i]; Pp[2] = -py[i]; Pp[3] = px[i];
      const double xi = (px[i]*Pq[0] + py[i]*Pq[1] + pz[i]*Pq[2] + pw[i]*Pq[3])*coef_xi/iz[i];
      const double xidt = xi*dthalf;
      const double cosxidt = cos(xidt);
      const double sinxidt = sin(xidt);
      qx[i] = qx[i]*cosxidt + Pq[0]*sinxidt;
      qy[i] = qy[i]*cosxidt + Pq[1]*sinxidt;
      qz[i] = qz[i]*cosxidt + Pq[2]*sinxidt;
      qw[i] = qw[i]*cosxidt + Pq[3]*sinxidt;

      px[i] = px[i]*cosxidt + Pp[0]*sinxidt;
      py[i] = py[i]*cosxidt + Pp[1]*sinxidt;
      pz[i] = pz[i]*cosxidt + Pp[2]*sinxidt;
      pw[i] = pw[i]*cosxidt + Pp[3]*sinxidt;
      if(((MODE>>TSHIFT)&MASK)==T_CONSTANT){
	sum_Ps_omp[nlane*thread+ii] += 0.5*xi*xi*iz[i];
      }
    }
  }
  for(int i=ive;i<ie[thread];i++){
    dvec4  Pq; Pq[0] = -qw[i]; Pq[1] = qz[i]; Pq[2] = -qy[i]; Pq[3] = qx[i];
    dvec4  Pp; Pp[0] = -pw[i]; Pp[1] = pz[i]; Pp[2] = -py[i]; Pp[3] = px[i];
    const double xi = (px[i]*Pq[0] + py[i]*Pq[1] + pz[i]*Pq[2] + pw[i]*Pq[3])*coef_xi/iz[i];
    const double xidt = xi*dthalf;
    const double cosxidt = cos(xidt);
    const double sinxidt = sin(xidt);
    qx[i] = qx[i]*cosxidt + Pq[0]*sinxidt;
    qy[i] = qy[i]*cosxidt + Pq[1]*sinxidt;
    qz[i] = qz[i]*cosxidt + Pq[2]*sinxidt;
    qw[i] = qw[i]*cosxidt + Pq[3]*sinxidt;

    px[i] = px[i]*cosxidt + Pp[0]*sinxidt;
    py[i] = py[i]*cosxidt + Pp[1]*sinxidt;
    pz[i] = pz[i]*cosxidt + Pp[2]*sinxidt;
    pw[i] = pw[i]*cosxidt + Pp[3]*sinxidt;
    if(((MODE>>TSHIFT)&MASK)==T_CONSTANT){
      const int ii = i%nlane;
      sum_Ps_omp[nlane*thread+ii] += 0.5*xi*xi*iz[i];
    }
  }
}

template <int MODE>
inline void Molecules::D5(){
#ifdef _OPENMP
  const int thread = omp_get_thread_num();
#else
  const int thread = 0;
#endif

  const dvec3  coef_v = L * tst->s * dthalf;
  const double coef_p = 2.0 * tst->s * dthalf;

  const int ivs = is[thread];
  const int ive = is[thread] + ((ie[thread] - is[thread])/nlane)*nlane;
  double* sum_Ps = &sum_Ps_omp[nlane*thread];
  for(int iv=ivs;iv<ive;iv+=nlane){
    #pragma omp simd
    for(int ii=0;ii<nlane;ii++){
      const int i = iv + ii;

      vx[i] += fx[i] * coef_v[0] * mi[i];
      vy[i] += fy[i] * coef_v[1] * mi[i];
      vz[i] += fz[i] * coef_v[2] * mi[i];
      dvec4 S0; S0[0] = qx[i]; S0[1] = -qy[i]; S0[2] = -qz[i]; S0[3] = -qw[i];
      dvec4 S1; S1[0] = qy[i]; S1[1] =  qx[i]; S1[2] = -qw[i]; S1[3] =  qz[i];
      dvec4 S2; S2[0] = qz[i]; S2[1] =  qw[i]; S2[2] =  qx[i]; S2[3] = -qy[i];
      dvec4 S3; S3[0] = qw[i]; S3[1] = -qz[i]; S3[2] =  qy[i]; S3[3] =  qx[i];
      px[i] += (S0[0]*nx[i] + S0[1]*ny[i] + S0[2]*nz[i] + S0[3]*nw[i])*coef_p;
      py[i] += (S1[0]*nx[i] + S1[1]*ny[i] + S1[2]*nz[i] + S1[3]*nw[i])*coef_p;
      pz[i] += (S2[0]*nx[i] + S2[1]*ny[i] + S2[2]*nz[i] + S2[3]*nw[i])*coef_p;
      pw[i] += (S3[0]*nx[i] + S3[1]*ny[i] + S3[2]*nz[i] + S3[3]*nw[i])*coef_p;
    }
  }
  for(int i=ive;i<ie[thread];i++){
    vx[i] += fx[i] * coef_v[0] * mi[i];
    vy[i] += fy[i] * coef_v[1] * mi[i];
    vz[i] += fz[i] * coef_v[2] * mi[i];
    dvec4 S0; S0[0] = qx[i]; S0[1] = -qy[i]; S0[2] = -qz[i]; S0[3] = -qw[i];
    dvec4 S1; S1[0] = qy[i]; S1[1] =  qx[i]; S1[2] = -qw[i]; S1[3] =  qz[i];
    dvec4 S2; S2[0] = qz[i]; S2[1] =  qw[i]; S2[2] =  qx[i]; S2[3] = -qy[i];
    dvec4 S3; S3[0] = qw[i]; S3[1] = -qz[i]; S3[2] =  qy[i]; S3[3] =  qx[i];
    px[i] += (S0[0]*nx[i] + S0[1]*ny[i] + S0[2]*nz[i] + S0[3]*nw[i])*coef_p;
    py[i] += (S1[0]*nx[i] + S1[1]*ny[i] + S1[2]*nz[i] + S1[3]*nw[i])*coef_p;
    pz[i] += (S2[0]*nx[i] + S2[1]*ny[i] + S2[2]*nz[i] + S2[3]*nw[i])*coef_p;
    pw[i] += (S3[0]*nx[i] + S3[1]*ny[i] + S3[2]*nz[i] + S3[3]*nw[i])*coef_p;
  }
#pragma omp single
  {
    if(((MODE>>TSHIFT)&MASK)==T_CONSTANT){
      tst->Ps -= prop.pot * dthalf;
      if(((MODE>>PSHIFT)&MASK)>0){
	tst->Ps -= (P*Molecules::GetVolume())*dthalf;
      }
    }
    if(((MODE>>PSHIFT)&MASK)==1){
      bst->Pv[0] += (sum(prop.vir)/(3.0*L[0]*L[1]*L[2]) - P)*tst->s*dthalf;
    }
    if(((MODE>>PSHIFT)&MASK)==2){
      bst->Pv[2] += tst->s * (prop.vir[2]/L[2] - P*Molecules::GetBottomArea()) * dthalf;
      //printf("(vir,P)= %lf %lf\n",prop.vir[2]/L[2],P*Molecules::GetBottomArea());
    }
    if(((MODE>>PSHIFT)&MASK)==3){
      bst->Pv[0] += (prop.vir[0]/L[0] - P*L[2]*param.wall_length)*tst->s*dthalf;
      bst->Pv[1] += (prop.vir[1]/L[1] - P*L[2]*param.wall_length)*tst->s*dthalf;
    }
  }
}

#else
template <int MODE>
inline void Molecules::D1(){
#ifdef _OPENMP
  const int thread = omp_get_thread_num();
#else
  const int thread = 0;
#endif

  const double si   = 1.0/tst->s;
  const dvec3  Li   = L.inv();
  const dvec3  coef_r  = Li*Li*si*dt;
  const double coef_xi = 0.25*si;
  //const double coef_tra = si*si*0.5;
  //double sum_tra = 0.0;
  //#pragma omp simd // need to check sum_Pv_omp summation is collect or not
  for(int i=is[thread];i<ie[thread];i++){
    Molecule& m = mlcl[i];
    //update coordinate
    m.r += m.v*coef_r;
    for(int d=0;d<3;d++){
      if(m.r[d] <  0.0) m.r[d] += 1.0;
      if(m.r[d] >= 1.0) m.r[d] -= 1.0;
    }

    //update angle and anglular velocisty
    const dvec4  Pq = P1(m.q);
    const dvec4  Pp = P1(m.p);
    const double xi = scalar_prod(m.p,Pq)*coef_xi/m.i[0];
    const double xidt = xi*dt;
    m.q = m.q*cos(xidt) + Pq*sin(xidt);
    m.p = m.p*cos(xidt) + Pp*sin(xidt);
    //update molecule
    //mlcl[i] = m;
    if constexpr (((MODE>>PSHIFT)&MASK)>0 || ((MODE>>TSHIFT)&MASK)==T_CONSTANT){
      const int ii = i%nlane;
      //sum_Pv += m.v*m.v*Li*Li*m.m;
      sum_Pv_omp[nlane*thread+ii][0] += 2.0 * m.m * m.v[0]*m.v[0];
      sum_Pv_omp[nlane*thread+ii][1] += 2.0 * m.m * m.v[1]*m.v[1];
      sum_Pv_omp[nlane*thread+ii][2] += 2.0 * m.m * m.v[2]*m.v[2];
    }
    if constexpr (((MODE>>TSHIFT)&MASK)==T_CONSTANT){
      const int ii = i%nlane;
      sum_Ps_omp[nlane*thread+ii] += xi*xi*m.i[0];
    }
  }
}

template <int MODE>
inline void Molecules::D2(){
#ifdef _OPENMP
  const int thread = omp_get_thread_num();
#else
  const int thread = 0;
#endif

  const double coef_xi = 0.25 / tst->s;
  //#pragma omp simd // make slower
  for(int i=is[thread];i<ie[thread];i++){
    Molecule m = mlcl[i];
    const dvec4  Pq = P2(m.q);
    const dvec4  Pp = P2(m.p);
    const double xi = scalar_prod(m.p,Pq)*coef_xi/m.i[1];
    const double xidt = xi*dthalf;
    m.q = m.q*cos(xidt) + Pq*sin(xidt);
    m.p = m.p*cos(xidt) + Pp*sin(xidt);
    mlcl[i] = m;
    if(((MODE>>TSHIFT)&MASK)==T_CONSTANT){
      const int ii = i%nlane;
      sum_Ps_omp[nlane*thread+ii] += 0.5*xi*xi*m.i[1];
    }
  }
}

template <int MODE>
inline void Molecules::D3(){
#ifdef _OPENMP
  const int thread = omp_get_thread_num();
#else
  const int thread = 0;
#endif
  const double coef_xi = 0.25 / tst->s;
  //#pragma omp simd // make slower
  for(int i=is[thread];i<ie[thread];i++){
    Molecule m = mlcl[i];
    const dvec4  Pq = P3(m.q);
    const dvec4  Pp = P3(m.p);
    const double xi = scalar_prod(m.p,Pq)*coef_xi/m.i[2];
    const double xidt = xi*dthalf;
    m.q = m.q*cos(xidt) + Pq*sin(xidt);
    m.p = m.p*cos(xidt) + Pp*sin(xidt);
    mlcl[i] = m;
    if(((MODE>>TSHIFT)&MASK)==T_CONSTANT){
      const int ii = i%nlane;
      sum_Ps_omp[nlane*thread+ii] += 0.5*xi*xi*m.i[2]; // 0.5 for dthalf
    }
  }
}
template <int MODE>
inline void Molecules::D5(){
#ifdef _OPENMP
  const int thread = omp_get_thread_num();
#else
  const int thread = 0;
#endif

  const dvec3  coef_v = L * tst->s * dthalf;
  const double coef_p = 2.0 * tst->s * dthalf;

  for(int i=is[thread];i<ie[thread];i++){
    Molecule m = mlcl[i];
    m.v += m.f * coef_v / m.m;
    m.p += scalar_prod(S(m.q),m.n)*coef_p;
    mlcl[i] = m;
  }
#pragma omp single
  {
    if(((MODE>>TSHIFT)&MASK)==T_CONSTANT){
      tst->Ps -= prop.pot * dthalf;
      if(((MODE>>PSHIFT)&MASK)>0){
	tst->Ps -= (P*Molecules::GetVolume())*dthalf;
      }
    }
    if(((MODE>>PSHIFT)&MASK)==1){
      bst->Pv[0] += (sum(prop.vir)/(3.0*L[0]*L[1]*L[2]) - P)*tst->s*dthalf;
    }
    if(((MODE>>PSHIFT)&MASK)==2){
      bst->Pv[2] += tst->s * (prop.vir[2]/L[2] - P*Molecules::GetBottomArea()) * dthalf;
      //printf("(vir,P)= %lf %lf\n",prop.vir[2]/L[2],P*Molecules::GetBottomArea());
    }
    if(((MODE>>PSHIFT)&MASK)==3){
      bst->Pv[0] += (prop.vir[0]/L[0] - P*L[2]*param.wall_length)*tst->s*dthalf;
      bst->Pv[1] += (prop.vir[1]/L[1] - P*L[2]*param.wall_length)*tst->s*dthalf;
    }
  }
}
#endif

template <int MODE>
void Molecules::D4(){
  #pragma omp single
  {
  //* for p constant
  if(((MODE>>PSHIFT)&MASK)==1){
    if(((MODE>>TSHIFT)&MASK)==1){
      tst->Ps -= bst->Pv[0]*bst->Pv[0]/bst->W*dthalf;
    }
    double V = Molecules::GetVolume();
    L = powf(V + bst->Pv[0]*tst->s/bst->W*dthalf,1./3.);
  }
  if(((MODE>>PSHIFT)&MASK)==2){
    if(((MODE>>TSHIFT)&MASK)==1){
      tst->Ps -= 0.5 * bst->Pv[2] * bst->Pv[2] / bst->W * dthalf;
    }
    L[2] += tst->s * bst->Pv[2] / bst->W * dthalf;
  }
  if(((MODE>>PSHIFT)&MASK)==3){
    if(((MODE>>TSHIFT)&MASK)==1){
      tst->Ps -= (bst->Pv[0]*bst->Pv[0] + bst->Pv[1]*bst->Pv[1])/bst->W*dthalf;
    }
    L[0] += bst->Pv[0]*tst->s/bst->W*dthalf;
    L[1] += bst->Pv[1]*tst->s/bst->W*dthalf;
  }
  }
  //printf("(L.z,Pv.z,Ps,s)= %lf %lf %lf %lf\n",L[2],bst->Pv[2],tst->Ps,tst->s);
  //*/
}

template <int MODE>
void Molecules::D6(){
  #pragma omp single
  {
    if(((MODE>>TSHIFT)&MASK)==T_CONSTANT){
      const double tmp = 1.0 + ( tst->Ps / (2.0*tst->Q) * dthalf );
      tst->s  *= tmp*tmp;
      tst->Ps /= tmp;
    }
  }
}

template<int MODE>
void Molecules::ExecuteStep(const bool doSort,const double mergin){
    prof.beg(Profiler::Total);
    prof.beg(Profiler::Integ);
    prof.beg(Profiler::D6);
    D6<MODE>();
    prof.end(Profiler::D6);
    prof.beg(Profiler::D5);
    D5<MODE>();
    prof.end(Profiler::D5);
    prof.beg(Profiler::D4);
    D4<MODE>();
    prof.end(Profiler::D4);
    prof.beg(Profiler::PsPv);
    init_D3D2D1D2D3_PsPv<MODE>();
    prof.end(Profiler::PsPv);
    prof.beg(Profiler::D3);
    D3<MODE>();
    prof.end(Profiler::D3);
    prof.beg(Profiler::D2);
    D2<MODE>();
    prof.end(Profiler::D2);
    prof.beg(Profiler::D1);
    D1<MODE>();
    prof.end(Profiler::D1);
    prof.beg(Profiler::D2);
    D2<MODE>();
    prof.end(Profiler::D2);
    prof.beg(Profiler::D3);
    D3<MODE>();
    prof.end(Profiler::D3);
    prof.beg(Profiler::PsPv);
    D3D2D1D2D3_PsPv<MODE>();
    prof.end(Profiler::PsPv);
    prof.beg(Profiler::D4);
    D4<MODE>();
    prof.end(Profiler::D4);
    prof.end(Profiler::Integ);
#ifdef ENABLE_AOS_TO_SOA_CONVERSION
    prof.beg(Profiler::SoAtoAoS);
    SoAtoAoS();
    prof.end(Profiler::SoAtoAoS);
#endif
    prof.beg(Profiler::CalcForce);
    CalcForcePot(doSort,mergin);
    prof.end(Profiler::CalcForce);
#ifdef ENABLE_AOS_TO_SOA_CONVERSION
    prof.beg(Profiler::AoStoSoA);
    AoStoSoA();
    prof.end(Profiler::AoStoSoA);
#endif
    prof.beg(Profiler::Integ);
    prof.beg(Profiler::D5);
    D5<MODE>();
    prof.end(Profiler::D5);
    prof.beg(Profiler::D6);
    D6<MODE>();
    prof.end(Profiler::D6);
    prof.end(Profiler::Integ);
    prof.end(Profiler::Total);
}

void Molecules::ExecuteSteps(const int sort_interval,const double mergin){
  #pragma omp parallel
  {
#ifdef ENABLE_AOS_TO_SOA_CONVERSION
    AoStoSoA();
#endif
    for(int s=0;s<param.interval;s++){
#pragma omp single
      {
	if(param.tconstant == 2){
#ifdef ENABLE_AOS_TO_SOA_CONVERSION
	  //SoAtoAoS();
#endif
	  VelocityScaling();
	  AngVelocityScaling();
#ifdef ENABLE_AOS_TO_SOA_CONVERSION
	  //AoStoSoA();
#endif
	}
      }
      //printf("mode= %d %d %d\n",(mode>>TSHIFT)&MASK, (mode>>PSHIFT)&MASK, (mode>>CSHIFT)&MASK);
      const bool doSort = s%sort_interval == 0;
      switch(mode){
      case NVE:
	ExecuteStep<NVE>(doSort,mergin);break;
      case NVT:
	ExecuteStep<NVT>(doSort,mergin);break;
      case NVTSCALE:
	ExecuteStep<NVE>(doSort,mergin);break;
      case NPH:
	ExecuteStep<NPH>(doSort,mergin);break;
      case NPT:
	ExecuteStep<NPT>(doSort,mergin);break;
      case NAxyPzH:
	ExecuteStep<NAxyPzH>(doSort,mergin);break;
      case NAxyPzT:
	ExecuteStep<NAxyPzT>(doSort,mergin);break;
      case NPxyLzH:
	ExecuteStep<NPxyLzH>(doSort,mergin);break;
      case NPxyLzT:
	ExecuteStep<NPxyLzT>(doSort,mergin);break;
      case NVE1D:
	ExecuteStep<NVE1D>(doSort,mergin);break;
      case NVT1D:
	ExecuteStep<NVT1D>(doSort,mergin);break;
      case NVTSCALE1D:
	ExecuteStep<NVE1D>(doSort,mergin);break;
      case NPH1D:
	ExecuteStep<NPH1D>(doSort,mergin);break;
      case NPT1D:
	ExecuteStep<NPT1D>(doSort,mergin);break;
      case NPTSCALE1D:
	ExecuteStep<NPH1D>(doSort,mergin);break;
      case NAxyPzH1D:
	ExecuteStep<NAxyPzH1D>(doSort,mergin);break;
      case NAxyPzT1D:
	ExecuteStep<NAxyPzT1D>(doSort,mergin);break;
      case NAxyPzTSCALE1D:
	ExecuteStep<NAxyPzH1D>(doSort,mergin);break;
      case NPxyLzH1D:
	ExecuteStep<NPxyLzH1D>(doSort,mergin);break;
      case NPxyLzT1D:
	ExecuteStep<NPxyLzT1D>(doSort,mergin);break;
      case NVE2D:
	ExecuteStep<NVE2D>(doSort,mergin);break;
      case NVT2D:
	ExecuteStep<NVT2D>(doSort,mergin);break;
      case NPH2D:
	ExecuteStep<NPH2D>(doSort,mergin);break;
      case NPT2D:
	ExecuteStep<NPT2D>(doSort,mergin);break;
      case NPTSCALE2D:
	ExecuteStep<NPH2D>(doSort,mergin);break;
      case NAxyPzH2D:
	ExecuteStep<NAxyPzH2D>(doSort,mergin);break;
      case NAxyPzT2D:
	ExecuteStep<NAxyPzT2D>(doSort,mergin);break;
      case NPxyLzH2D:
	ExecuteStep<NPxyLzH2D>(doSort,mergin);break;
      case NPxyLzT2D:
	ExecuteStep<NPxyLzT2D>(doSort,mergin);break;
      defalut:
	std::cerr << "error: undefined ensemble or confined dimention" << std::endl;
	exit(EXIT_FAILURE);
      } // switch
#pragma omp single
      {
	prop.time += unit_time * dt;
      } // omp single
    } // s loop
#ifdef ENABLE_AOS_TO_SOA_CONVERSION
    SoAtoAoS();
#endif
  } // omp parallel
  CalcProperties();
}
