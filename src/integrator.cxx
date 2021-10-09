#include "molecules.h"
#include "integrator.h"

template <int MODE>
inline void Molecules::init_D3D2D1D2D3_PsPv(){
#ifdef _OPENMP
  const int thread = omp_get_thread_num();
#else
  const int thread = 0;
#endif
  sum_Pv_omp[thread] = 0.0;
  sum_Ps_omp[thread] = 0.0;
}
template <int MODE>
inline void Molecules::D3D2D1D2D3_PsPv(){
  #pragma omp barrier
  #pragma omp single
  {
    if(((MODE>>PSHIFT)&MASK)>0){
      dvec3  sum_Pv  = 0.0;
      for(int i=0;i<nthreads;i++) sum_Pv  += sum_Pv_omp[i];
      //sum_Pv += m.v*m.v*Li*Li*m.m;
      sum_Pv[0] /= L[0]*L[0] * tst->s*tst->s;
      sum_Pv[1] /= L[1]*L[1] * tst->s*tst->s;
      sum_Pv[2] /= L[2]*L[2] * tst->s*tst->s;

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
    }

    double sum_Ps = 0.0;
    for(int i=0;i<nthreads;i++) sum_Ps += sum_Ps_omp[i];
    if(((MODE>>TSHIFT)&MASK)==T_CONSTANT){
      //std::cout << "D2 sum_Ps: " << 2.0*sum_Ps << " s: " << tst->s << ", Ps: " << tst->Ps << std::endl;
      tst->Ps += 2.0*sum_Ps*dthalf;
    }

    if(((MODE>>TSHIFT)&MASK)==T_CONSTANT){
      double sum_Ps = 0.0;
      for(int i=0;i<nthreads;i++) sum_Ps += sum_Ps_omp[i];
      tst->Ps += ( sum( TranslationalEnergy() ) + 2.0*sum_Ps - ( 1.0+log(tst->s) )*prop.gkT + prop.H0 ) * dthalf;
    }
  }
}

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
  #pragma omp simd
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
    if constexpr (((MODE>>PSHIFT)&MASK)>0){
      //sum_Pv += m.v*m.v*Li*Li*m.m;
      sum_Pv_omp[thread][0] += 2.0 * m.m * m.v[0]*m.v[0];
      sum_Pv_omp[thread][1] += 2.0 * m.m * m.v[1]*m.v[1];
      sum_Pv_omp[thread][2] += 2.0 * m.m * m.v[2]*m.v[2];
    }
    if constexpr (((MODE>>TSHIFT)&MASK)==T_CONSTANT){
      sum_Ps_omp[thread] += 2.0*xi*xi*m.i[0];
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
      sum_Ps_omp[thread] += xi*xi*m.i[1];
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
      sum_Ps_omp[thread] += xi*xi*m.i[2];
    }
  }
}

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
void Molecules::ExecuteStep(){
    prof.beg(Profiler::Total);
    prof.beg(Profiler::D6);
    D6<MODE>();
    prof.end(Profiler::D6);
    prof.beg(Profiler::D5);
    D5<MODE>();
    prof.end(Profiler::D5);
    prof.beg(Profiler::D4);
    D4<MODE>();
    prof.end(Profiler::D4);
    init_D3D2D1D2D3_PsPv<MODE>();
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
    D3D2D1D2D3_PsPv<MODE>();
    prof.beg(Profiler::D4);
    D4<MODE>();
    prof.end(Profiler::D4);
    prof.beg(Profiler::CalcForce);
    CalcForcePot();
    prof.end(Profiler::CalcForce);
    prof.beg(Profiler::D5);
    D5<MODE>();
    prof.end(Profiler::D5);
    prof.beg(Profiler::D6);
    D6<MODE>();
    prof.end(Profiler::D6);
    prof.end(Profiler::Total);

}

void Molecules::ExecuteSteps(){
  #pragma omp parallel
  {
  for(int s=0;s<param.interval;s++){
    #pragma omp single
    {
    if(param.tconstant == 2){
      VelocityScaling();
      AngVelocityScaling();
    }
    }
    //printf("mode= %d %d %d\n",(mode>>TSHIFT)&MASK, (mode>>PSHIFT)&MASK, (mode>>CSHIFT)&MASK);
    switch(mode){
    case NVE:
      ExecuteStep<NVE>();break;
    case NVT:
      ExecuteStep<NVT>();break;
    case NVTSCALE:
      ExecuteStep<NVE>();break;
    case NPH:
      ExecuteStep<NPH>();break;
    case NPT:
      ExecuteStep<NPT>();break;
    case NAxyPzH:
      ExecuteStep<NAxyPzH>();break;
    case NAxyPzT:
      ExecuteStep<NAxyPzT>();break;
    case NPxyLzH:
      ExecuteStep<NPxyLzH>();break;
    case NPxyLzT:
      ExecuteStep<NPxyLzT>();break;
    case NVE1D:
      ExecuteStep<NVE1D>();break;
    case NVT1D:
      ExecuteStep<NVT1D>();break;
    case NVTSCALE1D:
      ExecuteStep<NVE1D>();break;
    case NPH1D:
      ExecuteStep<NPH1D>();break;
    case NPT1D:
      ExecuteStep<NPT1D>();break;
    case NPTSCALE1D:
      ExecuteStep<NPH1D>();break;
    case NAxyPzH1D:
      ExecuteStep<NAxyPzH1D>();break;
    case NAxyPzT1D:
      ExecuteStep<NAxyPzT1D>();break;
    case NAxyPzTSCALE1D:
      ExecuteStep<NAxyPzH1D>();break;
    case NPxyLzH1D:
      ExecuteStep<NPxyLzH1D>();break;
    case NPxyLzT1D:
      ExecuteStep<NPxyLzT1D>();break;
    case NVE2D:
      ExecuteStep<NVE2D>();break;
    case NVT2D:
      ExecuteStep<NVT2D>();break;
    case NPH2D:
      ExecuteStep<NPH2D>();break;
    case NPT2D:
      ExecuteStep<NPT2D>();break;
    case NPTSCALE2D:
      ExecuteStep<NPH2D>();break;
    case NAxyPzH2D:
      ExecuteStep<NAxyPzH2D>();break;
    case NAxyPzT2D:
      ExecuteStep<NAxyPzT2D>();break;
    case NPxyLzH2D:
      ExecuteStep<NPxyLzH2D>();break;
    case NPxyLzT2D:
      ExecuteStep<NPxyLzT2D>();break;
    defalut:
      std::cerr << "error: undefined ensemble or confined dimention" << std::endl;
      exit(EXIT_FAILURE);
    } // switch
    #pragma omp single
    {
    prop.time += unit_time * dt;
    } // omp single
  } // s loop
  } // omp parallel
  CalcProperties();
}
