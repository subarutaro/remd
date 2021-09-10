#include "molecules.h"
#include "integrator.h"

template <int MODE>
void Molecules::D1(){
  const double si   = 1.0/tst->s;
  const dvec3  Li   = L.inv();
  const dvec3  coef_r  = Li*Li*si*dt;
  const double coef_xi = 0.25*si;
  //const double coef_tra = si*si*0.5;
  //double sum_tra = 0.0;
  double sum_rot = 0.0;
  dvec3 sum_Pv = 0.0;
  for(int i=0;i<nmol;i++){
    Molecule m = mlcl[i];
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
    mlcl[i] = m;
    if(((MODE>>PSHIFT)&MASK)>0){
      //sum_Pv += m.v*m.v*Li*Li*m.m;
      sum_Pv[0] += m.m * m.v[0]*m.v[0];
      sum_Pv[1] += m.m * m.v[1]*m.v[1];
      sum_Pv[2] += m.m * m.v[2]*m.v[2];
    }
    if(((MODE>>TSHIFT)&MASK)==T_CONSTANT){
      sum_rot += xi*xi*m.i[0];
    }
  }
  if(((MODE>>PSHIFT)&MASK)>0){
    //sum_Pv += m.v*m.v*Li*Li*m.m;
    sum_Pv[0] /= L[0]*L[0] * tst->s*tst->s;
    sum_Pv[1] /= L[1]*L[1] * tst->s*tst->s;
    sum_Pv[2] /= L[2]*L[2] * tst->s*tst->s;
  }

  if(((MODE>>TSHIFT)&MASK)==T_CONSTANT){
    tst->Ps += ( sum( TranslationalEnergy() ) + 2.0*sum_rot - ( 1.0+log(tst->s) )*prop.gkT + prop.H0 ) * dt;
  }
  if(((MODE>>PSHIFT)&MASK)==1){
    bst->Pv[0] += tst->s * sum(sum_Pv)/(Molecules::GetVolume()*3.0) * dt;
  }
  if(((MODE>>PSHIFT)&MASK)==2){
    bst->Pv[2] += tst->s * sum_Pv[2] / L[2] * dt;
    //printf("(kin)= %lff\n",sum_Pv[2] / L[2]);
  }
  if(((MODE>>PSHIFT)&MASK)==3){
    bst->Pv[0] += tst->s * sum_Pv[0] / L[0] * dt;
    bst->Pv[1] += tst->s * sum_Pv[1] / L[1] * dt;
  }
}

template <int MODE>
void Molecules::D2(){
  const double coef_xi = 0.25 / tst->s;
  double sum_Ps = 0.;
  for(int i=0;i<nmol;i++){
    Molecule m = mlcl[i];
    const dvec4  Pq = P2(m.q);
    const dvec4  Pp = P2(m.p);
    const double xi = scalar_prod(m.p,Pq)*coef_xi/m.i[1];
    const double xidt = xi*dthalf;
    m.q = m.q*cos(xidt) + Pq*sin(xidt);
    m.p = m.p*cos(xidt) + Pp*sin(xidt);
    mlcl[i] = m;
    if(((MODE>>TSHIFT)&MASK)==T_CONSTANT){
      sum_Ps += xi*xi*m.i[1];
    }
  }
  if(((MODE>>TSHIFT)&MASK)==T_CONSTANT){
    //std::cout << "D2 sum_Ps: " << 2.0*sum_Ps << " s: " << tst->s << ", Ps: " << tst->Ps << std::endl;
    tst->Ps += 2.0*sum_Ps*dthalf;
  }
}

template <int MODE>
void Molecules::D3(){
  const double coef_xi = 0.25 / tst->s;
  double sum_Ps = 0.0;
  for(int i=0;i<nmol;i++){
    Molecule m = mlcl[i];
    const dvec4  Pq = P3(m.q);
    const dvec4  Pp = P3(m.p);
    const double xi = scalar_prod(m.p,Pq)*coef_xi/m.i[2];
    const double xidt = xi*dthalf;
    m.q = m.q*cos(xidt) + Pq*sin(xidt);
    m.p = m.p*cos(xidt) + Pp*sin(xidt);
    mlcl[i] = m;
    if(((MODE>>TSHIFT)&MASK)==T_CONSTANT){
      sum_Ps += xi*xi*m.i[2];
    }
  }
  if(((MODE>>TSHIFT)&MASK)==T_CONSTANT){
    //std::cout << "D3 sum_Ps: " << 2.0*sum_Ps << " s: " << tst->s << ", Ps: " << tst->Ps << std::endl;
    tst->Ps += 2.0*sum_Ps*dthalf;
  }
}

template <int MODE>
void Molecules::D4(){
  //* for p constant
  if(((MODE>>PSHIFT)&MASK)==1){
    if(((MODE>>TSHIFT)&MASK)==1){
      tst->Ps -= 0.5*bst->Pv[0]*bst->Pv[0]/bst->W*dthalf;
    }
    double V = Molecules::GetVolume();
    L = powf(V + bst->Pv[0]*tst->s/bst->W*dthalf,1./3.);
  }
  if(((MODE>>PSHIFT)&MASK)==2){
    if(((MODE>>TSHIFT)&MASK)==1){
      tst->Ps -= 0.5 * bst->Pv[2] * bst->Pv[2] / bst->W * dthalf;
    }
    L[2] += 0.5 * tst->s * bst->Pv[2] / bst->W * dthalf;
  }
  if(((MODE>>PSHIFT)&MASK)==3){
    if(((MODE>>TSHIFT)&MASK)==1){
      tst->Ps -= 0.5*(bst->Pv[0]*bst->Pv[0] + bst->Pv[1]*bst->Pv[1])/bst->W*dthalf;
    }
    L[0] += bst->Pv[0]*tst->s/bst->W*dthalf;
    L[1] += bst->Pv[1]*tst->s/bst->W*dthalf;
  }
  //printf("(L.z,Pv.z,Ps,s)= %lf %lf %lf %lf\n",L[2],bst->Pv[2],tst->Ps,tst->s);
  //*/
}

template <int MODE>
void Molecules::D5(){
  const dvec3  coef_v = L * tst->s * dthalf;
  const double coef_p = 2.0 * tst->s * dthalf;
  //#pragma omp parallel for
  for(int i=0;i<nmol;i++){
    Molecule m = mlcl[i];
    m.v += m.f * coef_v / m.m;
    m.p += scalar_prod(S(m.q),m.n)*coef_p;
    mlcl[i] = m;
  }
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

template <int MODE>
void Molecules::D6(){
  if(((MODE>>TSHIFT)&MASK)==T_CONSTANT){
    const double tmp = 1.0 + ( tst->Ps / (2.0*tst->Q) * dthalf );
    tst->s  *= tmp*tmp;
    tst->Ps /= tmp;
  }
}

void Molecules::ExecuteSteps(){
  for(int s=0;s<param.interval;s++){
    if(param.tconstant == 2){
      VelocityScaling();
      AngVelocityScaling();
    }
#define INTEGRATE(x)    \
    D6<x>();		\
    D5<x>();		\
    D4<x>();		\
    D3<x>();		\
    D2<x>();		\
    D1<x>();		\
    D2<x>();		\
    D3<x>();		\
    D4<x>();		\
    CalcForcePot();	\
    D5<x>();		\
    D6<x>();
    //printf("mode= %d %d %d\n",(mode>>TSHIFT)&MASK, (mode>>PSHIFT)&MASK, (mode>>CSHIFT)&MASK);
    switch(mode){
    case NVE:
      INTEGRATE(NVE);break;
    case NVT:
      INTEGRATE(NVT);break;
    case NVTSCALE:
      INTEGRATE(NVE);break;
    case NPH:
      INTEGRATE(NPH);break;
    case NPT:
      INTEGRATE(NPT);break;
    case NAxyPzH:
      INTEGRATE(NAxyPzH);break;
    case NAxyPzT:
      INTEGRATE(NAxyPzT);break;
    case NPxyLzH:
      INTEGRATE(NPxyLzH);break;
    case NPxyLzT:
      INTEGRATE(NPxyLzT);break;
    case NVE1D:
      INTEGRATE(NVE1D);break;
    case NVT1D:
      INTEGRATE(NVT1D);break;
    case NVTSCALE1D:
      INTEGRATE(NVE1D);break;
    case NPH1D:
      INTEGRATE(NPH1D);break;
    case NPT1D:
      INTEGRATE(NPT1D);break;
    case NPTSCALE1D:
      INTEGRATE(NPH1D);break;
    case NAxyPzH1D:
      INTEGRATE(NAxyPzH1D);break;
    case NAxyPzT1D:
      INTEGRATE(NAxyPzT1D);break;
    case NAxyPzTSCALE1D:
      INTEGRATE(NAxyPzH1D);break;
    case NPxyLzH1D:
      INTEGRATE(NPxyLzH1D);break;
    case NPxyLzT1D:
      INTEGRATE(NPxyLzT1D);break;
    case NVE2D:
      INTEGRATE(NVE2D);break;
    case NVT2D:
      INTEGRATE(NVT2D);break;
    case NPH2D:
      INTEGRATE(NPH2D);break;
    case NPT2D:
      INTEGRATE(NPT2D);break;
    case NPTSCALE2D:
      INTEGRATE(NPH2D);break;
    case NAxyPzH2D:
      INTEGRATE(NAxyPzH2D);break;
    case NAxyPzT2D:
      INTEGRATE(NAxyPzT2D);break;
    case NPxyLzH2D:
      INTEGRATE(NPxyLzH2D);break;
    case NPxyLzT2D:
      INTEGRATE(NPxyLzT2D);break;
    defalut:
      std::cerr << "error: undefined ensemble or confined dimention" << std::endl;
      exit(EXIT_FAILURE);
    }
#undef INTEGRATE
    prop.time += unit_time * dt;
  }
  CalcProperties();
}
