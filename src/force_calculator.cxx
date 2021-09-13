#include "molecules.h"
#include "force_calculator.h"

ForceCalculator::ForceCalculator
(double _rcut,double _kcut,double _alpha,int _nthreads)
{
  rcut     = _rcut;
  kcut     = _kcut;
  alpha    = _alpha;
  nthreads = _nthreads;
  SetKvec();
}

ForceCalculator::ForceCalculator
(const Parameter _p){
  rcut     = _p.rcut;
  kcut     = _p.kcut;
  alpha    = _p.alpha;
  rswitch  = _p.rswitch;
  assert(rcut>=rswitch);

  prefix   = _p.input_prefix;

  confined_dim = _p.confined;
  nfwall       = _p.nfwall;
  wall_length  = _p.wall_length*1e-10 / unit_length;
  sgm_wall     = _p.sgm_wall * 1e-10  / unit_length;
  eps_wall     = _p.eps_wall * 1e3/Na / unit_energy;
  rho_wall     = _p.rho_wall;
  if(confined_dim != 0){
    if(nfwall<=0){
      std::cerr << "error: set nfwall for calculation of wall force" << std::endl;
      exit(EXIT_FAILURE);
    }
    fwall = (dvec2*)malloc(nfwall*sizeof(dvec2));
    CalcWallForce();
  }
  nthreads  = _p.nthreads;
  SetKvec();
}

ForceCalculator::~ForceCalculator(){
  if(gx  != nullptr) _mm_free(gx);
  if(gy  != nullptr) _mm_free(gy);
  if(gz  != nullptr) _mm_free(gz);
  if(gfx != nullptr) _mm_free(gfx);
  if(gfy != nullptr) _mm_free(gfy);
  if(gfz != nullptr) _mm_free(gfz);
  if(gvx != nullptr) _mm_free(gvx);
  if(gvy != nullptr) _mm_free(gvy);
  if(gvz != nullptr) _mm_free(gvz);

  if(ax  != nullptr) _mm_free(ax);
  if(ay  != nullptr) _mm_free(ay);
  if(az  != nullptr) _mm_free(az);
  if(afx != nullptr) _mm_free(afx);
  if(afy != nullptr) _mm_free(afy);
  if(afz != nullptr) _mm_free(afz);
  if(as  != nullptr) _mm_free(as);
  if(ae  != nullptr) _mm_free(ae);
  if(aq  != nullptr) _mm_free(aq);

  free(kvec);
  free(fwall);
}

void ForceCalculator::SetKvec(){
  nwave = 0;
  const int kmax = (int)kcut;
  int memsize_kvec = sizeof(ivec3) * (kmax*2+1)*(kmax*2+1)*(kmax*2+1);
  kvec = (ivec3*)malloc(memsize_kvec);

  for(int i=-kmax;i<kmax+1;i++){
  for(int j=-kmax;j<kmax+1;j++){
  for(int k=-kmax;k<kmax+1;k++){
    int tmp = i*i + j*j + k*k;
    if(0 < tmp && (double)tmp < kcut*kcut){
      kvec[nwave][0] = i;
      kvec[nwave][1] = j;
      kvec[nwave][2] = k;
      nwave++;
    }
  }}}
}

void ForceCalculator::CalcWallForce(){
  const double dr = wall_length / (double)nfwall;
  if(confined_dim == 1){
    const int    ninteg = 10000;
    const double dphi   = M_PI/(double)ninteg;
    const double sgm03  = powf(sgm_wall, 3);
    const double sgm06  = powf(sgm_wall, 6);
    const double R02    = wall_length*wall_length;

#pragma omp parallel for
    for(int i=0;i<nfwall;i++){
      const double r = (double)i*dr;
      const double r02 = r*r;
      double tmp, sum = 0.0;
      for(int j=0;j<ninteg;j++){
	const double phi = M_PI - (double)j*dphi; // sum from far side
	//const double rho_min = sqrt( R02 - r02*sin(phi)*sin(phi) ) - r*cos(phi);
	const double rho_min = sqrtf( R02 + r02 - 2.0*wall_length*r*cos(phi) );
#if 0 // for cylindrical nanopore
	const double rho_min03i = 1.0 / (rho_min * rho_min * rho_min);
	tmp  = 7.0/32.0 * sgm06 * rho_min03i * rho_min03i * rho_min03i;
	tmp -= rho_min03i;
	sum += tmp;
#else // for single wall nanotube
	const double rho_mini   = 1.0 / rho_min;
	const double rho_min05i = rho_mini * rho_mini * rho_mini * rho_mini * rho_mini;
	tmp  = 63.0/32.0 * sgm06 * rho_min05i * rho_min05i * rho_mini;
	tmp -= 3.0 * rho_min05i;
	//sum += tmp*rho_min;
	sum += tmp;
#endif
	//printf("(wall) %e %e %e\n",r, rho_min, M_PI*rho_wall*sgm06*eps_wall*tmp);
      }
      //printf("tmp= %lf\n",tmp*dphi);
      //fwall[i][1] = M_PI * rho_wall * eps_wall * sgm06 * sum * dphi;
      fwall[i][1] = M_PI * rho_wall * eps_wall * sgm06 * sum * wall_length * dphi;
    }
    
    for(int i=1;i<nfwall-1;i++){
      const double r = (double)i*dr;
      fwall[i][0] = -(fwall[i+1][1]-fwall[i-1][1])/(2.0*dr * r);
    }
    fwall[0][0] = fwall[1][0];
    fwall[nfwall-1][0] = (fwall[nfwall-2][0] - fwall[nfwall-3][0])*2.0 + fwall[nfwall-3][0];
  }
  
  if(confined_dim == 2){
    FILE *fp;
    fp = fopen("force.dat","w");
    for(int i=0;i<nfwall;i++){
      const double rd    = 0.5*wall_length + 0.5*(double)i*dr;
      const double rdi   = 1.0/rd;

      const double ru    = 0.5*wall_length - 0.5*(double)i*dr;
      const double rui   = 1.0/ru;
#if 0
      const double rsdi  = sgm_wall*rdi;
      const double rsd03i = rsdi*rsdi*rsdi;
      const double rsd06i = rsd03i*rsd03i;

      const double rsui  = sgm_wall*rui;
      const double rsu03i = rsui*rsui*rsui;
      const double rsu06i = rsu03i*rsu03i;

      fwall[i][0]  = 4.0/15.0*M_PI*rho_wall*eps_wall*rdi*rsd03i*(3.0*rsd06i - 1.0);
      fwall[i][0] -= 4.0/15.0*M_PI*rho_wall*eps_wall*rui*rsu03i*(3.0*rsu06i - 1.0);
      fwall[i][1]  = 4.0/45.0*M_PI*rho_wall*eps_wall*rsd03i*(rsd06i - 1.0);
      fwall[i][1] += 4.0/45.0*M_PI*rho_wall*eps_wall*rsu03i*(rsu06i - 1.0);
#else
      const double rd03i = rdi*rdi*rdi;
      const double rd06i = rd03i*rd03i;

      const double ru03i = rui*rui*rui;
      const double ru06i = ru03i*ru03i;

      fwall[i][0]  = (9.0 * eps_wall * rd06i * rd03i - 3.0 * sgm_wall * rd03i) * rdi;
      fwall[i][0] -= (9.0 * eps_wall * ru06i * ru03i - 3.0 * sgm_wall * ru03i) * rui;
      fwall[i][1]  = eps_wall * rd06i * rd03i - sgm_wall * rd03i;
      fwall[i][1] += eps_wall * ru06i * ru03i - sgm_wall * ru03i;
      fprintf(fp,"%e %e %e\n", 0.5*i*dr, fwall[i][0], fwall[i][1]);
#endif
    }
  }
}


void ForceCalculator::MakeParamList(MolTypeList mlist){
  unsigned int natomtype = 0;
  for(unsigned int i=0;i<mlist.size();i++){
    for(unsigned int j=0;j<mlist[i].a.size();j++){
      natomtype++;
  }}
  //std::cout << "natomtype " << natomtype << std::endl;

  q = (double*)malloc(sizeof(double)*natomtype);

  sgm = (double**)malloc(sizeof(double*)*natomtype);
  eps = (double**)malloc(sizeof(double*)*natomtype);
  for(unsigned int i=0;i<natomtype;i++){
    sgm[i] = (double*)malloc(sizeof(double)*natomtype);
    eps[i] = (double*)malloc(sizeof(double)*natomtype);
  }

  unsigned int count_i = 0,count_j;
  for(unsigned int ii=0;ii<mlist.size();ii++){
    for(unsigned int ij=0;ij<mlist[ii].a.size();ij++){
      count_j = 0;
      for(unsigned int ji=0;ji<mlist.size();ji++){
	for(unsigned int jj=0;jj<mlist[ji].a.size();jj++){
	  sgm[count_i][count_j] = 0.5*(mlist[ii].a[ij].s + mlist[ji].a[jj].s);
	  sgm[count_i][count_j] = sgm[count_i][count_j]*sgm[count_i][count_j];
	  eps[count_i][count_j] = sqrt(mlist[ii].a[ij].e * mlist[ji].a[jj].e);
	  count_j++;
      }}
      q[count_i] = mlist[ii].a[ij].q;
      count_i++;
  }}
  /*
  std::cout << natomtype << std::endl;
  for(unsigned int i=0;i<natomtype;i++){
    std::cout << "q " << i << " " << q[i]<< std::endl;
    for(unsigned int j=0;j<natomtype;j++){
      std::cout << "sgm " << i << "," << j << " " << sgm[i][j]<< std::endl;
      std::cout << "eps " << i << "," << j << " " << eps[i][j]<< std::endl;
    }
  }
  //*/
}

#ifdef NO_PBC
#define PBC(r,L,Lh)
#else
#define PBC(r,L,Lh)\
if(r[0] <  -Lh[0]) r[0] += L[0];\
if(r[0] >=  Lh[0]) r[0] -= L[0];\
if(r[1] <  -Lh[1]) r[1] += L[1];\
if(r[1] >=  Lh[1]) r[1] -= L[1];\
if(r[2] <  -Lh[2]) r[2] += L[2];\
if(r[2] >=  Lh[2]) r[2] -= L[2];
#endif


void ForceCalculator::LJ
(Molecule *mlcl,MolTypeList mtype,Property &pry,const int nmol,const dvec3 L,double *forcepot){
  const double Lh[3] = {0.5*L[0],0.5*L[1],0.5*L[2]};
  const double rc2 = rcut*rcut;

  //roop start
#ifdef DEBUG
  pry.vdw = 0.0;
#else
#pragma omp parallel for
#endif
  for(int i=0;i<nmol;i++){
    const Molecule mlcl_i = mlcl[i];
    const MolType mt_i = mtype[mlcl_i.type];
    for(unsigned int di=0;di<mt_i.a.size();di++){
      const dvec3 a_i = body_to_space(mlcl_i.q,mt_i.a[di].r);
      for(int j=0;j<nmol;j++){
	if(i==j) continue;
	const Molecule mlcl_j = mlcl[j];
	const MolType mt_j = mtype[mlcl_j.type];
	for(unsigned int dj=0;dj<mt_j.a.size();dj++){
	  const double eps = sqrt(mt_i.a[di].e * mt_j.a[dj].e);
	  if(eps == 0.0) continue;
	  const dvec3 a_j = body_to_space(mlcl_j.q,mt_j.a[dj].r);
	  //if(i==0 && di==0) std::cout << mlcl_i.r + a_i << " " << 3*j+dj << " " << mlcl_j.r + a_j << std::endl;
	  dvec3 r = (mlcl_i.r - mlcl_j.r)*L;
	  PBC(r,L,Lh);
	  r += a_i - a_j;
	  const double r02   = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
	  if(r02>rc2) continue;
	  const double sgm = 0.5*(mt_i.a[di].s + mt_j.a[dj].s);
	  const double r02i  = 1.0/r02;
	  const double rs02i = sgm*sgm*r02i;
	  const double rs06i = rs02i*rs02i*rs02i;
	  const dvec3 f = r * 24.0*eps*rs06i*(2.*rs06i - 1.)*r02i;
	  forcepot[4*(mlcl_i.id+di)+0] += f[0];
	  forcepot[4*(mlcl_i.id+di)+1] += f[1];
	  forcepot[4*(mlcl_i.id+di)+2] += f[2];
	  pry.vir += f*r;
	  const double u = 4.0*eps*rs06i*(rs06i - 1.0);
	  forcepot[4*(mlcl_i.id+di)+3] += 0.5*u;
	  //if(i==0 && di==0 && j==1) std::cout << r << " " <<  forcepot[4*(mlcl_i.id+di)+3] << std::endl;
	  //if(i==0 && di==0) std::cout << j << " " << r << " " <<  u << std::endl;
#ifdef DEBUG
	  pry.vdw += 0.5*u;
#endif
	}
      }
      //if(i==0 && di==0) std::cout << forcepot[4*(mlcl_i.id+di)+3] << std::endl;;
    }
  }

}


void ForceCalculator::Ewald
(const int nmol,const int natom, Molecule* mlcl,MolTypeList mtype,Property &pry, dvec3 L,double* forcepot){
  EwaldDirect(nmol,natom,mlcl,mtype,pry,L,forcepot);
  EwaldWave(nmol,natom,mlcl,mtype,pry,L,forcepot);
  //EwaldSelf(nmol,natom,mlcl,mtype,pry,forcepot);
  //EwaldIntra(nmol,natom,mlcl,mtype,pry,L,forcepot);
}

void ForceCalculator::EwaldDirect
(const int nmol,const int natom, Molecule* mlcl,MolTypeList mtype,Property &pry, dvec3 L, double* forcepot){
  //double etmp = 0.e0;
  //double ptmp;
  //double psum[6] = {0.e0, 0.e0, 0.e0, 0.e0, 0.e0, 0.e0};
  const double coef = 2.0/sqrt(M_PI);
  const double alpha2= alpha*alpha;
  const dvec3 Li = L.inv();
  const dvec3 Lh = L*0.5;

#ifdef DEBUG
  pry.drct = 0.0;
#else
#pragma omp parallel for
#endif
  for(int i=0;i<nmol;i++){
    const Molecule m_i = mlcl[i];
    const MolType mt_i = mtype[m_i.type];
    for(unsigned int di=0;di<mt_i.a.size();di++){
      const dvec3 r_i = m_i.r*L + body_to_space(m_i.q,mt_i.a[di].r);
      const double aqi = mt_i.a[di].q * alpha;
      for(int j=0;j<nmol;j++){  // without newton's 3rd raw
	if(i==j) continue;
	const Molecule m_j = mlcl[j];
	const MolType mt_j = mtype[m_j.type];
	for(unsigned int dj=0;dj<mt_j.a.size();dj++){
	  dvec3  rvec= r_i - (m_j.r*L + body_to_space(m_j.q,mt_j.a[dj].r));
	  PBC(rvec,L,Lh);
	  const double ar02= sum(rvec*rvec) * alpha2;
	  const double ar  = sqrt(ar02);
	  if(ar > rcut*alpha) continue;
	  const double ari = 1.0 / ar;
	  const double aqq  = aqi * mt_j.a[dj].q;
	  const double eari= erfc(ar)*ari;
	  const dvec3  df  = rvec * alpha2*aqq*(eari + coef*exp(-ar02))*ari*ari;
	  forcepot[4*(m_i.id + di)+0] += df[0];
	  forcepot[4*(m_i.id + di)+1] += df[1];
	  forcepot[4*(m_i.id + di)+2] += df[2];
	  forcepot[4*(m_i.id + di)+3] += 0.5*aqq*eari;
#ifdef DEBUG
	  pry.drct += 0.5*aqq*eari;
#endif
	}
      }
    }
  }
}

void ForceCalculator::EwaldWave
(const int nmol,const int natom, Molecule* mlcl,MolTypeList mtype,Property &pry,const dvec3 L, double* forcepot){
  double tmp3[3];
  const double p2 = 2.e0 * M_PI;
  const double pa2 = (M_PI/alpha)*(M_PI/alpha);
  //double etmp = 0.e0;
  //double ptmp;
  //double psum[6] = {0.e0, 0.e0, 0.e0, 0.e0, 0.e0, 0.e0};
  const dvec3 Li = L.inv();
  const dvec3 Lh = L*0.5;

  const double coef1 = 2.e0 * Li[0] * Li[1] * Li[2];
  const double coef2 = 1.0 / (2.0 * M_PI * L[0] * L[1] * L[2]);
#ifdef DEBUG
  pry.wave = 0.0;
#endif

  Time f("forward dft"), b("backward dft");
  //wave part
  for (int i = 0; i < nwave; i++){
    const dvec3 vec = (dvec3)kvec[i] * Li;
    double recvec0 = 0.e0;
    double recvec1 = 0.e0;

    //forward dft
    f.beg();
#pragma omp parallel for default(shared) reduction (+:recvec0,recvec1)
    for (unsigned int j = 0; j < nmol; j++){
      const Molecule m = mlcl[j];
      const MolType mt = mtype[m.type];
      for(unsigned int d = 0;d < mt.a.size(); d++){
	const double mr = p2 * sum(( m.r*L + body_to_space(m.q,mt.a[d].r) ) * vec);
	const double q = mt.a[d].q;
	recvec0 += q * sin(mr);
	recvec1 += q * cos(mr);
      }
    }
    f.end();
    tmp3[0]  = sum(vec*vec);
    tmp3[1]  = exp(-pa2*tmp3[0]) / tmp3[0];
    tmp3[2]  = tmp3[1] * (recvec0*recvec0+recvec1*recvec1);
    recvec0 *= tmp3[1];
    recvec1 *= tmp3[1];
    /*
    etmp       = etmp + tmp3[2];
    ptmp       = 2.e0 * ((1.e0 / tmp3[0]) + pa2);
    psum[0]   += tmp3[2] * (1.e0 - ptmp * vec[0] * vec[0]);
    psum[1]   += tmp3[2] * (- ptmp * vec[0] * vec[1]);
    psum[2]   += tmp3[2] * (- ptmp * vec[0] * vec[2]);
    psum[3]   += tmp3[2] * (1.e0 - ptmp * vec[1] * vec[1]);
    psum[4]   += tmp3[2] * (- ptmp * vec[1] * vec[2]);
    psum[5]   += tmp3[2] * (1.e0 - ptmp * vec[2] * vec[2]);
    //*/

    //inverse dft
    b.beg();
#ifdef DEBUG
#else
#pragma omp parallel for
#endif
    for (int j = 0; j < nmol; j++){
      const Molecule m = mlcl[j];
      const MolType mt = mtype[m.type];
      for(unsigned int d = 0;d < mt.a.size(); d++){
	const double mr = p2 * sum( (m.r*L + body_to_space(m.q,mt.a[d].r))*vec );
	const double q = mt.a[d].q;
	const double sinmr = sin(mr);
	const double cosmr = cos(mr);
	const double sincos = recvec1 * sinmr - recvec0 * cosmr;
	forcepot[4*(m.id + d)+0] += q * coef1 * sincos * vec[0];
	forcepot[4*(m.id + d)+1] += q * coef1 * sincos * vec[1];
	forcepot[4*(m.id + d)+2] += q * coef1 * sincos * vec[2];
	forcepot[4*(m.id + d)+3] += q * coef2 * (recvec1 * cosmr + recvec0 * sinmr);
#ifdef DEBUG
	pry.wave += q * coef2 * (recvec1 * cosmr + recvec0 * sinmr);
#endif
      }
    }
    b.end();
  }
  f.print();b.print();
  /*
  etmp *= tmp3[1] * 0.5e0;
  stress[0] += psum[0] * tmp3[1] * 0.5e0 * recip[0] * recip[4] * recip[8];
  stress[1] += psum[1] * tmp3[1] * 0.5e0 * recip[0] * recip[4] * recip[8];
  stress[2] += psum[2] * tmp3[1] * 0.5e0 * recip[0] * recip[4] * recip[8];
  stress[3] += psum[1] * tmp3[1] * 0.5e0 * recip[0] * recip[4] * recip[8];
  stress[4] += psum[3] * tmp3[1] * 0.5e0 * recip[0] * recip[4] * recip[8];
  stress[5] += psum[4] * tmp3[1] * 0.5e0 * recip[0] * recip[4] * recip[8];
  stress[6] += psum[2] * tmp3[1] * 0.5e0 * recip[0] * recip[4] * recip[8];
  stress[7] += psum[4] * tmp3[1] * 0.5e0 * recip[0] * recip[4] * recip[8];
  stress[8] += psum[5] * tmp3[1] * 0.5e0 * recip[0] * recip[4] * recip[8];
  //*/
}

void ForceCalculator::EwaldSelf
(const int nmol,const int natom, Molecule* mlcl,MolTypeList mtype,Property &pry, double* forcepot){
  //self part
  const double coef = alpha / sqrt(M_PI);
  double tmp = 0.0;
  for(int i=0;i<nmol;i++){
    MolType mt = mtype[mlcl[i].type];
    for(unsigned int d=0;d<mt.a.size();d++){
      tmp -= coef * mt.a[d].q * mt.a[d].q;
    }
  }
#ifdef DEBUG
  pry.self = tmp;
#endif
}

void ForceCalculator::EwaldIntra
(const int nmol,const int natom, Molecule* mlcl,MolTypeList mtype,Property &pry,const dvec3 L, double* forcepot){
  double tmp= 0.0;
  for(int j=0;j<nmol;j++){
    const Molecule m = mlcl[j];
    const MolType mt = mtype[m.type];
    for(unsigned int di=0;di<mt.a.size();di++){
      for(unsigned int dj=0;dj<mt.a.size();dj++){
	if(dj!=di){
	  const dvec3 rij = body_to_space(m.q,mt.a[di].r) - body_to_space(m.q,mt.a[dj].r);
	  const double r = sqrt(sum(rij*rij));
	  tmp -= mt.a[di].q*mt.a[dj].q*erf(alpha*r)/r;
	  //std::cout << "i,j = " << i << j << tmp << std::endl;
	}
      }
    }
  }
  tmp *= 0.5;
#ifdef DEBUG
  pry.intra = tmp;
#endif
}

/*
void ForceCalculator::CoulombDirect
(Molecule *mlcl,const long nmol,const dvec3 L){
  const int natom = 3;

  dvec3 ftmp[natom*nmol];
  for(int i=0;i<natom*nmol;i++){
    ftmp[i] = 0.;
  }

  SPCE spce;
  const dvec3 nmax = 2;

  for(int i=0;i<nmol;i++){
    const Molecule m_i = mlcl[i];
    const dvec3 r_i = m_i.r*L;
    for(int ia=0;ia<natom;ia++){
      const dvec3 a_i = body_to_space(m_i.q,spce.r[ia]) + r_i*L;
      for(int j=0;j<nmol;j++){
	const Molecule m_j = mlcl[j];
	for(int ja=0;ja<natom;ja++){
	  const dvec3 a_j = body_to_space(m_j.q,spce.r[ia]) + m_j.r*L;

	  for(int nx=-nmax[0];nx<nmax[0]+1;nx++){
	  for(int ny=-nmax[1];ny<nmax[1]+1;ny++){
	  for(int nz=-nmax[2];nz<nmax[2]+1;nz++){
	    if((i==j)&&(nx==0 && ny==0 && nz==0)) continue;
	    dvec3 n;
	    n[0] = nx; n[1] = ny; n[2] = nz;

	    const dvec3 r = a_i - (a_j + n*L);
	    const double r02 = scalar_prod(r,r);
	    const double r03 = r02*sqrt(r02);
	    const double r03i = 1./r03;

	    const dvec3 f = r * spce.q[ia]*spce.q[ja]*r03i;
	    ftmp[natom*i + ia] += f;
	  }}}
	}
      }
    }
  }

  for(int i=0;i<nmol;i++){
    for(int m=0;m<natom;m++){
      mlcl[i].f += ftmp[natom*i+m];
      const dvec3 fbody = space_to_body(mlcl[i].q,ftmp[natom*i+m]);
      mlcl[i].n[0] += scalar_prod(spce.r[m],fbody);
      mlcl[i].n[1] += spce.r[m][1]*fbody[2] - spce.r[m][2]*fbody[1];
      mlcl[i].n[2] += spce.r[m][2]*fbody[0] - spce.r[m][0]*fbody[2];
      mlcl[i].n[3] += spce.r[m][0]*fbody[1] - spce.r[m][1]*fbody[0];
    }
    std::cout << i << " " << mlcl[i].f << std::endl;
  }
}
//*/
void ForceCalculator::SR(const Molecule *m,Atom *a,const dvec3 L,const int natom){
  const dvec3 Lh = L*0.5;

  static const double rc2 = rcut*rcut;
  static const double alpha2 = alpha*alpha;
  static const double alphai = 1.0 / alpha;
  static const double coef = 2.0/sqrt(M_PI);

  Time SR("SR");
  SR.beg();
  #pragma omp parallel for num_threads(nthreads)
  for(int i=0;i<natom;i++){
    Atom  a_i = a[i];
    const dvec3 g_i = m[a_i.i].r;
    const double aq_i = alpha*q[a_i.t];
    dvec3  f_i = 0.0, v_i = 0.0;
    double e_i = 0.0;
    for(int j=0;j<natom;j++){
      const Atom  a_j = a[j];
      if(a_i.i == a_j.i) continue; //skip same molecule

      const dvec3 rg = g_i - m[a_j.i].r;
      dvec3 r = a_i.r - a_j.r;
#if 0
      // -0.5 <= r < 0.5
      if(rg[0] >= 0.5) r[0] -= L[0];
      if(rg[1] >= 0.5) r[1] -= L[1];
      if(rg[2] >= 0.5) r[2] -= L[2];
      if(rg[0] < -0.5) r[0] += L[0];
      if(rg[1] < -0.5) r[1] += L[1];
      if(rg[2] < -0.5) r[2] += L[2];
#else
      if(r[0] >= 0.5*L[0]) r[0] -= L[0];
      if(r[1] >= 0.5*L[1]) r[1] -= L[1];
      if(r[2] >= 0.5*L[2]) r[2] -= L[2];
      if(r[0] < -0.5*L[0]) r[0] += L[0];
      if(r[1] < -0.5*L[1]) r[1] += L[1];
      if(r[2] < -0.5*L[2]) r[2] += L[2];
#endif
      const double r02 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
      if(r02>rc2) continue;
      const double r01   = sqrt(r02);
      const double r01i  = 1.0 / r01;
      //Ewald part
      const double ar  = alpha * r01;
      const double ari = alphai * r01i;
      const double aqq = aq_i * q[a_j.t];
      const double eari= erfc(ar)*ari;
      double f = alpha2*aqq*(eari + coef*exp(-ar*ar))*ari*ari;
      double e = aqq*eari;
      //LJ part
      if(eps[a_i.t][a_j.t]!=0.0){
	const double r02i  = r01i*r01i;
	const double rs02i = sgm[a_i.t][a_j.t]*r02i;
	const double rs06i = rs02i*rs02i*rs02i;
	f += 24.0 * eps[a_i.t][a_j.t] * rs06i * (2.*rs06i - 1.0) * r02i;
	e +=  4.0 * eps[a_i.t][a_j.t] * rs06i * (rs06i - 1.0);
      }
      f_i += r*f;
      e_i += e;
      v_i += r*r*f;
    }
    a[i].f += f_i;
    a[i].e += e_i*0.5;
    a[i].v += v_i*0.5;
  }
  SR.end();
  //SR.print();
}

void ForceCalculator::LR(Atom *atom,const dvec3 L,const int natom){
  static const double p2 = 2.e0 * M_PI;
  static const double pa2 = (M_PI/alpha)*(M_PI/alpha);
  //double etmp = 0.e0;
  //double ptmp;
  //double psum[6] = {0.e0, 0.e0, 0.e0, 0.e0, 0.e0, 0.e0};
  const dvec3 Li = L.inv();
  const double coef1 = 2.e0 * Li[0] * Li[1] * Li[2];
  const double coef2 = 1.0 / (2.0 * M_PI * L[0] * L[1] * L[2]);

  Time f("forward dft"), b("backward dft"), LR("LR");
  //wave part
  LR.beg();
  for (int i = 0; i < nwave; i++){
    const dvec3 vec = (dvec3)kvec[i] * Li;
    double recvec0 = 0.e0;
    double recvec1 = 0.e0;

    //forward dft
    f.beg();
#pragma omp parallel for reduction (+:recvec0,recvec1) num_threads(nthreads)
    for (unsigned int j = 0; j < natom; j++){
      const Atom a = atom[j];
      const double mr = p2 * sum( a.r * vec);
      recvec0 += q[a.t] * sin(mr);
      recvec1 += q[a.t] * cos(mr);
    }
    f.end();
    const dvec3  vecvec = vec*vec;
    const double vec02  = sum(vec*vec);
    const double vec02i = 1.0/vec02;
    const double coef3 = exp(-pa2*vec02) * vec02i;
    const double coef4 = coef2 * coef3 * (recvec0*recvec0+recvec1*recvec1);
    recvec0 *= coef3;
    recvec1 *= coef3;

    double ptmp       = 2.e0 * ((1.e0 * vec02i) + pa2);
    // virial term in wave space is all added to atom[0]
    atom[0].v[0] += coef4 * (1.e0 - ptmp * vecvec[0]);
    atom[0].v[1] += coef4 * (1.e0 - ptmp * vecvec[1]);
    atom[0].v[2] += coef4 * (1.e0 - ptmp * vecvec[2]);

    //inverse dft
    b.beg();
#pragma omp parallel for num_threads(nthreads)
    for (int j = 0; j < natom; j++){
      Atom a = atom[j];
      const double mr = p2 * sum(a.r*vec);
      const double qj     = q[a.t];
      const double sinmr  = sin(mr);
      const double cosmr  = cos(mr);
      const double sincos = recvec1 * sinmr - recvec0 * cosmr;
      atom[j].f += vec * qj * coef1 * sincos;
      atom[j].e += qj * coef2 * (recvec1 * cosmr + recvec0 * sinmr);
    }
    b.end();
  }
  LR.end();
  //f.print();b.print();LR.print();
}

void ForceCalculator::Switching(Molecule *m,Atom *a,const MolTypeList mtl,const dvec3 L,const int nmol,Property& prop){
  const double rc = rcut;
  const double rl = rswitch;
  const double rlc = rl - rc;
  const double coef = 1.0 / (rlc*rlc*rlc*rlc*rlc);

  const double rcutSq = rcut*rcut;

  double r01i,r02i,rs06i;
  dvec3 f_lj,fclmb[3];
  prop.vir = 0.0;
  prop.lj = prop.clmb = 0.0;
  int acount_i = 0;
  for(int i=0;i<nmol;i++){
    const dvec3 gr_i = m[i].r;
    const int   na_i = mtl[m[i].type].a.size();
    m[i].f = 0.0;
    int acount_j = 0;
    for(int j=0;j<nmol;j++){
      const int na_j = mtl[m[i].type].a.size();
      if(i == j){acount_j += na_j;continue;} // skip index i == j
      const dvec3 gr_j = m[j].r;

      dvec3 gr = (gr_i - gr_j)*L;
      for(int d=0;d<3;d++){
	if(gr[d] >= 0.5*L[d]) gr[d] -= L[d];
	if(gr[d] < -0.5*L[d]) gr[d] += L[d];
      }

      double r02 = gr[0]*gr[0] + gr[1]*gr[1] + gr[2]*gr[2];
      if(r02 > rcutSq){acount_j += na_j;continue;} // skip molecule outside the cut-off radious
      double r01 = sqrt(r02);
      double sw,dsw;
      if(r01 > rl){
	const double r01c = r01 - rc;
	const double r01l = r01 - rl;
	sw  = coef*r01c*r01c*r01c*(10.0*r01l*r01l - 5.0*r01l*r01c + r01c*r01c);
	dsw = coef*30.0*r01c*r01c*r01l*r01l;
      }else{
	sw  = 1.0;
	dsw = 0.0;
      }
      const dvec3 dswr = gr*(dsw/r01);

      for(int ii=0;ii<na_i;ii++){
	const Atom ai = a[acount_i+ii];
	double e=0.0;
	dvec3 f = 0.0;
	for(int jj=0;jj<na_j;jj++){
	  const Atom aj = a[acount_j+jj];
	  const double qq = q[ai.t] * q[aj.t];
	  dvec3 r = ai.r - aj.r;
	  for(int d=0;d<3;d++){
	    if(r[d] >= 0.5*L[d]) r[d] -= L[d];
	    if(r[d] < -0.5*L[d]) r[d] += L[d];
	  }
	  r02 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
	  r02i = 1.0/r02;
	  r01i = sqrt(r02i);
	  rs06i = sgm[ai.t][aj.t]*r02i;
	  rs06i = rs06i*rs06i*rs06i;

	  const double elj = 4.0 * eps[ai.t][aj.t] * rs06i * (rs06i - 1.0);
	  const double ecl = qq*r01i;
	  prop.lj   += sw*elj;
	  prop.clmb += sw*ecl;
	  const double etmp = elj + ecl;
	  double ftmp = (48.0* eps[ai.t][aj.t] * rs06i * (rs06i - 0.5) * r02i + qq*r01i*r02i)*sw;
	  e += etmp;
	  f += r*ftmp;

	  a[acount_i+ii].f += r*ftmp;
	  a[acount_i+ii].e += etmp * sw * 0.5;
	  a[acount_i+ii].v += r*r*ftmp*0.5;
	}// jj roop
	f -= dswr*e;
	m[i].f += f;
	prop.vir += f*gr;
	a[acount_i].v -= dswr*gr*e;
      }// ii_roop
      acount_j += na_j;
    }// j roop
    acount_i += na_i;
  }// i roop
  prop.vir *= 0.5;

  prop.lj *= 0.5;
  prop.clmb *= 0.5;
}

#if 1

#define ZERO    0.0
#define HALF    0.5
#define ONE     1.0
#define FIVE    5.0
#define TEN    10.0
#define TWELVE 12.0
#define THIRTY 30.0

#else

#define ZERO    0.0f
#define HALF    0.5f
#define ONE     1.0f
#define FIVE    5.0f
#define TEN    10.0f
#define TWELVE 12.0f
#define THIRTY 30.0f

#endif

#define KERNEL_LJ_CLMB(a,b)						\
  {									\
    const FP qq = q_i[a] * aq[b];					\
    FP dx = dgx + ax_i[a] - ax[b];					\
    FP dy = dgy + ay_i[a] - ay[b];					\
    FP dz = dgz + az_i[a] - az[b];					\
    r02 = dx*dx + dy*dy + dz*dz;					\
    FP r02i = ONE / r02;						\
    FP r01i = sqrt(r02i);						\
    FP sigma = (s_i[a] + as[b]);					\
    FP rs06i = sigma * sigma * r02i;					\
    rs06i = rs06i*rs06i*rs06i;						\
    const FP epsilon = e_i[a] * ae[b];					\
    const FP elj = epsilon * rs06i * (rs06i - ONE);			\
    const FP ecl = qq*r01i;						\
    glj_i += sw*elj;							\
    gcl_i += sw*ecl;							\
    const FP etmp = elj + ecl;						\
    FP ftmp = (TWELVE * epsilon * rs06i * (rs06i - HALF) * r02i + qq*r01i*r02i)*sw; \
    e += etmp;								\
    fx += dx*ftmp;							\
    fy += dy*ftmp;							\
    fz += dz*ftmp;							\
    afx_i[a] += dx*ftmp;						\
    afy_i[a] += dy*ftmp;						\
    afz_i[a] += dz*ftmp;						\
  }
#define KERNEL_LJ(a,b)							\
  {									\
    FP dx = dgx + ax_i[a] - ax[b];					\
    FP dy = dgy + ay_i[a] - ay[b];					\
    FP dz = dgz + az_i[a] - az[b];					\
    r02 = dx*dx + dy*dy + dz*dz;					\
    FP r02i = ONE/r02;							\
    FP r01i = sqrt(r02i);						\
    FP sigma = (s_i[a] + as[b]);					\
    FP rs06i = sigma * sigma * r02i;					\
    rs06i = rs06i*rs06i*rs06i;						\
    const FP epsilon = e_i[a] * ae[b];					\
    const FP etmp = epsilon * rs06i * (rs06i - ONE);			\
    glj_i += sw*etmp;							\
    FP ftmp = TWELVE* epsilon * rs06i * (rs06i - HALF) * r02i * sw;	\
    e += etmp;								\
    fx += dx*ftmp;							\
    fy += dy*ftmp;							\
    fz += dz*ftmp;							\
    afx_i[a] += dx*ftmp;						\
    afy_i[a] += dy*ftmp;						\
    afz_i[a] += dz*ftmp;						\
  }
#define KERNEL_CLMB(a,b)			\
  {						\
    const FP qq = q_i[a] * aq[b];		\
    FP dx = dgx + ax_i[a] - ax[b];		\
    FP dy = dgy + ay_i[a] - ay[b];		\
    FP dz = dgz + az_i[a] - az[b];		\
    r02 = dx*dx + dy*dy + dz*dz;		\
    FP r02i = ONE/r02;				\
    FP r01i = sqrt(r02i);			\
    const FP etmp = qq*r01i;			\
    gcl_i += sw*etmp;				\
    FP ftmp = qq*r01i*r02i*sw;			\
    e += etmp;					\
    fx += dx*ftmp;				\
    fy += dy*ftmp;				\
    fz += dz*ftmp;				\
    afx_i[a] += dx*ftmp;			\
    afy_i[a] += dy*ftmp;			\
    afz_i[a] += dz*ftmp;			\
  }

#define INNER_MOST_LOOP							\
  if(i == j) continue;							\
  FP dgx = gx_i - gx[j];						\
  FP dgy = gy_i - gy[j];						\
  FP dgz = gz_i - gz[j];						\
  /* dgx -= std::round(dgx*Lxi)*Lx; */					\
  /* dgy -= std::round(dgy*Lyi)*Ly; */					\
  dgz -= std::round(dgz*Lzi)*Lz;					\
  FP r02 = dgx*dgx + dgy*dgy + dgz*dgz;					\
  r02 = std::min(rcutSq,r02);						\
  FP r01 = sqrt(r02);							\
  FP sw = ONE, dsw = ZERO;						\
  if(r01 > rl){								\
    const FP r01c = r01 - rc;						\
    const FP r01l = r01 - rl;						\
    sw  = coef*r01c*r01c*r01c*(TEN * r01l*r01l - FIVE * r01l*r01c + r01c*r01c);	\
    dsw = coef* THIRTY *r01c*r01c*r01l*r01l;				\
  }									\
  FP r01i = ONE / r01;							\
  const FP dswx = dgx*(dsw*r01i);					\
  const FP dswy = dgy*(dsw*r01i);					\
  const FP dswz = dgz*(dsw*r01i);					\
  FP fx,fy,fz,e;							\
  if constexpr (nsite == 3){						\
    fx = fy = fz = e = ZERO;						\
    KERNEL_LJ_CLMB(0,3*j+0);						\
    KERNEL_CLMB   (0,3*j+1);						\
    KERNEL_CLMB   (0,3*j+2);						\
    fx -= dswx*e;							\
    fy -= dswy*e;							\
    fz -= dswz*e;							\
    gfx_i += fx;							\
    gfy_i += fy;							\
    gfz_i += fz;							\
    gvx_i += fx*dgx;							\
    gvy_i += fy*dgy;							\
    gvz_i += fz*dgz;							\
    fx = fy = fz = e = ZERO;						\
    KERNEL_CLMB(1,3*j+0);						\
    KERNEL_CLMB(1,3*j+1);						\
    KERNEL_CLMB(1,3*j+2);						\
    fx -= dswx*e;							\
    fy -= dswy*e;							\
    fz -= dswz*e;							\
    gfx_i += fx;							\
    gfy_i += fy;							\
    gfz_i += fz;							\
    gvx_i += fx*dgx;							\
    gvy_i += fy*dgy;							\
    gvz_i += fz*dgz;							\
    fx = fy = fz = e = ZERO;						\
    KERNEL_CLMB(2,3*j+0);						\
    KERNEL_CLMB(2,3*j+1);						\
    KERNEL_CLMB(2,3*j+2);						\
    fx -= dswx*e;							\
    fy -= dswy*e;							\
    fz -= dswz*e;							\
    gfx_i += fx;							\
    gfy_i += fy;							\
    gfz_i += fz;							\
    gvx_i += fx*dgx;							\
    gvy_i += fy*dgy;							\
    gvz_i += fz*dgz;							\
  }									\
  if constexpr (nsite == 4){						\
  fx = fy = fz = e = ZERO;						\
  KERNEL_LJ(0,j+0*offset);							\
  fx -= dswx*e;								\
  fy -= dswy*e;								\
  fz -= dswz*e;								\
  gfx_i += fx;								\
  gfy_i += fy;								\
  gfz_i += fz;								\
  gvx_i += fx*dgx;							\
  gvy_i += fy*dgy;							\
  gvz_i += fz*dgz;							\
  fx = fy = fz = e = ZERO;						\
  KERNEL_CLMB(1,j+1*offset);							\
  KERNEL_CLMB(1,j+2*offset);							\
  KERNEL_CLMB(1,j+3*offset);							\
  fx -= dswx*e;								\
  fy -= dswy*e;								\
  fz -= dswz*e;								\
  gfx_i += fx;								\
  gfy_i += fy;								\
  gfz_i += fz;								\
  gvx_i += fx*dgx;							\
  gvy_i += fy*dgy;							\
  gvz_i += fz*dgz;							\
  fx = fy = fz = e = ZERO;						\
  KERNEL_CLMB(2,j+1*offset);							\
  KERNEL_CLMB(2,j+2*offset);							\
  KERNEL_CLMB(2,j+3*offset);							\
  fx -= dswx*e;								\
  fy -= dswy*e;								\
  fz -= dswz*e;								\
  gfx_i += fx;								\
  gfy_i += fy;								\
  gfz_i += fz;								\
  gvx_i += fx*dgx;							\
  gvy_i += fy*dgy;							\
  gvz_i += fz*dgz;							\
  fx = fy = fz = e = ZERO;						\
  KERNEL_CLMB(3,j+1*offset);						\
  KERNEL_CLMB(3,j+2*offset);						\
  KERNEL_CLMB(3,j+3*offset);						\
  fx -= dswx*e;								\
  fy -= dswy*e;								\
  fz -= dswz*e;								\
  gfx_i += fx;								\
  gfy_i += fy;								\
  gfz_i += fz;								\
  gvx_i += fx*dgx;							\
  gvy_i += fy*dgy;							\
  gvz_i += fz*dgz;							\
  }									\
   if constexpr (nsite == 5){						\
     fx = fy = fz = e = ZERO;						\
     KERNEL_LJ(0,5*j+0);						\
     fx -= dswx*e;							\
     fy -= dswy*e;							\
     fz -= dswz*e;							\
     gfx_i += fx;							\
     gfy_i += fy;							\
     gfz_i += fz;							\
     gvx_i += fx*dgx;							\
     gvy_i += fy*dgy;							\
     gvz_i += fz*dgz;							\
     fx = fy = fz = e = ZERO;						\
     KERNEL_CLMB(1,5*j+1);						\
     KERNEL_CLMB(1,5*j+2);						\
     KERNEL_CLMB(1,5*j+3);						\
     KERNEL_CLMB(1,5*j+4);						\
     fx -= dswx*e;							\
     fy -= dswy*e;							\
     fz -= dswz*e;							\
     gfx_i += fx;							\
     gfy_i += fy;							\
     gfz_i += fz;							\
     gvx_i += fx*dgx;							\
     gvy_i += fy*dgy;							\
     gvz_i += fz*dgz;							\
     fx = fy = fz = e = ZERO;						\
     KERNEL_CLMB(2,5*j+1);						\
     KERNEL_CLMB(2,5*j+2);						\
     KERNEL_CLMB(2,5*j+3);						\
     KERNEL_CLMB(2,5*j+4);						\
     fx -= dswx*e;							\
     fy -= dswy*e;							\
     fz -= dswz*e;							\
     gfx_i += fx;							\
     gfy_i += fy;							\
     gfz_i += fz;							\
     gvx_i += fx*dgx;							\
     gvy_i += fy*dgy;							\
     gvz_i += fz*dgz;							\
     fx = fy = fz = e = ZERO;						\
     KERNEL_CLMB(3,5*j+1);						\
     KERNEL_CLMB(3,5*j+2);						\
     KERNEL_CLMB(3,5*j+3);						\
     KERNEL_CLMB(3,5*j+4);						\
     fx -= dswx*e;							\
     fy -= dswy*e;							\
     fz -= dswz*e;							\
     gfx_i += fx;							\
     gfy_i += fy;							\
     gfz_i += fz;							\
     gvx_i += fx*dgx;							\
     gvy_i += fy*dgy;							\
     gvz_i += fz*dgz;							\
     fx = fy = fz = e = ZERO;						\
     KERNEL_CLMB(4,5*j+1);						\
     KERNEL_CLMB(4,5*j+2);						\
     KERNEL_CLMB(4,5*j+3);						\
     KERNEL_CLMB(4,5*j+4);						\
     fx -= dswx*e;							\
     fy -= dswy*e;							\
     fz -= dswz*e;							\
     gfx_i += fx;							\
     gfy_i += fy;							\
     gfz_i += fz;							\
     gvx_i += fx*dgx;							\
     gvy_i += fy*dgy;							\
     gvz_i += fz*dgz;							\
   }

template <int nsite>
void ForceCalculator::SwitchingTuning(Molecule* m,Atom* a,const MolTypeList mtl,const dvec3 L,const int nmol,Property& prop,const int* jstart,const int* jend,const int nlane){
  const FP rc = rcut;
  const FP rl = rswitch;
  const FP rlc = rl - rc;
  const FP coef = ONE / (rlc*rlc*rlc*rlc*rlc);

  const FP rcutSq = rcut*rcut;

  const FP Lx = L[0];
  const FP Ly = L[1];
  const FP Lz = L[2];
  const FP Lxi = ONE/L[0];
  const FP Lyi = ONE/L[1];
  const FP Lzi = ONE/L[2];
#ifdef INSERT_TIMER_FORCE
  prof.beg(Profiler::Pre);
#endif

  const int offset = ((nmol + nlane - 1)/nlane)*nlane;
  if(gx  == nullptr) gx  = (FP*)_mm_malloc(nmol*sizeof(FP),64);
  if(gy  == nullptr) gy  = (FP*)_mm_malloc(nmol*sizeof(FP),64);
  if(gz  == nullptr) gz  = (FP*)_mm_malloc(nmol*sizeof(FP),64);
  if(gfx == nullptr) gfx = (FP*)_mm_malloc(nmol*sizeof(FP),64);
  if(gfy == nullptr) gfy = (FP*)_mm_malloc(nmol*sizeof(FP),64);
  if(gfz == nullptr) gfz = (FP*)_mm_malloc(nmol*sizeof(FP),64);
  if(gvx == nullptr) gvx = (FP*)_mm_malloc(nmol*sizeof(FP),64);
  if(gvy == nullptr) gvy = (FP*)_mm_malloc(nmol*sizeof(FP),64);
  if(gvz == nullptr) gvz = (FP*)_mm_malloc(nmol*sizeof(FP),64);
  if(glj == nullptr) glj = (FP*)_mm_malloc(nmol*sizeof(FP),64);
  if(gcl == nullptr) gcl = (FP*)_mm_malloc(nmol*sizeof(FP),64);
  assert(gx  != nullptr);
  assert(gy  != nullptr);
  assert(gz  != nullptr);
  assert(gfx != nullptr);
  assert(gfy != nullptr);
  assert(gfz != nullptr);
  assert(gvx != nullptr);
  assert(gvy != nullptr);
  assert(gvz != nullptr);
  assert(glj != nullptr);
  assert(gcl != nullptr);
  for(int i=0;i<nmol;i++){
    gx[i] = m[i].r[0] * L[0];
    gy[i] = m[i].r[1] * L[1];
    gz[i] = m[i].r[2] * L[2];
    gfx[i] = gfy[i] = gfz[i] = ZERO;
    gvx[i] = gvy[i] = gvz[i] = ZERO;
    glj[i] = gcl[i] = ZERO;
  }

  const int natom = nsite * offset;
  if(ax  == nullptr) ax  = (FP*)_mm_malloc(natom*sizeof(FP),64);
  if(ay  == nullptr) ay  = (FP*)_mm_malloc(natom*sizeof(FP),64);
  if(az  == nullptr) az  = (FP*)_mm_malloc(natom*sizeof(FP),64);
  if(afx == nullptr) afx = (FP*)_mm_malloc(natom*sizeof(FP),64);
  if(afy == nullptr) afy = (FP*)_mm_malloc(natom*sizeof(FP),64);
  if(afz == nullptr) afz = (FP*)_mm_malloc(natom*sizeof(FP),64);
  if(as  == nullptr) as  = (FP*)_mm_malloc(natom*sizeof(FP),64);
  if(ae  == nullptr) ae  = (FP*)_mm_malloc(natom*sizeof(FP),64);
  if(aq  == nullptr) aq  = (FP*)_mm_malloc(natom*sizeof(FP),64);
  assert(ax  != nullptr);
  assert(ay  != nullptr);
  assert(az  != nullptr);
  assert(afx != nullptr);
  assert(afy != nullptr);
  assert(afz != nullptr);
  assert(as  != nullptr);
  assert(ae  != nullptr);
  assert(aq  != nullptr);
  for(int i=0;i<nmol;i++){
    for(int d=0;d<nsite;d++){
      const int index = nsite*i + d;
      ax[i+d*offset]  = a[index].r[0] - gx[i];
      ay[i+d*offset]  = a[index].r[1] - gy[i];
      az[i+d*offset]  = a[index].r[2] - gz[i];
      ae[i+d*offset]  = sqrt(4.0*eps[a[index].t][a[index].t]);
      as[i+d*offset]  = HALF*sqrt(sgm[a[index].t][a[index].t]);
      aq[i+d*offset]  = q[a[index].t];
      afx[index] = afy[index] = afz[index] = ZERO;
    }
  }
#ifdef INSERT_TIMER_FORCE
  prof.end(Profiler::Pre);
#endif 
  __assume(nlane%(64/sizeof(FP))==0);
  __assume_aligned(gx,64);
  __assume_aligned(gy,64);
  __assume_aligned(gz,64);
  __assume_aligned(gfx,64);
  __assume_aligned(gfy,64);
  __assume_aligned(gfz,64);
  __assume_aligned(gvx,64);
  __assume_aligned(gvy,64);
  __assume_aligned(gvz,64);
  __assume_aligned(glj,64);
  __assume_aligned(gcl,64);

  __assume_aligned(ax,64);
  __assume_aligned(ay,64);
  __assume_aligned(az,64);
  __assume_aligned(afx,64);
  __assume_aligned(afy,64);
  __assume_aligned(afz,64);
  __assume_aligned(as,64);
  __assume_aligned(ae,64);
  __assume_aligned(aq,64);

#ifdef INSERT_TIMER_FORCE
  prof.beg(Profiler::Force);
#endif
  for(int iv=0;iv<(nmol/nlane)*nlane;iv+=nlane){
    const int js = jstart[iv/nlane];
    int je = jend[iv/nlane];
    if(js>=je) je += nmol;
#pragma omp simd
    for(int ii=0;ii<nlane;ii++){
      const int i = iv + ii;

      const FP gx_i = gx[i];
      const FP gy_i = gy[i];
      const FP gz_i = gz[i];
      FP gfx_i = ZERO;
      FP gfy_i = ZERO;
      FP gfz_i = ZERO;
      FP gvx_i = ZERO;
      FP gvy_i = ZERO;
      FP gvz_i = ZERO;
      FP glj_i = ZERO;
      FP gcl_i = ZERO;

      FP ax_i[nsite];
      FP ay_i[nsite];
      FP az_i[nsite];
      FP q_i[nsite];
      FP e_i[nsite];
      FP s_i[nsite];
      FP afx_i[nsite];
      FP afy_i[nsite];
      FP afz_i[nsite];
      for(int d=0;d<nsite;d++){
	ax_i[d] = ax[i + d*offset];
	ay_i[d] = ay[i + d*offset];
	az_i[d] = az[i + d*offset];
	q_i[d] =  aq[i + d*offset];
	e_i[d] =  ae[i + d*offset];
	s_i[d] =  as[i + d*offset];
	afx_i[d] = afy_i[d] = afz_i[d] = ZERO;
      }

      for(int jj = js; jj<=je; jj++){
	int j = jj % nmol;
	INNER_MOST_LOOP;
      }
      for(int d=0;d<nsite;d++){
	afx[i+d*offset] = afx_i[d];
	afy[i+d*offset] = afy_i[d];
	afz[i+d*offset] = afz_i[d];
      }
      gfx[i] = gfx_i;
      gfy[i] = gfy_i;
      gfz[i] = gfz_i;
      gvx[i] = gvx_i;
      gvy[i] = gvy_i;
      gvz[i] = gvz_i;
      gcl[i] = gcl_i;
      glj[i] = glj_i;
    } // ii loop
  }// iv loop
  // tail loop
  for(int i=(nmol/nlane)*nlane;i<nmol;i++){
    const int js = jstart[i/nlane];
    int je = jend[i/nlane];
    if(js>=je) je += nmol;

    const FP gx_i = gx[i];
    const FP gy_i = gy[i];
    const FP gz_i = gz[i];
    FP gfx_i = ZERO;
    FP gfy_i = ZERO;
    FP gfz_i = ZERO;
    FP gvx_i = ZERO;
    FP gvy_i = ZERO;
    FP gvz_i = ZERO;
    FP glj_i = ZERO;
    FP gcl_i = ZERO;

    FP ax_i[nsite];
    FP ay_i[nsite];
    FP az_i[nsite];
    FP q_i[nsite];
    FP e_i[nsite];
    FP s_i[nsite];
    FP afx_i[nsite];
    FP afy_i[nsite];
    FP afz_i[nsite];
    for(int d=0;d<nsite;d++){
      ax_i[d] = ax[i + d*offset];
      ay_i[d] = ay[i + d*offset];
      az_i[d] = az[i + d*offset];
      q_i[d] =  aq[i + d*offset];
      e_i[d] =  ae[i + d*offset];
      s_i[d] =  as[i + d*offset];
      afx_i[d] = afy_i[d] = afz_i[d] = ZERO;
    }
    for(int jj=js; jj<=je; jj++){
      const int j = jj%nmol;
      INNER_MOST_LOOP;
    }// j loop
    for(int d=0;d<nsite;d++){
      afx[i+d*offset] = afx_i[d];
      afy[i+d*offset] = afy_i[d];
      afz[i+d*offset] = afz_i[d];
    }
    gfx[i] = gfx_i;
    gfy[i] = gfy_i;
    gfz[i] = gfz_i;
    gvx[i] = gvx_i;
    gvy[i] = gvy_i;
    gvz[i] = gvz_i;
    gcl[i] = gcl_i;
    glj[i] = glj_i;
  }// i loop
#ifdef INSERT_TIMER_FORCE
  prof.end(Profiler::Force);
#endif

#ifdef INSERT_TIMER_FORCE
  prof.beg(Profiler::Post);
#endif
  for(int i=0;i<nmol;i++){
    m[i].f[0] = gfx[i];
    m[i].f[1] = gfy[i];
    m[i].f[2] = gfz[i];
  }
  prop.vir = 0.0;
  prop.lj = prop.clmb = 0.0;
  for(int i=0;i<nmol;i++){
    prop.lj   += glj[i];
    prop.clmb += gcl[i];

    prop.vir[0] += gvx[i];
    prop.vir[1] += gvy[i];
    prop.vir[2] += gvz[i];
  }
  prop.lj   *= 0.5;
  prop.clmb *= 0.5;
  prop.vir  *= 0.5;

  for(int i=0;i<nmol;i++){
    for(int d=0;d<nsite;d++){
      a[nsite*i+d].f[0] = afx[i+d*offset];
      a[nsite*i+d].f[1] = afy[i+d*offset];
      a[nsite*i+d].f[2] = afz[i+d*offset];
    }
  }
#ifdef INSERT_TIMER_FORCE
  prof.end(Profiler::Post);
#endif

}
#undef KERNEL_LJ_CLMB
#undef KERNEL_CLMB

void ForceCalculator::Switching3site(Molecule *m,Atom *a,const MolTypeList mtl,const dvec3 L,const int nmol,Property& prop,const int* js,const int* je,const int nlane){
  SwitchingTuning<3>(m,a,mtl,L,nmol,prop,js,je,nlane);
}
void ForceCalculator::Switching4site(Molecule *m,Atom *a,const MolTypeList mtl,const dvec3 L,const int nmol,Property& prop,const int* js,const int* je,const int nlane){
  SwitchingTuning<4>(m,a,mtl,L,nmol,prop,js,je,nlane);
}
void ForceCalculator::Switching5site(Molecule *m,Atom *a,const MolTypeList mtl,const dvec3 L,const int nmol,Property& prop,const int* js,const int* je,const int nlane){
  SwitchingTuning<5>(m,a,mtl,L,nmol,prop,js,je,nlane);
}


void ForceCalculator::Confined(Molecule* m,Atom *a,const MolTypeList mtl,const dvec3 L,const int nmol,Property& prop){
  const double dr = wall_length / (double)nfwall;
  int a_count = 0;
  prop.wall = 0.0;
  for(int i=0;i<nmol;i++){
    for(int ii=0;ii<mtl[m[i].type].a.size();ii++){
      if(a[a_count].t==0){ // only for oxygen
	if(confined_dim == 1){
	  dvec3  r_i = a[a_count].r - L*0.5;
#ifdef TAKAIWA
	  double r01 = sqrt(r_i[2]*r_i[2] + r_i[1]*r_i[1]);
#else
	  double r01 = sqrt(r_i[0]*r_i[0] + r_i[1]*r_i[1]);
#endif
	  int index = (int)(r01/dr);
	  if(index>=nfwall) index = nfwall -1;
	  double f = (fwall[index+1][0]-fwall[index][0])/dr*(r01 - dr*index) + fwall[index][0];
	  double e = (fwall[index+1][1]-fwall[index][1])/dr*(r01 - dr*index) + fwall[index][1];
#ifdef TAKAIWA
	  a[a_count].f[2] += f*r_i[2];
	  a[a_count].f[1] += f*r_i[1];
#ifdef SWITCHING
	  m[i].f[2] += f*r_i[2];
	  m[i].f[1] += f*r_i[1];
#endif
#else
	  a[a_count].f[0] += f*r_i[0];
	  a[a_count].f[1] += f*r_i[1];
#ifdef SWITCHING
	  m[i].f[0] += f*r_i[0];
	  m[i].f[1] += f*r_i[1];
#endif
#endif
	  a[a_count].e    += e;
	  prop.wall += e;
	}else if(confined_dim == 2){
	  int index = abs((int)( (a[a_count].r[2] - 0.5*L[2]) / dr));
	  if(index>=nfwall) index = nfwall -1;
	  a[a_count].f[2] += ((fwall[index+1][0]-fwall[index][0])/dr*(a[a_count].r[2] - dr*index) + fwall[index][0])*a[i].r[2];
	  a[a_count].e    += (fwall[index+1][1]-fwall[index][1])/dr*(a[a_count].r[2] - dr*index) + fwall[index][1];
	}
      }
      a_count++;
    }
  }
}

