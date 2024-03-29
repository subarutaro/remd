#include "molecules.h"
#include "unit.h"

#ifdef SORT_MOLECULES
#include <algorithm>
#endif

void Molecules::MakeMolTypeList()
{ // for spce
  // making MolTypeList
  MolType water;
  double m,s,e,q;
  dvec3  r;
  int id = 0;
  std::string name;

  //data from oplsaa in gromacs
  //SPCE
  r[0] = 6.517829e-2;
  r[1] = r[2] = 0.0;
  m = 15.9994e-3 / Na / unit_mass;
  s = 3.16557e-10     / unit_length;
  e = 0.650194e3 / Na / unit_energy;
  q = -0.8476 * ele   / unit_coulomb;
  name = "   OW";
  water.AddAtom(m,r,s,e,q,id++,name);

  r[0] = -5.172994e-01;
  r[1] =  8.128467e-01;
  r[2] =  0.0;
  m = 1.00800e-3 / Na / unit_mass;
  s = 0.0;
  e = 0.0;
  q = 0.4238 * ele    / unit_coulomb;
  name = "  HW1";
  water.AddAtom(m,r,s,e,q,id++,name);
  r[1] = -r[1];
  name = "  HW2";
  water.AddAtom(m,r,s,e,q,id++,name);
  water.name = "SOL  ";
  mlist.push_back(water);
  mlist[0].CalcIner();

  /*
    std::cout << "=== MolTypeList ===" << std::endl;
    std::cout << mlist[0].name << std::endl;
    std::cout << mlist[0].a[0].name << std::endl;
    std::cout << mlist[0].a[0].m << std::endl;
    std::cout << mlist[0].a[0].r << std::endl;
    std::cout << mlist[0].a[0].s << std::endl;
    std::cout << mlist[0].a[0].e << std::endl;
    std::cout << mlist[0].a[0].q << std::endl;
    std::cout << mlist[0].a[1].name << std::endl;
    std::cout << mlist[0].a[1].m << std::endl;
    std::cout << mlist[0].a[1].r << std::endl;
    std::cout << mlist[0].a[1].s << std::endl;
    std::cout << mlist[0].a[1].e << std::endl;
    std::cout << mlist[0].a[1].q << std::endl;
    std::cout << mlist[0].a[2].name << std::endl;
    std::cout << mlist[0].a[2].m << std::endl;
    std::cout << mlist[0].a[2].r << std::endl;
    std::cout << mlist[0].a[2].s << std::endl;
    std::cout << mlist[0].a[2].e << std::endl;
    std::cout << mlist[0].a[2].q << std::endl;
    std::cout << mlist[0].I << std::endl;
  //*/
}

void Molecules::MakeMolTypeList(std::string filename){
  std::ifstream ifs(filename.c_str());
  if(ifs.fail()){
    std::cerr << "error: open molecule type file failed" << std::endl;
    exit(EXIT_FAILURE);
  }
  double m,s,e,q;
  dvec3  r;
  int aid = 0,mid = 0;
  int natom;
  std::string MolName,AtomName;
  std::string mline,aline;

  while(std::getline(ifs,mline)){
    MolType mtype;
    std::istringstream mss(mline);
    mss >> MolName;
    while(MolName.size()<5) MolName += " ";
    if(mss.fail()) continue;
    if(!isalpha(MolName[0])) continue;
    mss >> natom;
    for(int i=0;i<natom;i++){
      std::getline(ifs,aline);
      std::istringstream ass(aline);
      //reading parameters from file
      ass >> AtomName;
      while(AtomName.size()<5) AtomName = " " + AtomName;
      ass >> r[0] >> r[1] >> r[2];
      ass >> m >> s >> e >> q;
      // normalize parameters
      r *= 1e-10      / unit_length;
      m *= 1e-03 / Na / unit_mass;
      s *= 1e-10      / unit_length;
      e *= 1e+03 / Na / unit_energy;
      q *= ele        / unit_coulomb;

      mtype.AddAtom(m,r,s,e,q,aid++,AtomName);
    }
    //mtype.name = MolName;
    mtype.name = "SOL  ";
    mlist.push_back(mtype);
    mlist[mid++].CalcIner();
#if 1
    std::cout << " " << mlist[mid-1].name;
    std::cout << " " << mlist[mid-1].i[0];
    std::cout << " " << mlist[mid-1].i[1];
    std::cout << " " << mlist[mid-1].i[2];
    std::cout << std::endl;
    for(int i=0;i<natom;i++){
      MolType tmp = mlist[mid-1];
      std::cout << " " << tmp.a[i].name;
      std::cout << " " << tmp.a[i].r[0];
      std::cout << " " << tmp.a[i].r[1];
      std::cout << " " << tmp.a[i].r[2];
      std::cout << " " << tmp.a[i].m;
      std::cout << " " << tmp.a[i].s;
      std::cout << " " << tmp.a[i].e;
      std::cout << " " << tmp.a[i].q;
      std::cout << std::endl;
    }
#endif
  }
}

void Molecules::SetCubicFCC(){
  long n=1;
  while(nmol > 4*n*n*n) n++;
  const int nlattice = n;
  dvec3 l;
  l = 1.0/(double)nlattice;
  //std::cout << "nlattice= " << nlattice << std::endl;
  //std::cout << "l= " << l << std::endl;

  dvec3 unit[4];
  unit[0][0] = 0.;      unit[0][1] = 0.;      unit[0][2] = 0.;
  unit[1][0] = 0.5*l[0];unit[1][1] = 0.5*l[1];unit[1][2] = 0.;
  unit[2][0] = 0.5*l[0];unit[2][1] = 0.;      unit[2][2] = 0.5*l[2];
  unit[3][0] = 0.;      unit[3][1] = 0.5*l[1];unit[3][2] = 0.5*l[2];

  int nskip = 4*n*n*n - nmol;
  const int prob = (int)((double)RAND_MAX*(double)(4*n*n*n-nmol)/(double)(4*n*n*n));
  int count = 0;
  for(int z=0;z<nlattice;z++){
  for(int y=0;y<nlattice;y++){
  for(int x=0;x<nlattice;x++){
    dvec3 tmp;
    tmp[0] = x*l[0];tmp[1] = y*l[1];tmp[2] = z*l[2];
    for(int d=0;d<4;d++){
      if(nskip>0){
	if(rand()>prob){
	  nskip--;continue;
	}
      }
      mlcl[count++].r = unit[d]+tmp;
      if(count == nmol)	return;
    }
  }}}
}

void Molecules::KillMomentum(){
  dvec3  sum = 0.;
  double m = 0.;
  for(int i=0;i<nmol;i++){
    sum += mlcl[i].v *mlcl[i].m;
    m   += mlcl[i].m;
  }
  sum /= m;
  if(param.confined == 1){
    for(int i=0;i<nmol;i++) mlcl[i].v[2] -= sum[2];
    KillAngularMomentumInTube();
  }else if(param.confined == 2){
    for(int i=0;i<nmol;i++){
      mlcl[i].v[0] -= sum[0];
      mlcl[i].v[1] -= sum[1];
    }
  }else{
    for(int i=0;i<nmol;i++) mlcl[i].v -= sum;
  }
}

void Molecules::KillAngularMomentumInTube(){
  double mom  = 0.;
  double msum = 0.;
  for(int i=0;i<nmol;i++){
    dvec3 r = (mlcl[i].r - 0.5) * L;
    dvec3 v = mlcl[i].v / (L*tst->s);
    double m = mlcl[i].m;
    mom  += m*(r[0]*v[1] - r[1]*v[0]);
    msum += m;
  }
  mom /= msum;
  for(int i=0;i<nmol;i++){
    dvec3  r = (mlcl[i].r - 0.5) * L;
    double m = mlcl[i].m;
    double r02 = r[0]*r[0] + r[1]*r[1];
    //if(r02 != 0.0){
    if(r02 > 1.0){
      mlcl[i].v[0] += mom * r[1] / r02 * L[0] * tst->s;
      mlcl[i].v[1] -= mom * r[0] / r02 * L[1] * tst->s;
    }
  }
}

// set vel -0.5 to 0.5;
void Molecules::RandomVel(){
  dvec3 sum = 0.;
  for(int i=0;i<nmol;i++){
    for(int d=0;d<3;d++){
      mlcl[i].v[d] = (double)rand() / (double)RAND_MAX;
      mlcl[i].v[d] *= L[d];
    }
  }
  KillMomentum();
  if(((mode>>CSHIFT)&MASK)==1) KillAngularMomentumInTube();
}

void Molecules::SetMassCenterToSystemCenter(){
  dvec3  tmp  = 0.0;
  double mass = 0.0;
#ifdef TAKAIWA
  for(int i=0;i<nmol;i++){
    mlcl[i].r[1] += 0.5;
    mlcl[i].r[2] += 0.5;
  }
#else
  for(int i=0;i<nmol;i++){
    tmp  += mlcl[i].r * mlcl[i].m;
    mass += mlcl[i].m;
  }
  tmp /= mass;

  std::cout << "Total mass = " << mass << std::endl;
  std::cout << "move gravity center from (" << tmp << ") to (Lx/2, Ly/2,Lz/2)" << std::endl;
  tmp[0] = 0.5 - tmp[0];
  tmp[1] = 0.5 - tmp[1];
  tmp[2] = 0.5 - tmp[2];

  for(int i=0;i<nmol;i++){
    mlcl[i].r += tmp;
    for(int d=0;d<3;d++){
      if(mlcl[i].r[d] >= 1.0) mlcl[i].r[d] -= 1.0;
      if(mlcl[i].r[d] <  0.0) mlcl[i].r[d] += 1.0;
    }
    /*
    std::cout << " " << mlcl[i].r[0] << " " << mlcl[i].r[1] << " " << mlcl[i].r[2];
    std::cout << " " << mlcl[i].q[0] << " " << mlcl[i].q[1] << " " << mlcl[i].q[2] << mlcl[i].q[3];
    std::cout << std::endl;
    //*/
  }
#endif
}

void Molecules::VelocityScaling(){
  //std::cout << "velocity scaling" << std::endl;
  const double tra = sum(TranslationalEnergy());
  if(tra==0.){
    fprintf(stderr,"warning: rotational energy is zero. velocity scaling for translational velocity scaling is not done\n");
    return;
  }
  const double coef = sqrt(1.5*nmol*T/tra);
  for(int i=0;i<nmol;i++) mlcl[i].v *= coef;
  KillMomentum();
  if(((mode>>CSHIFT)&MASK)==1) KillAngularMomentumInTube();
  tst->s=1.0;tst->Ps=0.0;bst->Pv=0.0;
}

void Molecules::ZeroAngularVel(){
  dvec4 sum = 0.;
  for(int i=0;i<nmol;i++){
    for(int d=0;d<4;d++){
      mlcl[i].p[d] = 0.f;
    }
  }
}

void Molecules::RandomAngularVel(){
  dvec4 sum = 0.;
  for(int i=0;i<nmol;i++){
    for(int d=0;d<4;d++){
      mlcl[i].p[d] = (double)rand() / (double)RAND_MAX;
    }
    sum += mlcl[i].p;
  }
  sum /= (double)nmol;
  for(int i=0;i<nmol;i++) mlcl[i].p -= sum;
}

void Molecules::AngVelocityScaling(){
  //std::cout << "velocity scaling" << std::endl;
  double rot = RotationalEnergy();
  if(rot==0.){
    fprintf(stderr,"warning: rotational energy is zero. velocity scaling for angular velocity scaling is not done\n");
    return;
  }
  const double coef = sqrt(1.5*nmol*T/rot);
  dvec4 sum = 0.;
  for(int i=0;i<nmol;i++){
    mlcl[i].p *= coef;
    sum += mlcl[i].p;
  }
  sum /= (double)nmol;
  //for(int i=0;i<nmol;i++) mlcl[i].p -= sum;
}

#include <algorithm>

void Molecules::Sort(){
  std::sort(mlcl,&mlcl[nmol-1],[](Molecule m0,Molecule m1){return m0.r[2] < m1.r[2];});
}

void Molecules::ConvertToAtoms(){
#ifdef _OPENMP
  const int thread = omp_get_thread_num();
#else
  const int thread = 0;
#endif

#if defined(FOUR_SITE) || defined(THREE_SITE) || defined(FIVES_ITE)
  const MolType mt = mlist[mlcl[0].type];
  const int nsite = mt.a.size();

  for(unsigned int i=is[thread];i<ie[thread];i++){
    const Molecule& m = mlcl[i];
    for(unsigned int d=0;d<nsite;d++){
      Atom& a = atom[nsite*i + d];
      a.t = mt.a[d].id;
      a.i = i;
      dvec3 b2s = body_to_space(m.q, mt.a[d].r);
      a.r[0] = m.r[0] * L[0] + b2s[0];
      a.r[1] = m.r[1] * L[1] + b2s[1];
      a.r[2] = m.r[2] * L[2] + b2s[2];
#ifndef SWITCHING
      if(a.r[0] >= L[0]) a.r[0] -= L[0];
      if(a.r[1] >= L[1]) a.r[1] -= L[1];
      if(a.r[2] >= L[2]) a.r[2] -= L[2];
      if(a.r[0] <  0)    a.r[0] += L[0];
      if(a.r[1] <  0)    a.r[1] += L[1];
      if(a.r[2] <  0)    a.r[2] += L[2];
#endif
      //if(m.id==0) std::cout << "(body) r of atom 0:" << mt.a[d].r << std::endl;
      //if(m.id==0) std::cout << "(space)r of atom 0:" << body_to_space(m.q, mt.a[d].r) << std::endl;
      a.f[0] = 0.0;
      a.f[1] = 0.0;
      a.f[2] = 0.0;
#ifndef SWITCHING
      a.e = 0.0;
      a.v[0] = 0.0;
      a.v[1] = 0.0;
      a.v[2] = 0.0;
#endif
     }
  }
#else
  #pragma omp single
  {
    int atom_count = 0;
#pragma omp simd
    for(unsigned int i=0;i<nmol;i++){
      const Molecule m = mlcl[i];
      const MolType mt = mlist[m.type];
      for(unsigned int d=0;d<mt.a.size();d++){
	//assert(natom > atom_count);
	atom[atom_count].t = mt.a[d].id;
	atom[atom_count].i = i;
	atom[atom_count].r = m.r * L + body_to_space(m.q, mt.a[d].r);
#if 1
	if(atom[atom_count].r[0] >= L[0]) atom[atom_count].r[0] -= L[0];
	if(atom[atom_count].r[1] >= L[1]) atom[atom_count].r[1] -= L[1];
	if(atom[atom_count].r[2] >= L[2]) atom[atom_count].r[2] -= L[2];
	if(atom[atom_count].r[0] <  0)    atom[atom_count].r[0] += L[0];
	if(atom[atom_count].r[1] <  0)    atom[atom_count].r[1] += L[1];
	if(atom[atom_count].r[2] <  0)    atom[atom_count].r[2] += L[2];
#endif
	atom[atom_count].f = 0.;
	atom[atom_count].e = 0.;
	atom[atom_count].v = 0.;
	atom_count++;
	//if(m.id==0) std::cout << "(body) r of atom 0:" << mt.a[d].r << std::endl;
	//if(m.id==0) std::cout << "(space)r of atom 0:" << body_to_space(m.q, mt.a[d].r) << std::endl;
      }
    }
  }
#endif
}

void Molecules::ConvertFromAtoms(){
#ifdef _OPENMP
  const int thread = omp_get_thread_num();
#else
  const int thread = 0;
#endif

#if defined(FOUR_SITE) || defined(THREE_SITE) || defined(FIVES_ITE)
  double etmp = 0.0;dvec3 vtmp = 0.0;
  const MolType& mt = mlist[mlcl[0].type];
  const int nsite = mt.a.size();
  for(int i=is[thread];i<ie[thread];i++){
    Molecule& m = mlcl[i];
    dvec3 f = 0.0;
    dvec4 n = 0.0;
    for(unsigned int d = 0;d<nsite;d++){
      const Atom& a = atom[4*i + d];
#ifndef SWITCHING
      f += a.f;
#endif
      const dvec3 fbody = space_to_body(m.q,a.f);
      n[0] += scalar_prod(mt.a[d].r,fbody);
      n[1] += mt.a[d].r[1]*fbody[2] - mt.a[d].r[2]*fbody[1];
      n[2] += mt.a[d].r[2]*fbody[0] - mt.a[d].r[0]*fbody[2];
      n[3] += mt.a[d].r[0]*fbody[1] - mt.a[d].r[1]*fbody[0];
#ifndef SWITCHING
      etmp += a.e;
      vtmp += a.v;
#endif
    }
#ifndef SWITCHING
    m.f = f;
#endif
    m.n = n;
  }
  //std::cout << "pot= " << tmp << std::endl;
  #pragma omp master
  {
#ifdef SWITCHING
    prop.pot = prop.lj + prop.clmb + prop.wall;
#else
    assert(false);
    prop.pot = etmp;
    prop.vir = vtmp;
#endif
  }
#else
  #pragma omp single
  {
    int atom_count = 0;
    double etmp = 0.0;dvec3 vtmp = 0.0;
    for(int i=0;i<nmol;i++){
      dvec3 f = 0.0;
      dvec4 n = 0.0;
      const Molecule m = mlcl[i];
      MolType mt = mlist[m.type];
      for(unsigned int d = 0;d<mt.a.size();d++){
	assert(natom > atom_count);
	Atom a = atom[atom_count++];
	f += a.f;
	const dvec3 fbody = space_to_body(m.q,a.f);
	n[0] += scalar_prod(mt.a[d].r,fbody);
	n[1] += mt.a[d].r[1]*fbody[2] - mt.a[d].r[2]*fbody[1];
	n[2] += mt.a[d].r[2]*fbody[0] - mt.a[d].r[0]*fbody[2];
	n[3] += mt.a[d].r[0]*fbody[1] - mt.a[d].r[1]*fbody[0];
	etmp += a.e;
	vtmp += a.v;
      }
#ifndef SWITCHING
      mlcl[i].f = f;
#endif
      mlcl[i].n = n;
    }
    //std::cout << "pot= " << tmp << std::endl;
#ifdef SWITCHING
    prop.pot = prop.lj + prop.clmb + prop.wall;
#else
    prop.pot = etmp;
    prop.vir = vtmp;
#endif
  }
#endif
}

auto compare = [](const Molecule& a,const Molecule& b) -> bool {return a.r[2] <= b.r[2];};

template <class T,class Compare>
inline void parallel_sort(T* m,const int n,Compare c){
#ifdef _OPENMP
  const int nthreads = omp_get_num_threads();
 #pragma omp single
  {
    parallel_sort_body(m,0,n-1,c,n/nthreads);
  }
#else
  std::sort(m,m+n,c);
#endif
}

template <class T,class Compare>
void parallel_sort_body(T* m,const int left,const int right,Compare c,const int threshold = 32){
#ifdef _OPENMP
  const int tid = omp_get_thread_num();
#else
  const int tid = 0;
#endif
  //printf("%d: %llu, %d %d\n",tid,m,left,right);
  fflush(stdout);

  if(left >= right) return;
  if(right - left <= threshold){
    std::sort(&m[left],&m[right+1],c);
    return;
  }
  // partitioning
  int ii = left;
  int jj = right;
  auto pivot = m[(left+right)/2];
  while(true){
    //printf("  %d: ii= %d, jj= %d\n",tid,ii,jj);
    while(c(m[ii],pivot)) ii++;
    while(c(pivot,m[jj])) jj--;
    if(ii >= jj) break;
    auto tmp = m[ii];
    m[ii] = m[jj];
    m[jj] = tmp;
    ii++;
    jj--;
  }
  //printf("    %d: ii= %d, jj= %d\n",tid,ii,jj);
  #pragma omp task
  parallel_sort_body(m,left,jj,c,threshold);
  #pragma omp task
  parallel_sort_body(m,ii,right,c,threshold);
  #pragma omp taskwait
}

void Molecules::CalcForcePot(const bool doSort,const double mergin){
#if 0//def ENABLE_AOS_TO_SOA_CONVERSION
  prof.beg(Profiler::SoAtoAoS);
#pragma omp single
  {
    SoAtoAoS();
  }
  prof.beg(Profiler::SoAtoAoS);
#endif

#ifdef SWITCHING
#if defined(FOUR_SITE) || defined(THREE_SITE) || defined(FIVE_SITE)
  prof.beg(Profiler::Sort);
#pragma omp single
  if(doSort){
    std::sort(mlcl,mlcl+nmol,compare);
  } // omp single
  prof.end(Profiler::Sort);
  prof.beg(Profiler::Neigh);
  if(doSort){
#ifdef _OPENMP
    const int thread = omp_get_thread_num();
#else
    const int thread = 0;
#endif
    const int ivs = (is[thread] / nlane) * nlane;
    const int ive = (ie[thread] == nmol) ? ((nmol+nlane-1)/nlane)*nlane : (ie[thread]/nlane)*nlane;
    for(int i=ivs;i<ive;i+=nlane){
      int s = i;
      double lb = mlcl[s].r[2] - (fclcltr->rcut + mergin)/L[2];
      while(mlcl[s].r[2] >= lb){
	s--;
	if(s<0){
	  s += nmol;
	  lb += 1.0;
	}
      }
      jstart[i/nlane] = s;

      int e = std::min(i+nlane-1,nmol-1);
      double ub = mlcl[e].r[2] + (fclcltr->rcut + mergin)/L[2];
      while(mlcl[e].r[2] <= ub){
	e++;
	if(e>=nmol){
	  e -= nmol;
	  ub -= 1.0;
	}
      }
      jend[i/nlane] = e;
      assert(0 <= s && s < nmol);
      assert(0 <= e && e < nmol);
      //printf("(js,je)= %d %d\n",s,e);
    }
  }
  prof.end(Profiler::Neigh);
  //printf("gen neighbor list finished\n");
#endif
#endif

  prof.beg(Profiler::M2A);
  ConvertToAtoms();
  prof.end(Profiler::M2A);
  //printf("M2A finished\n");

#ifdef SWITCHING
  prof.beg(Profiler::Switching);
#ifdef THREE_SITE
  fclcltr->Switching3site(mlcl,atom,mlist,L,nmol,prop,jstart,jend,nlane);
#elif FOUR_SITE
  fclcltr->Switching4site(mlcl,atom,mlist,L,nmol,prop,jstart,jend,nlane,is,ie);
#elif FIVE_SITE
  fclcltr->Switching5site(mlcl,atom,mlist,L,nmol,prop,jstart,jend,nlane);
#else
  fclcltr->Switching(mlcl,atom,mlist,L,nmol,prop);
#endif
  prof.end(Profiler::Switching);
#else
  fclcltr->SR(mlcl,atom,L,natom);
  fclcltr->LR(atom,L,natom);
#endif
  //printf("Non-bond force finished\n");
  prof.beg(Profiler::Wall);
  if(((mode>>CSHIFT)&MASK)>0) fclcltr->Confined(mlcl,atom,mlist,L,nmol,prop,is,ie);
  prof.end(Profiler::Wall);
  //printf("Wall force finished\n");
  prof.beg(Profiler::A2M);
  ConvertFromAtoms();
  prof.end(Profiler::A2M);
  //printf("A2M finished\n");

#if 0//def ENABLE_AOS_TO_SOA_CONVERSION
  prof.beg(Profiler::AoStoSoA);
  #pragma omp single
  {
    AoStoSoA();
  }
  prof.end(Profiler::AoStoSoA);
#endif
}

dvec3 Molecules::TranslationalEnergy(){
  const dvec3 coef = (L*tst->s).inv();

  dvec3 kin = 0.0;
  for(int i=0;i<nmol;i++){
    const Molecule m = mlcl[i];
    kin += (m.v*m.v)*m.m;
  }
  kin *= coef*coef*0.5;
  return kin;
}

double Molecules::RotationalEnergy(){
  const dvec3 coef = 0.25/(tst->s);
  dvec3 rot = 0.0;
  for(int i=0;i<nmol;i++){
    const Molecule m = mlcl[i];
    dvec3 xi;
    xi[0] = scalar_prod(m.p,P1(m.q));
    xi[1] = scalar_prod(m.p,P2(m.q));
    xi[2] = scalar_prod(m.p,P3(m.q));
    xi   *= m.i.inv()*coef;

    rot += m.i*xi*xi*2.0;
  }
  return sum(rot);
}

double Molecules::ThermostatEnergy(const double gkT){
  return 0.5*tst->Ps*tst->Ps/tst->Q + gkT*log(tst->s);
}

double Molecules::BarostatEnergy(const double P){
  if(((mode>>PSHIFT)&MASK) > 0){
    if(((mode>>CSHIFT)&MASK)==0)
      return 0.5*scalar_prod(bst->Pv,bst->Pv)/bst->W + P*L[0]*L[1]*L[2];
    if(((mode>>CSHIFT)&MASK)==1){
      return 0.5*bst->Pv[2]*bst->Pv[2]/bst->W + P*Molecules::GetVolume();
    }
    if(((mode>>CSHIFT)&MASK)==2){
      double wl = param.wall_length;
      return 0.5*(bst->Pv[0]*bst->Pv[0] * bst->Pv[1]*bst->Pv[1])/bst->W + P*Molecules::GetVolume();
    }
  }else{
    return 0.0;
  }
}

dvec3 Molecules::Momentum(){
  const dvec3 coef = (L*tst->s).inv();
  dvec3 tmp = 0.0;
  for(int i=0;i<nmol;i++){
    tmp += mlcl[i].v * mlcl[i].m * coef;
  }
  return tmp;
}

dvec3 Molecules::RotationalMomentum(){
  dvec3 tmp = 0.0;
  for(int i=0;i<nmol;i++){
    const Molecule m = mlcl[i];
    const dvec44 St = S_T(m.q);
    dvec3 a;
    a[0] = scalar_prod(St[1],m.p);
    a[1] = scalar_prod(St[2],m.p);
    a[2] = scalar_prod(St[3],m.p);
    tmp += a;
  }
  return tmp;
}

void Molecules::CalcProperties(){
  // pry.pot is updated by CalcForce
  //CalcForcePot();
  dvec3 tmp = TranslationalEnergy();
  prop.tra = sum(tmp);
  prop.rot = RotationalEnergy();
  prop.kin = prop.tra + prop.rot;
  prop.T   = prop.kin * unit_temp / (3.*nmol);
  prop.Tave += prop.T;

  prop.tsm = 0.5 * tst->Ps * tst->Ps / tst->Q;
  prop.tsp = prop.gkT * log(tst->s);
  prop.tst = prop.tsm + prop.tsp;

  prop.bst = BarostatEnergy(P);

  prop.tot = prop.pot + prop.kin;
  if(((mode>>TSHIFT)&MASK) == 1) prop.tot += prop.tst;
  if(((mode>>PSHIFT)&MASK) >  0) prop.tot += prop.bst;
  //prop.ham = tst->s * (prop.tot - prop.H0);
  prop.ham = prop.tot - prop.H0;

  prop.tmo = Momentum();
  prop.rmo = RotationalMomentum();

  const double volume = Molecules::GetVolume();
  prop.prs[0] = (tmp[0]*2.0 + prop.vir[0]) / volume * unit_press;
  prop.prs[1] = (tmp[1]*2.0 + prop.vir[1]) / volume * unit_press;
  prop.prs[2] = (tmp[2]*2.0 + prop.vir[2]) / volume * unit_press;
  prop.prs[3] = sum(tmp*2.0 + prop.vir) / (3.*volume) * unit_press;
  prop.Pave += prop.prs;
  prop.nave++;
}

void Molecules::PrintAll(std::ostream &s){
  s << "#Coordinate( (0 0 0) < (x,y,z) < (" << L << ") )"<< std::endl;
  for(int i=0;i<nmol;i++) s << i << " " << mlcl[i].r*L << std::endl;
  s << "#Velocity" << std::endl;
  for(int i=0;i<nmol;i++) s << i << " " << mlcl[i].v << std::endl;
  s << "#Force" << std::endl;
  for(int i=0;i<nmol;i++) s << i << " " << mlcl[i].f << std::endl;
  s << "#Angle(q0,q1,q2,q3,scalar_prod(q,q)=1.0)" << std::endl;
  for(int i=0;i<nmol;i++) s << i << " " << mlcl[i].q << " " << scalar_prod(mlcl[i].q,mlcl[i].q)<< std::endl;
  s << "#Angular velocity" << std::endl;
  for(int i=0;i<nmol;i++) s << i << " " << mlcl[i].p << std::endl;
  s << "#Torque" << std::endl;
  for(int i=0;i<nmol;i++) s << i << " " << mlcl[i].n << std::endl;
  s << "#Mass" << std::endl;
  for(int i=0;i<nmol;i++) s << i << " " << mlcl[i].m << std::endl;
  s << "#Coordinate of atoms" << std::endl;
  for(int i=0;i<natom;i++) s << i << " " << atom[i].r << std::endl;
  s << "#Force of atoms" << std::endl;
  for(int i=0;i<natom;i++) s << i << " " << atom[i].f << std::endl;

  s << "#Thermostat (s, Ps, Q)" << std::endl;
  s << tst->s << " " << tst->Ps << " " << tst->Q << std::endl;
  s << "#Barostat (Pv, W)" << std::endl;
  s << bst->Pv << " " << bst->W << std::endl;
}
