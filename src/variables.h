#pragma once
#include <iostream>
#include "vector.h"
#include "unit.h"

struct Property{
  double time;// simulation time (sec)

#ifdef DEBUG
  double vdw;  // van der waals potential energy

  double clm;  // total coulombic potential energy (= drct + wave + self + intra)
  double drct; // diretct part of coulombic potential energy
  double wave; // wave part of coulombic potential energy
  double self; // self part of coulombic potential energy
  double intra;// intra part of coulombic potential energy

  double wall; // potential between confined moelcules and wall
#endif
  double pot; // vdw + clmb;
  double tra; // translational kinetic energy
  double rot; // rotational kinetic energy
  double kin; // tra + rot

  double tsm; // momen?t of thermostat
  double tsp; // potential? of thermostat
  double tst; // energy of thermostat

  double bst; // energy of barostat

  double tot; // pot + kin + tst + bst
  double ham; // s * (tot - H0);

  dvec3  tmo; // tlanslational momentum
  dvec3  rmo; // rotational momentum
  dvec3  vir; // 3D virial
  dvec4  prs; // P,Pxx,Pyy,Pzz

  double Tave;
  dvec4  Pave;
  long   nave;

  double gkT;
  double H0;

  Property(){time=0.;Pave=Tave=0.0;nave=0;};
  void flush(){time=0.;Pave=Tave=0.0;nave=0;};

  friend std::ostream &operator<<(std::ostream &s, const Property &p){
    s << std::setw(10) << std::scientific << p.time << " " ;
    s << std::setw(10) << p.pot  << " " ;
    s << std::setw(10) << p.tra  << " " ;
    s << std::setw(10) << p.rot  << " " ;
    s << std::setw(10) << p.kin  << " " ;
    s << std::setw(10) << p.tsm  << " " ;
    s << std::setw(10) << p.tsp  << " " ;
    s << std::setw(10) << p.tst  << " " ;
    s << std::setw(10) << p.bst  << " " ;
    s << std::setw(10) << p.tot  << " " ;
    s << std::setw(10) << p.ham  << " " ;
    s << p.tmo << " " ;
    if(p.nave!=0){
      s << p.Tave/(double)p.nave << " " ;
    }else{
      s << p.Tave << " " ;
    }
    s << p.prs << " " ;
    s << p.vir << " " ;
    if(p.nave!=0){
      s << p.Pave/(double)p.nave << " " ;
    }else{
      s << p.Pave << " " ;
    }
    return s;
  }
};

struct Atom{
  int    t; // atom type
  int    i; // molecule index
  dvec3  r; // coordinate
  dvec3  f; // force
  double e; // energy
  dvec3  v; // virial
};

struct Molecule{
  dvec3  r;   //coordinate
  int    type;//molecule type
  dvec3  v;   //velocity
  double m;   //mass
  dvec3  f;   //force
  int    id;  //index of first atom in this moleculet
  dvec3  i;   //inertia
  dvec4  q;   //angle
  dvec4  p;   //angular velocity
  dvec4  n;   //torque

  Molecule(){
    r = v = f = 0.;
    p = q = n = 0.;
    q[0] = 1.;
    m = 0.0;
    type = 0;
  }
  ~Molecule(){};

  friend std::ostream &operator<<(std::ostream &s, const Molecule &m){
    s << "r: " << m.r << std::endl;
    s << "v: " << m.v << std::endl;
    s << "f: " << m.f << std::endl;
    s << "q: " << m.q << std::endl;
    s << "p: " << m.p << std::endl;
    s << "n: " << m.n << std::endl;
    s << "type: " << m.type << std::endl;
    s << "id:   " << m.id;
    return s;
  };
};

struct Thermostat{
  double s;
  double Ps;
  double Q;
  Thermostat(const double _Q){
    Q = _Q;
    s  = 1.0;
    Ps = 0.0;
  };
  void flush(){s=1.0;Ps=0.0;};
};

struct Barostat{
  dvec3  Pv;
  double W;
  Barostat(const double _W){
    W = _W;
    Pv = 0.;
  }
  void flush(){Pv=0.0;};
};

struct AtomType{
  double m;
  dvec3  r;
  double s;
  double e;
  double q;
  int    id;
  std::string name;
  AtomType
  (const double _m, const dvec3 _r,
   const double _s, const double _e,
   const double _q, const int _id,
   const std::string _name)
  //: m(_m), r(_r), s(_s), e(_e), q(_q), id(_id),name(_name) {};
  //*
  {
    m=_m;r=_r;s=_s;e=_e;q=_q;id=_id;name=_name;
  }
  //*/
};

struct MolType{
  std::string name;
  std::vector<AtomType> a;
  dvec3 i;
  void AddAtom
  (const double _m,const dvec3 _r,const double _s,const double _e,const double _q,const int _id,const std::string _name){
    AtomType tmp(_m,_r,_s,_e,_q,_id,_name);
    a.push_back(tmp);
  };
  void CalcIner(){
    double xx,yy,zz,xy,yz,zx;
    xx = yy = zz = xy = yz = zx = 0.0;
    for(unsigned int j=0;j<a.size();j++){
      double m = a[j].m;
      double x = a[j].r[0];double y = a[j].r[1];double z = a[j].r[2];
      xx += m*(y*y + z*z);
      yy += m*(z*z + x*x);
      zz += m*(x*x + y*y);
      xy -= m*(y*z);
      yz -= m*(z*x);
      zx -= m*(x*y);
    }
    assert(xy < 1e-10 && yz < 1e-10 && zx < 1e-10); // xy, yz, and zx must be 0
    //std::cout << "xy yz zx = " << xy << " " << yz << " " << zx << std::endl;
    i[0] = xx;i[1] = yy; i[2] = zz;
    //std::cout << "I= " << I << std::endl;
  }
};

typedef std::vector<MolType> MolTypeList;
