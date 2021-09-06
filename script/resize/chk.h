#include <iostream>
#include <vector>
#include <string>

class Molecule{
public:
  double r[3]; //coordinate
  int    type; //molecule type
  double v[3]; //velocity
  double m;    //mass
  int    id;   //index of first atom in this moleculet
  double i[3]; //inertia
  double q[4]; //angle
  double p[4]; //angular velocity

  Molecule(){}
  Molecule(const Molecule &_m){
    for(int d=0;d<3;d++) r[d] = _m.r[d];
    type = _m.type;
    for(int d=0;d<3;d++) v[d] = _m.v[d];
    m  = _m.m;
    id = _m.id;
    for(int d=0;d<3;d++) i[d] = _m.i[d];
    for(int d=0;d<4;d++) q[d] = _m.q[d];
    for(int d=0;d<4;d++) p[d] = _m.p[d];
  }
  ~Molecule(){}

  friend std::ostream &operator<<(std::ostream &os, const Molecule &m){
    os << " " << m.r[0] << " " << m.r[1] << " " << m.r[2];
    os << " " << m.v[0] << " " << m.v[1] << " " << m.v[2];
    os << " " << m.q[0] << " " << m.q[1] << " " << m.q[2] << " " << m.q[3];
    os << " " << m.p[0] << " " << m.p[1] << " " << m.p[2] << " " << m.p[3];
    os << " " << m.m;
    os << " " << m.i[0] << " " << m.i[1] << " " << m.i[2];
    os << " " << m.type;
    os << " " << m.id;
    return os;
  }
  friend std::istream &operator>>(std::istream &is, Molecule &m){
    is >> m.r[0] >> m.r[1] >> m.r[2];
    is >> m.v[0] >> m.v[1] >> m.v[2];
    is >> m.q[0] >> m.q[1] >> m.q[2] >> m.q[3];
    is >> m.p[0] >> m.p[1] >> m.p[2] >> m.p[3];
    is >> m.m;
    is >> m.i[0] >> m.i[1] >> m.i[2];
    is >> m.type;
    is >> m.id;
    return is;
  }
};

class Thermostat{
public:
  double s;
  double Ps;
  double m;
  Thermostat(){}
  Thermostat(const Thermostat &t){s=t.s;Ps=t.Ps;m=t.m;}
  Thermostat(const double _s,const double _Ps,const double _m) : s(_s),Ps(_Ps),m(_m) {}

  friend std::ostream &operator<<(std::ostream &os, const Thermostat &t){
    os << " " << t.s << " " << t.Ps << " " << t.m;
    return os;
  }
  friend std::istream &operator>>(std::istream &is, Thermostat &t){
    is >> t.s >> t.Ps >> t.m;
    return is;
  }
};

class Barostat{
public:
  double x,y,z;
  double m;
  Barostat(){}
  Barostat(const Barostat &b){x=b.x;y=b.y;z=b.z;m=b.m;}
  Barostat(const double _x,const double _y,const double _z,const double _m) : x(_x),y(_y),z(_z),m(_m) {}
  friend std::ostream &operator<<(std::ostream &os, const Barostat &b){
    os << " " << b.x << " " << b.y << " " << b.z << " " << b.m;
    return os;
  }
  friend std::istream &operator>>(std::istream &is, Barostat &b){
    is >> b.x >> b.y >> b.z >> b.m;
    return is;
  }
};

class Length{
public:
  double x,y,z;
  Length(){};
  Length(const Length &l){x=l.x;y=l.y;z=l.z;}
  Length(const double _x,const double _y,const double _z) : x(_x),y(_y),z(_z) {}

  friend std::ostream &operator<<(std::ostream &os, const Length &l){
    os << " " << l.x << " " << l.y << " " << l.z;
    return os;
  }
  friend std::istream &operator>>(std::istream &is, Length &l){
    is >> l.x >> l.y >> l.z;
    return is;
  }
};

class CheckPointFile{
private:
  int nmol,natom;
  std::vector<Molecule> m;
  Thermostat t;
  Barostat   b;
  Length     l;

public:
  void Resize(const Length);
  void Resize(const int,const int);

  friend std::ostream &operator<<(std::ostream &os, CheckPointFile &chk){
    os << " " << chk.nmol << " " << chk.natom << std::endl;
    for(auto &it : chk.m) os << it << std::endl;
    os << chk.t << std::endl;
    os << chk.b << std::endl;
    os << chk.l << std::endl;
    return os;
  }
  friend std::istream &operator>>(std::istream &is, CheckPointFile &chk){
    is >> chk.nmol >> chk.natom;
    for(int i=0;i<chk.nmol;i++){
      Molecule _m;
      is >> _m;
      chk.m.push_back(_m);
    }
    is >> chk.t;
    is >> chk.b;
    is >> chk.l;
    return is;
  }
};
