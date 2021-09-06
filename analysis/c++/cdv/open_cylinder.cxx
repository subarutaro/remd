#include "open_cylinder.h"

double distance(double r0[3],double r1[3]){
  const double dx = r0[0] - r1[0];
  const double dy = r0[1] - r1[1];

  return sqrt(dx*dx + dy*dy);
}

double angle(double r[3], double center[3],const double outer){
  const double x0 = outer;
  const double y0 = 0.0;

  const double x1 = r[0] - center[0];
  const double y1 = r[1] - center[1];

  const double a = sqrt(x0*x0);
  const double b = sqrt(x1*x1 + y1*y1);
  const double c = sqrt( (x0-x1)*(x0-x1) + (y0-y1)*(y0-y1) );

  const double angle = acos( (a*a + b*b - c*c) / (2.*a*b) );
  if(y1>0) return  angle;
  else     return -angle;
}

void calc_outer_inner(vector<CDV> cdv, double center[3], double &outer, double &inner){
  outer = 0.0;inner = 100.0;
  for(vector<CDV>::iterator it = cdv.begin();it != cdv.end(); ++it){
    const double d = distance(it->r,center);
    if(d > outer) outer = d;
    if(d < inner) inner = d;
  }
}

void define_outer_inner(vector<CDV> &cdv, double center[3],const double outer, const double inner){
  for(vector<CDV>::iterator it = cdv.begin();it != cdv.end(); ++it){
    if(it->t != 0) continue;
    if(distance(it->r,center) < (outer+inner)*0.5)
      for(int i=0;i<=2;i++) (it+i)->t += 4;
  }
}

void open_cylinder(vector<CDV> &cdv, map< pair<int,int>, int> &bond, double center[3], const double outer){
  //cerr << "outer:" << outer << " inter:" << inner << endl;
  for(vector<CDV>::iterator it = cdv.begin();it != cdv.end(); ++it){
    //if(it->t > 3) continue;
    const double d = distance(it->r,center);
    const double o = angle(it->r,center,outer);
    //it->r[0] = outer * o + center[0];
    it->r[0] = d * o + center[0];
    it->r[1] = outer - d + center[1];

  }
  delete_periodic_bond(cdv,bond,center,outer);
}

void delete_periodic_bond(vector<CDV> cdv,map< pair<int,int>, int> &bond,double center[3],double outer){
  //cerr << bond.size() << endl;
  for(map<pair<int,int>, int >::iterator it = bond.begin();it != bond.end(); ++it){
    pair<int,int> index = it->first;
    double dx = cdv[index.first].r[0] - cdv[index.second].r[0];
    double dy = cdv[index.first].r[1] - cdv[index.second].r[1];
    double dz = cdv[index.first].r[2] - cdv[index.second].r[2];
    double d  = sqrt(dx*dx + dy*dy + dz*dz);
    if(d > M_PI*outer){
      //cerr << "d= " << d << ", erace " << index.first << " and " << index.second << endl;
      bond.erase(index);
    }
  }
  //cerr << bond.size() << endl;
}

