#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <string>
#include <utility>

#include <cassert>

#include "vec.hxx"

const int XMAX = 4000;
const int YMAX = 4000;
const int CMAX = 1000;

class DoS{
 public:
  Vec3   value[XMAX][YMAX];
  Vec3   g[XMAX][YMAX];
  Vec3   gs[XMAX][YMAX]; // smoothed g
  double slope[XMAX][YMAX];
  double dosmin  = 1e30;
  Vec3  gmin = Vec3(0,0,1e30);
  const double goffset = 10.;
  int ncontour = 0;
  std::vector<Vec3> contour[CMAX];

  std::vector<Vec3> minima;
  std::vector<Vec3> extrema;
  std::vector<Vec3> cluster;
  std::vector<Vec3> border;
  std::vector<Vec3> cross_section;
  std::vector<Vec3> shortest_path;
  class Peak{
  public:
    int x,y;
    double z;

    Peak(const int _x,const int _y,const double _z) : x(_x),y(_y),z(_z) {}
    const bool operator>(const Peak p) const {return z > p.z;}
    const bool operator<(const Peak p) const {return z < p.z;}
  };
  std::vector<Peak> peak;
  Vec3 saddle;
  int gx,gy;
  int nmol;
  int xmax,ymax;
  int xmin,ymin;

  DoS();
  DoS(const int _nmol) : nmol(_nmol){}
  ~DoS();

  void Read(const std::string &filename);
  void CalcFreeEnergy(const double temperature,const double pressure);
  void Smoothing();
  void Contouring();
  void Clustering();
  void CrossSection();
  void SearchSaddlePoint();
  void SearchLocalMinimum();
  void ShortestPath();

  void Output(std::ostream&);
};
