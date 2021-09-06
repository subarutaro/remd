#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <string>
#include <utility>

#include <cassert>

#include "vec.h"

const int XMAX = 4000;
const int YMAX = 4000;
const int CMAX = 1000;

class DoS{
 private:
  Vec3   value[XMAX][YMAX];
  double g[XMAX][YMAX];
  double slope[XMAX][YMAX];
  double min  = 1e30;
  double gmin = 1e30;
  int ncontour = 0;
  std::vector<Vec3> contour[CMAX];
  std::vector<std::pair<int,int> > local;
  std::vector<Vec3> cluster;
  std::vector<Vec3> border;
  int gx,gy;
  int nmol;
  int xmax,ymax;

 public:
  DoS();
  DoS(const int _nmol) : nmol(_nmol){}
  ~DoS();

  void Read(const std::string &filename);
  void CalcFreeEnergy(const double temperature,const double pressure);
  void Smoothing();
  void Contouring();
  void Clustering();
  void SearchSeddlePoint();

  void Output(std::ostream&);
};
