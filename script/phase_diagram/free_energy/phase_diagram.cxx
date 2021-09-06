#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <string>
#include <utility>

#include "dos.hxx"

using namespace std;

int main(int argc,char** argv){
  if(argc != 5){
    cerr << "usage: " << argv[0] << " [# of molecules] [phys value file] [dos file] [output file]" << endl;
    exit(EXIT_FAILURE);
  }
  const int nmol = atoi(argv[1]);
  const string phvfile = argv[2];
  const string dosfile = argv[3];
  const string output  = argv[4];

  vector<pair<double,double> > coex;

  ifstream phv(phvfile);
  assert(!phv.fail());
  pair<double,double> eq;
  double lastc = 0.0;
  double lastp = 0.0;
  for(string line;std::getline(phv,line);){
    if(line[0]=='#') continue;
    if(line.empty()) continue;

    stringstream strs(line);
    double p,t,c,dummy;
    strs >> p >> t >> dummy >> c;
    if(lastp < p){
      if(lastp != 0) coex.push_back(eq);
      lastc = 0.0;
    }
    if(lastc < c){
      lastc = c;
      eq.first  = p;
      eq.second = t;
    }
    lastp = p;
  }

  ofstream pd(output);
  assert(!pd.fail());
  
  DoS *dos;
  dos = new DoS(nmol);
  dos->Read(dosfile);
  for(auto &it : coex){
    dos->CalcFreeEnergy(it.second,it.first);

    //for(int i=0;i<10;i++) dos->Smoothing();
    //dos->SearchSaddlePoint();
    dos->SearchLocalMinimum();
    dos->ShortestPath();
    
    stringstream strs;
    strs << 'P' << it.first << 'T' << it.second;
    const string fesfile = strs.str() + ".fes";
    ofstream fes(fesfile);
    dos->Output(fes);

    cout << "# of minima is " << dos->minima.size() << endl;
    
    const string crsfile = strs.str() + ".crs";
    ofstream crs(crsfile);
    dos->CrossSection();
    crs << "# 1:Pressure 2:Coordinate 3:RelativeFreeEnergy 4:Volume 5:Potential 6:FreeEnegy" << endl;
    for(auto &it_crs : dos->cross_section){
      const Vec3   v0 = dos->minima[0] - dos->minima[1];
      const double r0 = sqrt(v0.x*v0.x + v0.y*v0.y);
      const Vec3 diff = it_crs - dos->cross_section[0];
      const double r = sqrt(diff.x*diff.x + diff.y*diff.y) / r0;
      crs << it.first << ' ' << r << ' ' << it_crs.z - dos->gmin.z << it_crs << endl;
    }

    pd << it.first << ' ' << it.second;
    for(auto it_ex = dos->extrema.begin()+1;it_ex<dos->extrema.end();it_ex++) pd << ' ' << it_ex->z - dos->extrema[0].z;
    pd << endl;
  }
}
