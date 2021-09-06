#include "seddle_point.h"

#include <algorithm>

using namespace std;

DoS::DoS(){}
DoS::~DoS(){}

void DoS::Read(const string &filename){
  ifstream ifs(filename);
  assert(!ifs.fail());

  int ix=0,iy=0;
  for(string line;getline(ifs,line);){
    if(line.empty()){
      ix++;
      ymax = iy;
      iy=0;
      continue;
    }
    if(!isdigit(line[0])) continue;

    stringstream strs(line);
    Vec3 v;
    strs >> v;
    value[ix][iy++] = v;
    if(v.z < min) min = v.z;
  }
  xmax = ix;
}

void DoS::CalcFreeEnergy
(const double temperature,
 const double pressure){
  for(int x=0;x<xmax;x++){
    for(int y=0;y<ymax;y++){
      const Vec3 v = value[x][y];
      g[x][y] = (v.z > min) ? -v.z + (v.y + v.x*pressure)*nmol/temperature : 0;
      if(g[x][y] < gmin){ gmin = g[x][y]; gx = x; gy = y;}
      //cout << v.x << ' ' << v.y << ' ' << g[x][y] << endl;
    }
    //cout << endl;
  }
}

void DoS::Smoothing(){
  double g_tmp[xmax][ymax];
  for(int x=0;x<xmax;x++)
    for(int y=0;y<ymax;y++) g_tmp[x][y] = g[x][y];
      
  for(int x=1;x<xmax-1;x++){
    for(int y=1;y<ymax-1;y++){
      double sum = g_tmp[x][y];
      double count = 1.0;
      const double g0 = g_tmp[x-1][y];
      const double g1 = g_tmp[x+1][y];
      const double g2 = g_tmp[x][y-1];
      const double g3 = g_tmp[x][y+1];

      if(g0 != 0){ sum += g0; count += 1.0;}
      if(g1 != 0){ sum += g1; count += 1.0;}
      if(g2 != 0){ sum += g2; count += 1.0;}
      if(g3 != 0){ sum += g3; count += 1.0;}

      g[x][y] = sum / count;
    }
  }
}

void DoS::SearchSeddlePoint(){
  const Vec3 k(0.0,0.0,1.0);
  for(int x=0;x<xmax-1;x++){
    for(int y=0;y<ymax-1;y++){
      const Vec3 v = value[x][y];
      const Vec3 v0 = v - value[x+1][y];
      const Vec3 v1 = v - value[x][y+1];
      const Vec3 n = cross(v0,v1);
      slope[x][y] = (g[x][y] != 0.0) ? dot(n/norm(n),k) : 0;
    }
  }
}

void DoS::Contouring(){
  const double dz = 0.2;
  int i=0;
  for(double z = gmin;z < gmin+20.;z += dz,i++){
    for(int x=0;x<xmax;x++){
      for(int y=0;y<ymax;y++){
	if(z <= g[x][y] && g[x][y] < z+dz){
	  Vec3 v = value[x][y];
	  v.z = z + 0.5*dz;
	  contour[i].push_back(v);
	}
      }
    }
  }
  ncontour = i;
}

void DoS::Clustering(){
  cout << "--- start clustering ---" << endl;
  pair<int,int> prev(0,0);
  for(int x=1;x<xmax-1;x++){
    for(int y=1;y<ymax-1;y++){
      Vec3 v = value[x][y];
      v.z = g[x][y];
      if(g[x][y] == 0) continue;
      int ix = x, iy = y;
      while(true){
	if(g[ix-1][iy] < g[ix][iy]){
	  ix--; continue;
	}else if(g[ix+1][iy] < g[ix][iy]){
	  ix++; continue;
	}else if(g[ix][iy-1] < g[ix][iy]){
	  iy--; continue;
	}else if(g[ix][iy+1] < g[ix][iy]){
	  iy++; continue;
	}else{
	    break;
	}
      }
      pair<int,int> p(ix,iy);
      auto it = find(local.begin(),local.end(),p);
      if(p != prev) border.push_back(v);
      if(it == local.end() && g[ix][iy] < 0.5*gmin){
	local.push_back(p);
      }
      prev = p;
    }
  }
  cout << "--- writing local_minima ---" << endl;
  ofstream lm("local_minima.dat");
  for(auto &it : local)
    lm << value[it.first][it.second] << ' ' << g[it.first][it.second] << endl;

  ofstream bd("border.dat");
  for(auto &it : border) bd << it << endl;
    
  
  cout << "--- writing cluster ---" << endl;
  ofstream cl("cluster.dat");
  if(cl.fail()) cerr << "error: opening cluster.dat failed ---" << endl;
  for(auto &it : local){
    for(int x=1;x<xmax-1;x++){
      for(int y=1;y<ymax-1;y++){
	Vec3 v = value[x][y];
	v.z = g[x][y];
	if(g[x][y] == 0) continue;
	int ix = x, iy = y;
	while(true){
	  if(g[ix-1][iy] < g[ix][iy]){
	    ix--; continue;
	  }else if(g[ix+1][iy] < g[ix][iy]){
	    ix++; continue;
	  }else if(g[ix][iy-1] < g[ix][iy]){
	    iy--; continue;
	  }else if(g[ix][iy+1] < g[ix][iy]){
	    iy++; continue;
	  }else{
	    pair<int,int> p(ix,iy);
	    if(it == p) cl << v << endl;
	    break;
	  }
	}
      }
    }
    cl << endl << endl;
  }
 
}

void DoS::Output(ostream &os){
  for(int x=0;x<xmax;x++){
    for(int y=0;y<ymax;y++){
      os << value[x][y] << ' ' << g[x][y] << ' ' << slope[x][y] << endl;
    }
    os << endl;
  }
  os << endl << endl;
}

#define DOS_TEST
#ifdef DOS_TEST

#include <cstdlib>
int main(int argc,char** argv){
  if(argc != 6){
    cerr << "usage: " << argv[0] << " [input file] [nmol] [temperature] [pressure] [output file]" << endl;
    exit(EXIT_FAILURE);
  }

  const string input       = argv[1];
  const int    nmol        = atoi(argv[2]);
  const double temperature = atof(argv[3]);
  const double pressure    = atof(argv[4]);
  const string output      = argv[5];

  DoS *dos;
  dos = new DoS(nmol);
  dos->Read(input);
  dos->CalcFreeEnergy(temperature,pressure);
  dos->Contouring();
  //dos->SearchSeddlePoint();
  ofstream ofs(output);
  dos->Output(ofs);
  for(int i=0;i<10;i++) dos->Smoothing();
  dos->Clustering();

}

#endif
