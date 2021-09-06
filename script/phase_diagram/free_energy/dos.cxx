#ifndef H_DOS
#define H_DOS

#include "dos.hxx"

#include <algorithm>

using namespace std;

DoS::DoS(){}
DoS::~DoS(){}

void DoS::Read(const string &filename){
  ifstream ifs(filename);
  assert(!ifs.fail());

  int ix=0,iy=0;
  int xmaxtmp,ymaxtmp;
  for(string line;getline(ifs,line);){
    if(line.empty()){
      ix++;
      ymaxtmp = iy;
      iy=0;
      continue;
    }
    if(!isdigit(line[0])) continue;

    stringstream strs(line);
    Vec3 v;
    strs >> v;
    value[ix][iy++] = v;
    if(v.z < dosmin) dosmin = v.z;
  }
  xmaxtmp = ix;

  xmax = ymax = 0;
  xmin = XMAX;
  ymin = YMAX;
  for(int x=0;x<xmaxtmp;x++)
    for(int y=0;y<ymaxtmp;y++)
      if(value[x][y].z != dosmin){
	if(x > xmax) xmax = x;
	if(x < xmin) xmin = x;
	if(y > ymax) ymax = y;
	if(y < ymin) ymin = y;
      }
}

void DoS::CalcFreeEnergy
(const double temperature,
 const double pressure){
  gmin.z = 1e32;

  for(int x=xmin;x<=xmax;x++){
    for(int y=ymin;y<=ymax;y++){
      const Vec3 v = value[x][y];
      const double gtmp = (v.z > dosmin) ? -v.z + (v.y + v.x*pressure)*nmol/temperature : 0;
      g[x][y] = Vec3(v.x,v.y,gtmp);
      if(gtmp < gmin.z) gmin = g[x][y];
    }
  }
  for(int x=xmin;x<=xmax;x++)
    for(int y=ymin;y<ymax;y++)
      gs[x][y] = g[x][y];
}

void DoS::Smoothing(){
  double g_tmp[xmax][ymax];
  double gmin_tmp = 0.0;
  for(int x=xmin;x<=xmax;x++)
    for(int y=ymin;y<=ymax;y++) g_tmp[x][y] = gs[x][y].z;

  for(int x=xmin+1;x<=xmax-1;x++){
    for(int y=xmin+1;y<=ymax-1;y++){
      double sum = 4.0 * g_tmp[x][y];
      //if(sum == 0.0) continue;
      double count = 4.0;
      const double g0 = g_tmp[x-1][y];
      const double g1 = g_tmp[x+1][y];
      const double g2 = g_tmp[x][y-1];
      const double g3 = g_tmp[x][y+1];

      if(g0 != 0){ sum += g0; count += 1.0;}
      if(g1 != 0){ sum += g1; count += 1.0;}
      if(g2 != 0){ sum += g2; count += 1.0;}
      if(g3 != 0){ sum += g3; count += 1.0;}

      const double gtmp = sum / count;
      gs[x][y].z = gtmp;
      if(gtmp < gmin_tmp) gmin_tmp = gtmp;
    }
  }
  const double scaler = gmin.z / gmin_tmp;
  for(int x=xmin;x<=xmax;x++)
    for(int y=ymin;y<=ymax;y++) gs[x][y].z *= scaler;
}

void DoS::SearchSaddlePoint(){
  extrema.erase(extrema.begin(),extrema.end());

  const Vec3 zi(0.0,0.0,1.0);
  for(int x=0;x<=xmax-1;x++){
    for(int y=0;y<=ymax-1;y++){
      double sum = 0.0, count = 0.0;
#if 1
      const Vec3 v[4] = {gs[x][y] - gs[x][y+1],
			 gs[x][y] - gs[x-1][y],
			 gs[x][y] - gs[x][y-1],
			 gs[x][y] - gs[x+1][y]};
      Vec3 n[4] = {cross(v[0],v[1]),
		   cross(v[1],v[2]),
		   cross(v[2],v[3]),
		   cross(v[3],v[0])};
      for(int i=0;i<4;i++) n[i] = n[i] / norm(n[i]);
      for(int i=0;i<4;i++) sum += dot(n[i],zi);
      count = 4.0;
#else
      const Vec3 v = gs[x][y];
      for(int i=-1;i<=1;i++){
	for(int j=-1;j<=1;j++){
	  const Vec3 v0 = v - gs[x+i][y+j];
	  if(v0.z == 0) continue;
	  for(int k=-1;k<=1;k++){
	    for(int l=-1;l<=1;l++){
	      const Vec3 v1 = v - gs[x+k][y+l];
	      if(v0 == v1 || v1.z == 0) continue;
	      const Vec3 n = cross(v0,v1);
	      sum += abs(dot(n/norm(n),zi));
	      count += 1.0;
	}}}}
#endif
      slope[x][y] = (g[x][y].z != 0.0) ? sum/count : 0;
    }
  }

  const int range = 5;
  for(int x=range;x<=xmax-range;x++){
    for(int y=range;y<=ymax-range;y++){
      const double s = slope[x][y];
      if(s == 0) continue;
      bool isLarger = true;
      for(int i=-range;i<=range;i++)
	for(int j=-range;j<=range;j++)
	  isLarger = isLarger && (s >= slope[x+i][y+j]);
      if(isLarger) peak.push_back(Peak(x,y,s));
    }
  }
  sort(peak.begin(),peak.end(),std::greater<Peak>());

  for(auto &it : peak)
    extrema.push_back(g[it.x][it.y]);
}

void DoS::SearchLocalMinimum(){
  minima.erase(minima.begin(),minima.end());
  const int range = 10;
  const Vec3 k(0.0,0.0,1.0);
  for(int x=range;x<=xmax-range;x++){
    for(int y=range;y<=ymax-range;y++){
      if(g[x][y].z == 0) continue;
      bool isSmaller = true;
      for(int i=-range;i<=range;i++)
	for(int j=-range;j<=range;j++)
	  isSmaller = isSmaller && (gs[x+i][y+j].z >= gs[x][y].z);
      if(isSmaller) minima.push_back(g[x][y]);
    }
  }
  sort(minima.begin(),minima.end());
}
/*
void DoS::Contouring(){
  const double dz = 0.2;
  int i=0;
  for(double z = gmin;z < gmin+goffset;z += dz,i++){
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
//*/

void DoS::Clustering(){
  cout << "--- start clustering ---" << endl;
  cluster.erase(cluster.begin(),cluster.end());
  border.erase(border.begin(),border.end());
  minima.erase(minima.begin(),minima.end());
  saddle.z = 1e32;

  Vec3 prev = gs[0][0];
  for(int x=1;x<=xmax-1;x++){
    for(int y=1;y<=ymax-1;y++){
      //if(gs[x][y].z < gmin.z || (gmin.z+goffset) <= gs[x][y].z) continue;
      if(gs[x][y].z == 0.0) continue;
      Vec3 v = gs[x][y];
      int ix = x, iy = y;
      while(true){
	if(gs[ix-1][iy].z < gs[ix][iy].z){
	  ix--; continue;
	}else if(gs[ix+1][iy].z < gs[ix][iy].z){
	  ix++; continue;
	}else if(gs[ix][iy-1].z < gs[ix][iy].z){
	  iy--; continue;
	}else if(gs[ix][iy+1].z < gs[ix][iy].z){
	  iy++; continue;
	}else{
	    break;
	}
      }
      auto it = find(minima.begin(),minima.end(),gs[ix][iy]);
      if(gs[ix][iy] != prev && g[x][y].z != 0){
	if(v.z < saddle.z) saddle = g[x][y];
	//if(gmin.z <= prevg && prevg < (gmin.z+goffset))
	border.push_back(g[x][y]);
      }
      //if(it == local.end() && gs[ix][iy].z < 0.5*gmin.z)
      if(it == minima.end())
	minima.push_back(g[x][y]);
      prev = gs[x][y];
    }
  }
  cout << "--- writing local_minima ---" << endl;
  ofstream lm("local_minima.dat");
  for(auto &it : minima)
    lm << it << endl;

  ofstream bd("border.dat");
  for(auto &it : border) bd << it << endl;

  cout << "--- writing cluster ---" << endl;
  ofstream cl("cluster.dat");
  if(cl.fail()) cerr << "error: opening cluster.dat failed ---" << endl;
  for(auto &it : minima){
    for(int x=1;x<=xmax-1;x++){
      for(int y=1;y<=ymax-1;y++){
	//if(gs[x][y].z < gmin.z || (gmin.z+goffset) <= gs[x][y].z) continue;
	Vec3 v = gs[x][y];
	if(v.z == 0) continue;
	int ix = x, iy = y;
	while(true){
	  if(gs[ix-1][iy].z < gs[ix][iy].z){
	    ix--; continue;
	  }else if(gs[ix+1][iy].z < gs[ix][iy].z){
	    ix++; continue;
	  }else if(gs[ix][iy-1].z < gs[ix][iy].z){
	    iy--; continue;
	  }else if(gs[ix][iy+1].z < gs[ix][iy].z){
	    iy++; continue;
	  }else{
	    if(it == g[ix][iy]) cl << g[x][y] << endl;
	    break;
	  }
	}
      }
    }
    cl << endl << endl;
  }
}

const bool isOnPlane(const Vec3 p,const Vec3 n,const Vec3 p0,const double offset){
  if(n.z != 0){
    const Vec3 c = {-n.x/n.z, -n.y/n.z, p0.z};
    const double v = c.x*(p.x-p0.x) + c.y*(p.y-p0.y) + c.z;
    return (v-0.5*offset < p.z && p.z < v+0.5*offset);
  }else{
    const Vec3 c = {-n.x/n.y, p0.y, 0.0};
    const double v = c.x * (p.x-p0.x) + c.y;
    return (v-0.5*offset < p.y && p.y < v+0.5*offset);
  }
}

void DoS::CrossSection(){
  cross_section.erase(cross_section.begin(),cross_section.end());

  const double offset = 0.1;
#if 0
  const Vec3 v0 = extrema[0] - extrema[1];
  const Vec3 v1 = extrema[0] - extrema[2];
#else
  const Vec3 v0 = minima[0] - minima[1];
  const Vec3 v1 = Vec3(0,0,1);
#endif
  Vec3 n  = cross(v0/norm(v0),v1/norm(v1));

  for(int x=xmin;x<=xmax;x++){
    for(int y=ymin;y<=ymax;y++){
      if(g[x][y].z >= 0) continue;
      const Vec3   gtmp  = g[x][y];
      if(isOnPlane(gtmp,n,minima[0],offset))
	if(min(minima[0].x,minima[1].x) <= gtmp.x && gtmp.x <= max(minima[0].x,minima[1].x))
	  if(min(minima[0].y,minima[1].y) <= gtmp.y && gtmp.y <= max(minima[0].y,minima[1].y))
	    cross_section.push_back(gtmp);
    }
  }
  sort(cross_section.begin(),cross_section.end(),
       [](const Vec3 a, const Vec3 b){return a.x > b.x;});
}

void DoS::ShortestPath(){
  shortest_path.erase(shortest_path.begin(),shortest_path.end());
  extrema.erase(extrema.begin(),extrema.end());
  if(minima.size() < 2) return;

  shortest_path.push_back(minima[0]);
  const int npoint = 40;
  const auto gm = minima[0];
  const auto lm = minima[1];
  Vec3 v = (lm - gm) / (double)npoint;
  v.z = 0;
  const double offset = norm(v);
  for(int i=0; i<npoint; i++){
    const Vec3 p = gm + v * i;
    vector<Vec3> plane;
    for(int x=xmin;x<=xmax;x++)
      for(int y=ymin;y<=ymax;y++)
	if(isOnPlane(g[x][y],v,p,offset) && g[x][y].z!=0) plane.push_back(g[x][y]);
    sort(plane.begin(),plane.end());
    if(!plane.empty()) shortest_path.push_back(plane[0]);
  }
  shortest_path.push_back(minima[1]);

  Vec3 max = gmin;;
  for(auto &it : shortest_path) max = (max<it) ? it : max;
  saddle = max;
  extrema.push_back(max);
}

void DoS::Output(ostream &os){
  for(int x=xmin;x<=xmax;x++){
    for(int y=ymin;y<=ymax;y++){
      os << value[x][y] << ' ' << g[x][y].z << ' ' << gs[x][y].z << ' ' << slope[x][y] << ' ' << g[x][y].z - gmin.z << endl;
    }
    os << endl;
  }
  os << endl << endl;
  os << "#local minimum" << endl;
  for(auto &it : minima) os << it << ' ' << it.z - gmin.z << endl;
  os << endl << endl;

  os << "#saddle point" << endl;
  for(auto &it : extrema)
    if(find(minima.begin(),minima.end(),it) == minima.end())
      os << it << ' ' << it.z - gmin.z << endl;
  os << endl << endl;

  os << "#shortest path" << endl;
  for(auto &it : shortest_path)
      os << it << ' ' << it.z - gmin.z << endl;
  os << endl << endl;
}

//#define DOS_TEST
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
  //dos->SearchSaddlePoint();
  ofstream ofs(output);
  dos->Output(ofs);
  for(int i=0;i<10;i++) dos->Smoothing();
  dos->Clustering();

}

#endif

#endif
