#include "chk.h"

#include <cmath>
#include <cstdlib>
#include <cassert>

#include <fstream>

void CheckPointFile::Resize(const Length _l){
  const double cx = _l.x / l.x;
  const double cy = _l.y / l.y;
  const double cz = _l.z / l.z;
  for(int i=0;i<nmol;i++){
    m[i].r[0] *= cx;
    m[i].r[1] *= cy;
    m[i].r[2] *= cz;
  }
  l = _l;
}

void CheckPointFile::Resize(const int _nmol,const int dim = 3){
  assert(0 < dim && dim <= 3);

  const double c = pow((double)_nmol/(double)nmol, 1./(double)dim);
  const int n = (int)c + 1;
  const double tail = 2.0;
  const int na = natom / nmol;
  
  Length _l = l;
  const double rate = 1.2;
  if(dim != 1){
    _l.x *= c*rate;
    _l.y *= c*rate;
  }
  if(dim != 2) _l.z *= c * rate;

  int count = nmol;
  for(int nx=0;nx<n;nx++){
    for(int ny=0;ny<n;ny++){
      for(int nz=0;nz<n;nz++){
	for(int i=0;i<nmol;i++){
	  if(((nx == 0) && (ny == 0) && (nz == 0)) || (count >= _nmol)) break;
	  Molecule _m = m[i];
	  _m.r[0] += l.x * nx;
	  _m.r[1] += l.y * ny;
	  _m.r[2] += l.z * nz;
	  if((_m.r[0] < _l.x) && (_m.r[1] < _l.y) && (_m.r[2] < _l.z)){
	    for(int j=0;j<nmol;j++){
	      const double x = m[j].r[0] - _m.r[0];
	      const double y = m[j].r[1] - _m.r[1];
	      const double z = m[j].r[2] - _m.r[2];
	      const double r = sqrt(x*x + y*y + z*z);
	      if(r < 0.1){
		std::cerr << "warning: molecules are too close to each other!(" << i << ":" << nx << "," << ny << "," << nz << ")" << std::endl;
		std::cerr << "before:" << m[j] << std::endl;
		std::cerr << "after: " << _m << std::endl;
	      }
	    }
	    _m.id = count*na;
	    m.push_back(_m);
	    count++;
	  }
	}
      }}}

  if(count != _nmol){
    std::cerr << "error: nmol is not equal to " << _nmol << "(" << count << ")" << std::endl;
    std::cerr << "       Adjust parameter!" << std::endl;
    exit(EXIT_FAILURE);
  }

  if(dim != 1){
    _l.x += tail;
    _l.y += tail;
  }
  if(dim != 2) _l.z += tail;
  l = _l;
  nmol  = count;
  natom = count * na;
}

int main(int argc,char **argv){
  if(argc != 5){
    std::cerr << "usage: " << argv[0] << " [nmol] [dim] [input] [output]" << std::endl;
    exit(EXIT_FAILURE);
  }
  const int nmol = atoi(argv[1]);
  const int dim  = atoi(argv[2]);
  std::string input  = argv[3];
  std::string output = argv[4];

  std::ifstream ifs(input);
  if(ifs.fail()){
    std::cerr << "error: failed to open " << input  << std::endl;
    exit(EXIT_FAILURE);
  }

  CheckPointFile chk;
  ifs >> chk;

  chk.Resize(nmol,dim);
  
  std::ofstream ofs(output);
  if(ofs.fail()){
    std::cerr << "error: failed to open " << output << std::endl;
    exit(EXIT_FAILURE);
  }
  ofs << chk;
}
