#ifndef H_CDV_FORMAT
#define H_CDV_FORMAT

#include <cmath>
#include <vector>

#define MAX_ATOMTYPE 4
#define MAX_BONDTYPE 3

using namespace std;

class CDV{
 public:
  int i; // index of atom
  int m; // index of molecule
  int t;
  double r[3];
  void Print(ostream &os){
    os << i << " " << t << " " << r[0] << " " << r[1] << " " << r[2] << endl;
  }
};

class Color{
 public:
  double r;
  double g;
  double b;

  Color(double _r,double _g,double _b){set(_r,_g,_b);};
  Color(){r=1.0;g=1.0;b=0.0;}; // default color is yellow
  void set(double _r,double _g,double _b){r=_r;g=_g;b=_b;};
};

class Box{
 public:
  double s;
  double e;
};

class Header{
 public:
  Box x;
  Box y;
  Box z;

  Color  c[MAX_ATOMTYPE];
  double r[MAX_ATOMTYPE];

  Color  b[MAX_BONDTYPE];
  double w[MAX_BONDTYPE];

  Header(){
    x.s = 0.0;
    x.e = 100.0;
    y.s = 0.0;
    y.e = 100.0;
    z.s = 0.0;
    z.e = 100.0;

    for(int i=0;i<MAX_ATOMTYPE;i++){
      c[i].set(1.0,0.0,0.0);
      r[i] = 0.5;
    }

    for(int i=0;i<MAX_BONDTYPE;i++){
      b[i].set(1.0,1.0,0.0);
      w[i] = 0.1;
    }
  };
  void Print(std::ostream &os){
    os << "'";
    os << " box_sx=" << x.s << " box_ex=" << x.e;
    os << " box_sy=" << y.s << " box_ey=" << y.e;
    os << " box_sz=" << z.s << " box_ez=" << z.e;
    for(int i=0;i<MAX_ATOMTYPE;i++){
      os << " r" << i << "=" << r[i];
      os << " c" << i << "=(" << c[i].r << "," << c[i].g << "," << c[i].b << ")";
    }
    for(int i=0;i<MAX_BONDTYPE;i++){
      os << " bond" << i << "_c=(" << b[i].r << "," << b[i].g << "," << b[i].b << ")";
      os << " bond" << i << "_wt=" << w[i];
    }
    os << endl;
  };
};


#endif
