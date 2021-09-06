#include <iosfwd>
#include <iostream>
#include <sstream>
#include <string>
#include <cassert>
#include <cmath>

#include <vector>
#include <map>

#include "io.h"
#include "bond.h"
#include "open_cylinder.h"
#include "cdv_format.h"

using namespace std;

int main(){
  Header header;
  vector<CDV>  cdv;
  //vector<Bond> bond;
  map< pair<int,int>, int> bond;

  read_cdv(cin,header,cdv,bond);

  make_bonds(cdv,bond);

  double center[3];
  center[0] = 0.5*(header.x.e - header.x.s);
  center[1] = 0.5*(header.y.e - header.y.s);
  center[2] = 0.5*(header.z.e - header.z.s);
  open_cylinder(cdv,bond,center);

  write_cdv(cout,header,cdv,bond);

  return 0;
}
