#ifndef H_OPEN_CYLINDER
#define H_OPEN_CYLINDER

#include <iostream>
#include <fstream>
#include <map>
#include <vector>

#include "cdv_format.h"

using namespace std;

void calc_outer_inner(vector<CDV> cdv, double center[3], double &outer, double &inner);
void define_outer_inner(vector<CDV> &cdv, double center[3],const double outer, const double inner);
void open_cylinder(vector<CDV> &cdv, map< pair<int,int>, int> &bond,double center[3],const double outer);
void delete_periodic_bond(vector<CDV> cdv,map< pair<int,int>, int> &bond,double center[3],const double outer);
void separate_inter_outer(vector<CDV> cdv,double center[3]);
#endif
