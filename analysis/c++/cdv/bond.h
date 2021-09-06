#ifndef H_BOND
#define H_BOND

#include <iostream>
#include <vector>
#include <map>
#include "cdv_format.h"

using namespace std;

void make_bonds(vector<CDV> cdv, map< pair<int,int>, int> &bond);
#endif
