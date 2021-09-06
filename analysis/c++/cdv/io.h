#ifndef H_IO
#define H_IO

#include <iosfwd>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <cassert>
#include <cmath>

#include <vector>
#include <map>

#include "cdv_format.h"

using namespace std;

void read_header();
//void read_cdv(std::istream &is, Header &header, vector<CDV> &cdv, vector<Bond> &bond);
//void write_cdv(std::ostream &os, Header header,vector<CDV> cdv, vector<Bond> bond);
void read_cdv(std::istream &is, Header &header, vector<CDV> &cdv, map< pair<int,int>, int> &bond);
void write_cdv(std::ostream &os, Header header,vector<CDV> cdv, map< pair<int,int>, int> bond);
#endif
