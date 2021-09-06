#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <string>
#include <utility>

using namespace std;

int main(int argc,char** argv){
  if(argc != 3){
    cerr << "usage: " << argv[0] << " [input file] [output file]" << endl;
    exit(EXIT_FAILURE);
  }
  const string input  = argv[1];
  const string output = argv[2];

  vector<pair<double,double> > coex;

  ifstream ifs(input);
  pair<double,double> eq;
  double lastc = 0.0;
  double lastp = 0.0;
  for(string line;std::getline(ifs,line);){
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

  ofstream ofs(output);
  for(auto &it : coex) ofs << it.first << ' ' << it.second << endl;
}
