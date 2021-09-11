#ifndef PROFILER_H
#define PROFILER_H

#include <string>
#include <chrono>

class Profiler{
 public:
  enum {
    D1 = 0,
    D2,
    D3,
    D4,
    D5,
    D6,
    CalcForce,
    M2A,
    Switching,
    Wall,
    A2M,
    Total,
    NumProf
  };
  char name[256][NumProf] = {
    "D1    ",
    "D2    ",
    "D3    ",
    "D4    ",
    "D5    ",
    "D6    ",
    "Force ",
    "  M2A   ",
    "  Nb    ",
    "  Wall  ",
    "  A2M   ",
    "Total  "
  };

  std::chrono::time_point<std::chrono::system_clock> s[NumProf],e[NumProf];
  std::chrono::duration<double> elapsed[NumProf];
  unsigned long long count[NumProf];

  Profiler(){
    clear();
  }
  ~Profiler(){
    print();
  }
  void clear(){
    for(int i=0;i<NumProf;i++){
      elapsed[i] = s[i] - s[i];
      count[i] = 0;
    }
  }
  void beg(const int i){
    s[i] = std::chrono::system_clock::now();
  }
  void end(const int i){
    e[i] = std::chrono::system_clock::now();
    elapsed[i] += e[i] - s[i];
    count[i]++;
  }

  void print(FILE* fp = stdout){
    fprintf(fp,"--- TimeProfile in sec ---\n");
    for(int i=0;i<NumProf;i++){
      fprintf(fp,"%s%e\n",name[i],elapsed[i].count());
    }
  }
};


#endif
