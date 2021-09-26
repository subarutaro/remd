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
    Sort,
    Neigh,
    M2A,
    Switching,
#ifdef INSERT_TIMER_FORCE
    Pre,
    Force,
    Post,
#endif
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
    "  Sort  ",
    "  Neigh ",
    "  M2A   ",
    "  Nb    ",
#ifdef INSERT_TIMER_FORCE
    "    Pre  ",
    "    Force ",
    "    Post  ",
#endif
    "  Wall  ",
    "  A2M   ",
    "Total  "
  };

  std::chrono::time_point<std::chrono::system_clock> s[NumProf],e[NumProf];
  std::chrono::duration<double> elapsed[NumProf];
  unsigned long long count[NumProf];

  Profiler(){
#ifdef INSERT_TIMER_FORCE
    fprintf(stderr,"WARNING: inserting profiler in force can cause significant performance drop. You should remove INSERT_TIMER_FORCE macro.\n");
#endif
    clear();
  }
  ~Profiler(){
    //print();
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

  Profiler operator+(const Profiler& rhs) const {
    Profiler ret;
    for(int i=0;i<NumProf;i++){
      ret.elapsed[i] = this->elapsed[i] + rhs.elapsed[i];
    }
    return ret;
  }
};


#endif
