#ifndef PROFILER_H
#define PROFILER_H

#include <string>
#include <chrono>

#ifdef _OPENMP
#include <omp.h>
constexpr int NTHREAD_MAX = 256;
#else
constexpr int NTHREAD_MAX = 1;
#endif

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
    Head,
    Tail,
    Post,
#endif
    Wall,
    A2M,
#ifdef ENABLE_AOS_TO_SOA_CONVERSION
    AoStoSoA,
    SoAtoAoS,
#endif
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
    "  Sort     ",
    "  Neigh    ",
    "  M2A      ",
    "  Nb       ",
#ifdef INSERT_TIMER_FORCE
    "    Pre      ",
    "    Force    ",
    "      Head   ",
    "      Tail   ",
    "    Post   ",
#endif
    "  Wall     ",
    "  A2M      ",
#ifdef ENABLE_AOS_TO_SOA_CONVERSION
    "  AoStoSoA ",
    "  SoAtoAoS ",
#endif
    "Total  "
  };

  std::chrono::time_point<std::chrono::system_clock> s[NumProf][NTHREAD_MAX],e[NumProf][NTHREAD_MAX];
  std::chrono::duration<double> elapsed[NumProf][NTHREAD_MAX];
  unsigned long long count[NumProf][NTHREAD_MAX];
  int nthreads;

  Profiler(){
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#else
    nthreads = 1;
#endif

#ifdef INSERT_TIMER_FORCE
    fprintf(stderr,"WARNING: inserting profiler in force can cause significant performance drop. You should remove INSERT_TIMER_FORCE macro.\n");
#endif
    clear();
  }
  ~Profiler(){
    //print();
  }
  void clear(){
    for(int thread=0;thread<NTHREAD_MAX;thread++){
      for(int i=0;i<NumProf;i++){
	elapsed[i][thread] = s[i][thread] - s[i][thread];
	count[i][thread] = 0;
      }
    }
  }

  void beg(const int i){
#ifdef _OPENMP
    const int thread = omp_get_thread_num();
#else
    const int thread = 0;
#endif
    s[i][thread] = std::chrono::system_clock::now();
  }
  void end(const int i){
#ifdef _OPENMP
    const int thread = omp_get_thread_num();
#else
    const int thread = 0;
#endif
    e[i][thread] = std::chrono::system_clock::now();
    elapsed[i][thread] += e[i][thread] - s[i][thread];
    count[i][thread]++;
  }

  void print(FILE* fp = stdout){
    fprintf(fp,"--- TimeProfile in sec ---\n");
    for(int i=0;i<NumProf;i++){
      fprintf(fp,"%s",name[i]);
      for(int thread=0;thread<nthreads;thread++) fprintf(fp," %e",elapsed[i][thread].count());
      fprintf(fp,"\n");
    }
  }

  Profiler operator+(const Profiler& rhs) const {
    Profiler ret;
    for(int t=0;t<NTHREAD_MAX;t++){
      for(int i=0;i<NumProf;i++){
	ret.elapsed[i][t] = this->elapsed[i][t] + rhs.elapsed[i][t];
      }
    }
    return ret;
  }
};


#endif
