#ifndef TIMER_H
#define TIMER_H
//copied from HERMITE by Dr. Nitadori

#include <sys/time.h>
#if 0
struct Profile{
	enum{
		LJ = 0,
		EWALDSR,
		EWALDLR,
		INTEG,
		INIT,
		IO,
		// don't touch the below
		TOTAL,
		MISC,
		NUM_ELEM,
	};

	static const char *name(const int i){
		static const char *strs[NUM_ELEM] = {
			"LJ      ",
			"EWALD/SR",
			"EWALD/LR",
			"integ   ",
			"init    ",
			"I/O     ",
			"total   ",
			"misc    ",
		};
		return strs[i];
	}

	double time[NUM_ELEM];
	double tprev;

	long num_step_prev;
	long num_bstep_prev;
	
	void flush(){
		for(int i=0; i<NUM_ELEM; i++) time[i] = 0.0;
	}
	void beg(const int elem, const bool reuse = false){
		if(reuse) time[elem] -= tprev;
		else      time[elem] -= wtime();
	}
	void end(const int elem){
		// time[elem] += wtime();
		tprev = wtime();
		time[elem] += tprev;
	}
	void show(
		  //const long nbody,
		  const long num_step_tot,
		  const long num_bstep_tot,
		  FILE *fp = stderr,
		  const char *fmt = " %s : %e sec, %e usec : %6.2f %%\n")
	{
	  //*
		const long ns  = num_step_tot  - num_step_prev;
		const long nbs = num_bstep_tot - num_bstep_prev;
		num_step_prev  = num_step_tot;
		num_bstep_prev = num_bstep_tot;
	  //*/
		time[MISC] = time[TOTAL];
		for(int i=0; i<NUM_ELEM-2; i++){
			time[MISC] -= time[i];
		}
		for(int i=0; i<NUM_ELEM; i++){
			fprintf(fp, fmt, name(i), 
					time[i], 
					time[i] * (1.e6/nbs),
					100.0 * time[i] / time[TOTAL]);
		}
		/*
		const double nact   = double(ns) / double(nbs);
		const double wtime  = time[TOTAL];
		const double Gflops = (((1.0e-9 * ns) * nbody) * Particle::flops) / wtime;
		fprintf(stdout, "## nbody  wtime(sec)  Gflops  wtime/block(usec)  nact\n");
		fprintf(stdout, "## %ld %e %e %e %e\n",
				nbody, wtime, Gflops, wtime/nbs*1.e6, nact);
		//*/
	}
};
#endif

#define PROFILE
struct Time{
  double s,e;
  char type[256];
#ifdef PROFILE
  Time(const char* _type){ sprintf(type,"%s",_type);s = e = 0.0; }
  void beg(){ s+=wtime();}
  void end(){ e+=wtime();}
  void print(){std::cout << type << ":" << e - s << "[s]" << std::endl;}
  static double wtime(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return double(tv.tv_sec) + 1.e-6*double(tv.tv_usec);
  }
#else
  Time(const char* _type){}
  void beg(){}
  void end(){}
  void print(){}
  double wtime(){return 0.0;}
#endif
};

#endif
