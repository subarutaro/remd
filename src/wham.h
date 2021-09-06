#ifndef WHAM_H
#define WHAM_H

#include "typedef.h"
#include "mdunit.h"
//#include "remd.h"
#include "remdinfo.h"

#define NumTemp 200
#define NumPress 10

class WHAM{
 private:
  D **dens_state;
  D dos_min;
  D **fe_surface;
  D *g;
  I pot_num,vol_num;
  I pot_max = 0, pot_min = 0;
  I vol_max = 0, vol_min = 0;

  D ent[NumPress][NumTemp];
  D hc[NumPress][NumTemp];
  D fe[NumPress][NumTemp];
  D potential[NumPress][NumTemp];
  D volume[NumPress][NumTemp];
  D pressure[NumPress][NumTemp];
  D temperature[NumPress][NumTemp];
  Average average[NumPress][NumTemp];

 public:
  WHAM(){};
  WHAM(REMDInfo *remdinfo);
  ~WHAM(){
    /*
    for(int i=0;i<vol_num;i++){
      free(dens_state[i]);
    }
    free(dens_state);
    free(g);
    //*/
  };

  void CalcDensState(REMDInfo *remdinfo);
  void CalcG(REMDInfo *remdinfo);
  void CalcPhysValue(REMDInfo *remdinfo);
  void Output(REMDInfo *remdinfo);
};

#endif
