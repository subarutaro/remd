
#ifndef REMD_H
#define REMD_H

#include "typedef.h"
//#include <md.h>
#include "molecules.h"
//typedef Molecules MD;
#include "io_manager.h"

#include "remdinfo.h"
#include "wham.h"
//#include <time_measure.h>
#ifdef ENABLE_GPU
#include "remdgpu.h"
#endif
#include "integrator.h"

class REMD{
 private:
  Parameter param;
  int ngpus;
  unsigned long step;

  D *pot = nullptr;
  D *vol = nullptr;
  D *vir = nullptr;

  D *press = nullptr;
  D *temp  = nullptr;
  Average *ave = nullptr;

  REMDInfo  *remdinfo = nullptr;
  Molecules **md      = nullptr;
  IOManager **iomngr  = nullptr;
  WHAM      *wham     = nullptr;
#ifdef ENABLE_GPU
  REMDGPU   *mdgpu    = nullptr;
#endif
  std::ostream *os_ene = nullptr;
  std::fstream *binary = nullptr;

  //function for replica exchange
  D    CalcExchangeProb(I type,I ind);
  void ReplicaExchange();
  void DesignedWalkReplicaExchange();
  //template<uint N>
  void ReplicaPermutation();
  void IncrementHistogram();

  //output
  int  snap_interval;
  int  chk_interval;
  int  energy_interval;
  void OutputHistogram();
  void OutputIndex();
  void OutputAcceptanceRatio();
  void OutputBackUp();
  void OutputEnergy();

  void OutputBinary();

public:
  REMD();
//REMD(char *inputfile);
  ~REMD();

  void PrecalcREMD();
  //void ExecuteREMDcpu(unsigned long step_max,unsigned long interval);
  //void ExecuteREMDgpuVelScale(unsigned long step_max,unsigned long interval);
  //void ExecuteREMDgpuNVE(unsigned long step_max,unsigned long interval);
  //void ExecuteREMDgpuNVT(unsigned long step_max,unsigned long interval);
  //void ExecuteREMDgpuNPT(unsigned long step_max,unsigned long interval);
  //void ExecuteREMDgpuManySample(I step_max,I interval,I ngpu);
  void ExecuteREMD();
  void ExecuteMDs();
  int CheckVarIncludedinRange();

  void ExecuteWHAM();

  //debug function
  void OutputDebugInfo();
};

#endif






