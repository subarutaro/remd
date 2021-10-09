#ifndef REMDINFO_H
#define REMDINFO_H

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <iostream>
#include <fstream>
#include <iomanip>

#include "parameter.h"
#include "unit.h"

#include "mdunit.h"
#include "typedef.h"

#include "integrator.h"

//#define HISTOGRAM1D

#define NUM_POT 500
#define NUM_VOL 500

//#define LOG_SAMPLING

class Average{
 public:
  enum{
    LJ = 0,CLMB,WALL,
    VOL,
    TRA,ROT,
    T,P,
    Px,Py,Pz,
    TOT,DRF,
    TSTs,TSTv,
    BSTx,BSTy,BSTz,
    // don't touch the below
    DUMMY,
    NUM_ELEM,
  };
  static const char *name(const int i){
    static const char *strs[NUM_ELEM] = {
      "lj          ","coulomb     ","wall        ",
      "pot         ",
      "vol         ",
      "temperature ",
      "pressure x  ","pressure y  ","pressure z  ",
      "tra         ","rot         ",
      "total       ","drift       ",
      "tst_s       ","tst_v       ",
      "bst_x       ","bst_y       ","bst_z       ",
    };
    return strs[i];
  }

  D sum[NUM_ELEM];

  void operator+=(Average a){
    for(int i=0;i<NUM_ELEM;i++) this->sum[i] += a.sum[i];
  }
  void operator+=(D s){
    for(int i=0;i<NUM_ELEM;i++) this->sum[i] += s;
  }
  void operator-=(Average a){
    for(int i=0;i<NUM_ELEM;i++) this->sum[i] -= a.sum[i];
  }
  void operator-=(D s){
    for(int i=0;i<NUM_ELEM;i++) this->sum[i] -= s;
  }
  Average operator*(D s){
    Average tmp;
    for(int i=0;i<NUM_ELEM;i++) tmp.sum[i] = this->sum[i]*s;
    return tmp;
  }
  Average operator/(D s){
    Average tmp;
    for(int i=0;i<NUM_ELEM;i++) tmp.sum[i] = this->sum[i]/s;
    return tmp;
  }

  void flush(){
    for(int i=0;i<NUM_ELEM;i++) sum[i] = 0.0;
  };
  void show(FILE *fp = stderr,
	    const char *fmt = " %s : %e \n"){
    for(int i=0;i<NUM_ELEM;i++){
      fprintf(fp,fmt,name(i),sum[i]);
    }
  };
};

class Histogram{
 private:
  int vmax, vmin;
  int pmax, pmin;
  int vnum,pnum;
  int sum;

  int **histogram;
  Average **ave;

 public:
  Histogram(int vnum,int pnum);
  ~Histogram();
  void Increment(int v, int p);
  //void Output(char *filename);

  int GetSum(){return sum;};
  int GetHist(int v,int p){return histogram[v][p];};

  void KeepAverages(int v,int p,Average _ave){ave[v][p] += _ave;};
  Average GetAverages(int v,int p){
    if(histogram[v][p]>0) return ave[v][p]/(D)histogram[v][p];
    Average tmp;tmp.flush();
    return tmp;
  };

  int GetVolMin(){return vmin;};
  int GetVolMax(){return vmax;};
  int GetPotMin(){return pmin;};
  int GetPotMax(){return pmax;};
  void SetVolMin(int _vmin){vmin = _vmin;};
  void SetVolMax(int _vmax){vmax = _vmax;};
  void SetPotMin(int _pmin){pmin = _pmin;};
  void SetPotMax(int _pmax){pmax = _pmax;};
  void SetHist(int v,int p,int h){histogram[v][p] = h;};
  void SetAve(int v,int p,Average a){ave[v][p] = a * (double)histogram[v][p];};
  void SetSum(int s){sum = s;};
  
  void flush(){
    for(int v=vmin;v<=vmax;v++){
      for(int p=pmin;p<=pmax;p++){
	histogram[v][p] = 0;
	sum = 0;
	ave[v][p].flush();
    }}
  };
};

class AcceptRate{
 private:
  int sum_accept;
  int sum_reject;
 public:
  AcceptRate(){
    sum_accept=0;sum_reject=0;
  };
  ~AcceptRate(){};
  void Accept(){sum_accept++;};
  void Reject(){sum_reject++;};
  D GetRate(){return (D)sum_accept/(D)(sum_accept+sum_reject);};
  void flush(){sum_accept = sum_reject = 0;};
};

typedef struct{
  int a;
  int b;
}Pair;

class ExchangeList{
 private:
  Pair **list;
  AcceptRate **rate;
  int    dim_temp;
  int    dim_press;
  int    length[4];
 public:
  ExchangeList(int in_dim_temp,int in_dim_press);
  ~ExchangeList();
  int GetLength(int type){
    if(type>3){
      fprintf(stderr,"error: index of GetLength(=%d) must be less than 4\n",type);
      exit(0);
    }
    return length[type];
  };
  Pair GetPair(int type,int ind){return list[type][ind];};
  void Accept(int type,int ind){rate[type][ind].Accept();};
  void Reject(int type,int ind){rate[type][ind].Reject();};
  D GetRate(int type,int ind){return rate[type][ind].GetRate();};
  void FlushRate(){
    for(int type=0;type<4;++type){
      for(int ind=0;ind<length[type];ind++){
	rate[type][ind].flush();
      }
    }
  };
};

class TunnelCount{
 private:
  int count;
  bool max,min;

  bool isTunnelHappen(){
    if(max && min) return true;
    else return false;
  }
 public:
  TunnelCount(){flush();}
  void flush(){max = min = false;count = 0;};
  void ReachMax(){
    max = true;
    if(isTunnelHappen()){
      count++;
      min = false;
    }
  }
  void ReachMin(){
    min = true;
    if(isTunnelHappen()){
      count++;
      max = false;
    }
  }
  int GetCount(){return count;};
};

class REMDInfo{
 private:
  unsigned long step_max;
  unsigned long interval;

  int nmol;
  int nreplica;
  int nreplica_global; // for MPI

  D temperature_max;
  D temperature_min;

  D pressure_max;
  D pressure_min;

  D delta_energy;
  D energy_max;
  D energy_min;

  D delta_volume;
  D volume_max;
  D volume_min;

  int ensemble;

  int nproc;
  int rank;
  int ngpu;

  int mode;

  //char header[256];
  char input_dir[256];
  char output_dir[256];
  int  restart;

  int nenergy;
  int nbkup;

  int bkup_ninterval;
  int output_ninterval;

  bool *isExchanged;
  Histogram   **histogram;

  int dim_temp,dim_press;
  int ec_type;
  ExchangeList *pairlist;
  TunnelCount **tc;

  void SetConditionGeometrical(int dim_temp,int dim_press);
  void SetConditionArithmetic(int dim_temp,int dim_press);
  void SetConditionFromFile(const char*,int,int);
  void SetConditionFromHeatCapacity(std::string,int,int,const bool);
  void SetConditionFromHeatCapacity2(std::string,int,int);

  void ExchangeIndex(Pair pair){
    int tmp=index[pair.a];
    index[pair.a]=index[pair.b];
    index[pair.b]=tmp;
  };
  void ExchangeTemperature(Pair pair){
    D tmp=temperature[index[pair.a]];
    temperature[index[pair.a]]=temperature[index[pair.b]];
    temperature[index[pair.b]]=tmp;
  };
  void ExchangePressure(Pair pair){
    D tmp=pressure[index[pair.a]];
    pressure[index[pair.a]]=pressure[index[pair.b]];
    pressure[index[pair.b]]=tmp;
  };

 public:
  int *index;
  D *temperature;
  D *pressure;

  REMDInfo(const Parameter param);
  ~REMDInfo();
  void IncrementHistogram(D vol,D pot,int ind){
    if(vol <  volume_min ||
       vol >= volume_max ||
       pot <  energy_min ||
       pot >= energy_max ){
      fprintf(stderr,"error: potential energy or volume is not in the histogram range\n");
      fprintf(stderr,"       energy must be %lf <= %lf < %lf\n",energy_min,pot,energy_max);
      fprintf(stderr,"       volume must be %lf <= %lf < %lf\n",volume_min,vol,volume_max);
      fflush(stderr);
      assert(vol >= volume_min && vol < volume_max);
      assert(pot >= energy_min && pot < energy_max);
      exit(EXIT_FAILURE);
    }

#ifdef LOG_SAMPLING
    int v = (int)(log(vol - volume_min) / delta_volume);
    int p = (int)(log(pot - energy_min) / delta_energy);
#else
    const int v = (int)((vol - volume_min) / delta_volume);
    const int p = (int)((pot - energy_min) / delta_energy);
#endif
    histogram[ind]->Increment(v,p);
  };
  void KeepAverages(D vol,D pot,int ind, Average a){
#ifdef LOG_SAMPLING
    const int v = (int)(log(vol - volume_min) / delta_volume);
    const int p = (int)(log(pot - energy_min) / delta_energy);
#else
    const int v = (int)((vol - volume_min) / delta_volume);
    const int p = (int)((pot - energy_min) / delta_energy);
#endif
    histogram[ind]->KeepAverages(v,p,a);
  }

#ifdef LOG_SAMPLING
  D GetAbsoluteVolume(int idx){return exp(delta_volume*((D)idx + 0.5)) + volume_min;};
  D GetRelativeVolume(int idx){return exp(delta_volume*((D)idx + 0.5));};
  D GetAbsoluteEnergy(int idx){return exp(delta_energy*((D)idx + 0.5)) + energy_min;};
  D GetRelativeEnergy(int idx){return exp(delta_energy*((D)idx + 0.5));};
#else
  D GetAbsoluteVolume(int idx){return delta_volume*((D)idx + 0.5) + volume_min;};
  D GetRelativeVolume(int idx){return delta_volume*((D)idx + 0.5);};
  D GetAbsoluteEnergy(int idx){return delta_energy*((D)idx + 0.5) + energy_min;};
  D GetRelativeEnergy(int idx){return delta_energy*((D)idx + 0.5);};
#endif

  void ExchangeAccepted(int type,int ind){
    Pair pair = GetPair(type,ind);
    isExchanged[index[pair.a]] = true;
    isExchanged[index[pair.b]] = true;
    ExchangeTemperature(pair);
    ExchangePressure(pair);
    ExchangeIndex(pair);
    pairlist->Accept(type,ind);
  };
  void ExchangeRejected(int type,int ind){
    Pair pair = GetPair(type,ind);
    isExchanged[index[pair.a]] = false;
    isExchanged[index[pair.b]] = false;
    pairlist->Reject(type,ind);
  };
  void MoveToNextList(){
    // only for 2 dimension
    printf("ec_type changed from %d",ec_type);
    if(dim_temp==1 || dim_press==1){
      ec_type += 2;
    }else{
      ec_type += 1;
    }
    if(ec_type > 3) ec_type -= 4;
    printf(" to %d\n",ec_type);
  };  

  void Permutate(int ind[],int ilist[],D t[],int n,int N){
    /*
    for(int i=0;i<n;i++){
      printf("%d ind= %d,t= %e\n",i,ind[i],t[i]);
    }
    for(int i=0;i<n;i++){
      printf("%d ind= %d,t= %e,ilist= %d\n",i,ind[ilist[i]],t[ilist[i]],ilist[i]);
    }
    //*/
    for(int i=0;i<n;i++){
      temperature[N+i]=t[ilist[i]];
    }
    for(int i=0;i<n;i++){
      index[N+i] = ind[ilist[i]];
    }
  };

  void CheckTunneling(){
    tc[index[0]         ]->ReachMin();
    tc[index[nreplica-1]]->ReachMax();
  }

  void ReadHistogramFromBackUp();

  bool ConditionAdjustor();

  void OutputAcceptRate();
  void OutputHistogram();
  void OutputTunnelCount();

  //setter
  void SetTemperature(int ind,D T){temperature[ind]=T;};
  void SetPressure(int ind,D P){pressure[ind]=P;};
  void SetAllIsExchanged(bool b){for(int i=0;i<nreplica;i++) isExchanged[i]=b;};
  void SetAllIsExchangedGlobal(bool b){for(int i=0;i<nreplica_global;i++) isExchanged[i]=b;};
  
  //getter
  int GetNmol(){return nmol;};
  unsigned long GetStepMax(){return step_max;};
  unsigned long GetInterval(){return interval;};
  int GetNumReplica(){return nreplica;};
  int GetNumReplicaGlobal(){return nreplica_global;};
  int GetReplicaOffset(){ return rank * nreplica; }
  D GetTemperatureMax(){return temperature_max;};
  D GetTemperatureMin(){return temperature_min;};
  D GetPressureMax(){return pressure_max;};
  D GetPressureMin(){return pressure_min;};
  D GetDeltaEnergy(){return delta_energy;};
  D GetEnergyMax(){return energy_max;};
  D GetEnergyMin(){return energy_min;};
  D GetVolumeMax(){return volume_max;};
  D GetVolumeMin(){return volume_min;};
  int GetNumProc(){return nproc;};
  int GetRank(){return rank;};
  int GetNumGPU(){return ngpu;};
  int GetMode(){return mode;};
  int GetEnsemble(){return ensemble;};
  char* GetInputDir(){return input_dir;};
  char* GetOutputDir(){return output_dir;};
  int GetNumBackup(){return nbkup;};
  int GetNumPointEnergy(){return nenergy;};

  int GetListType(){return ec_type;};
  int GetListLength(int type){return pairlist->GetLength(type);};//get length of type th exchange list
  Pair GetPair(int type,int ind){return pairlist->GetPair(type,ind);};//get exchange ind th pair of type th pair list
  int GetIndex(int ind){return index[ind];};            //get replica number of ind th temperature
  D GetTemperature(int ind){return temperature[ind];};//get temperature of ind th replica
  D GetPressure(int ind){return pressure[ind];};      //get pressure    of ind th replica
  bool GetIsExchanged(int ind){return isExchanged[ind];};//get isExchanged of ind th replica
  int GetSumHist(int ind){return histogram[ind]->GetSum();};
  int GetHistogram(int ind,int v,int p){return histogram[ind]->GetHist(v,p);};
  D GetVirialAve(int v,int p);
  //D GetPressAve(int v,int p);
  //D GetTempAve(int v,int p);
  Average GetAverages(int,int);

  void ReadTemperatureFromFile(char *filename);

  void BroadcastConditionAndIndex();

  void ShowAll();
};

#endif





