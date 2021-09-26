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
  I vmax, vmin;
  I pmax, pmin;
  I vnum,pnum;
  I sum;

  I **histogram;
  Average **ave;

 public:
  Histogram(I vnum,I pnum);
  ~Histogram();
  void Increment(I v, I p);
  //void Output(char *filename);

  I GetSum(){return sum;};
  I GetHist(I v,I p){return histogram[v][p];};

  void KeepAverages(I v,I p,Average _ave){ave[v][p] += _ave;};
  Average GetAverages(I v,I p){
    if(histogram[v][p]>0) return ave[v][p]/(D)histogram[v][p];
    Average tmp;tmp.flush();
    return tmp;
  };

  I GetVolMin(){return vmin;};
  I GetVolMax(){return vmax;};
  I GetPotMin(){return pmin;};
  I GetPotMax(){return pmax;};
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
  I sum_accept;
  I sum_reject;
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
  I a;
  I b;
}Pair;

class ExchangeList{
 private:
  Pair **list;
  AcceptRate **rate;
  I    dim_temp;
  I    dim_press;
  I    length[4];
 public:
  ExchangeList(I in_dim_temp,I in_dim_press);
  ~ExchangeList();
  I GetLength(I type){
    if(type>3){
      fprintf(stderr,"error: index of GetLength(=%d) must be less than 4\n",type);
      exit(0);
    }
    return length[type];
  };
  Pair GetPair(I type,I ind){return list[type][ind];};
  void Accept(I type,I ind){rate[type][ind].Accept();};
  void Reject(I type,I ind){rate[type][ind].Reject();};
  D GetRate(I type,I ind){return rate[type][ind].GetRate();};
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

  I nmol;
  I nreplica;
  I nreplica_global; // for MPI

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

  I ensemble;

  I nproc;
  I rank;
  I ngpu;

  I mode;

  //char header[256];
  char input_dir[256];
  char output_dir[256];
  int  restart;

  I nenergy;
  I nbkup;

  I bkup_ninterval;
  I output_ninterval;

  bool *isExchanged;
  Histogram   **histogram;

  I dim_temp,dim_press;
  I ec_type;
  ExchangeList *pairlist;
  TunnelCount **tc;

  void SetConditionGeometrical(I dim_temp,I dim_press);
  void SetConditionArithmetic(I dim_temp,I dim_press);
  void SetConditionFromFile(const char*,I,I);
  void SetConditionFromHeatCapacity(std::string,I,I,const bool);
  void SetConditionFromHeatCapacity2(std::string,I,I);

  void ExchangeIndex(Pair pair){
    I tmp=index[pair.a];
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
  I *index;
  D *temperature;
  D *pressure;

  REMDInfo(const Parameter param);
  ~REMDInfo();
  void IncrementHistogram(D vol,D pot,I ind){
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
    int v = (I)(log(vol - volume_min) / delta_volume);
    int p = (I)(log(pot - energy_min) / delta_energy);
#else
    const int v = (I)((vol - volume_min) / delta_volume);
    const int p = (I)((pot - energy_min) / delta_energy);
#endif
    histogram[ind]->Increment(v,p);
  };
  void KeepAverages(D vol,D pot,I ind, Average a){
#ifdef LOG_SAMPLING
    const int v = (I)(log(vol - volume_min) / delta_volume);
    const int p = (I)(log(pot - energy_min) / delta_energy);
#else
    const int v = (I)((vol - volume_min) / delta_volume);
    const int p = (I)((pot - energy_min) / delta_energy);
#endif
    histogram[ind]->KeepAverages(v,p,a);
  }

#ifdef LOG_SAMPLING
  D GetAbsoluteVolume(I idx){return exp(delta_volume*((D)idx + 0.5)) + volume_min;};
  D GetRelativeVolume(I idx){return exp(delta_volume*((D)idx + 0.5));};
  D GetAbsoluteEnergy(I idx){return exp(delta_energy*((D)idx + 0.5)) + energy_min;};
  D GetRelativeEnergy(I idx){return exp(delta_energy*((D)idx + 0.5));};
#else
  D GetAbsoluteVolume(I idx){return delta_volume*((D)idx + 0.5) + volume_min;};
  D GetRelativeVolume(I idx){return delta_volume*((D)idx + 0.5);};
  D GetAbsoluteEnergy(I idx){return delta_energy*((D)idx + 0.5) + energy_min;};
  D GetRelativeEnergy(I idx){return delta_energy*((D)idx + 0.5);};
#endif

  void ExchangeAccepted(I type,I ind){
    Pair pair = GetPair(type,ind);
    isExchanged[index[pair.a]] = true;
    isExchanged[index[pair.b]] = true;
    ExchangeTemperature(pair);
    ExchangePressure(pair);
    ExchangeIndex(pair);
    pairlist->Accept(type,ind);
  };
  void ExchangeRejected(I type,I ind){
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

  void Permutate(I ind[],I ilist[],D t[],I n,I N){
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
  void SetTemperature(I ind,D T){temperature[ind]=T;};
  void SetPressure(I ind,D P){pressure[ind]=P;};
  void SetAllIsExchanged(bool b){for(int i=0;i<nreplica;i++) isExchanged[i]=b;};
  void SetAllIsExchangedGlobal(bool b){for(int i=0;i<nreplica_global;i++) isExchanged[i]=b;};
  
  //getter
  I GetNmol(){return nmol;};
  unsigned long GetStepMax(){return step_max;};
  unsigned long GetInterval(){return interval;};
  I GetNumReplica(){return nreplica;};
  I GetNumReplicaGlobal(){return nreplica_global;};
  I GetReplicaOffset(){ return rank * nreplica; }
  D GetTemperatureMax(){return temperature_max;};
  D GetTemperatureMin(){return temperature_min;};
  D GetPressureMax(){return pressure_max;};
  D GetPressureMin(){return pressure_min;};
  D GetDeltaEnergy(){return delta_energy;};
  D GetEnergyMax(){return energy_max;};
  D GetEnergyMin(){return energy_min;};
  D GetVolumeMax(){return volume_max;};
  D GetVolumeMin(){return volume_min;};
  I GetNumProc(){return nproc;};
  I GetRank(){return rank;};
  I GetNumGPU(){return ngpu;};
  I GetMode(){return mode;};
  I GetEnsemble(){return ensemble;};
  char* GetInputDir(){return input_dir;};
  char* GetOutputDir(){return output_dir;};
  I GetNumBackup(){return nbkup;};
  I GetNumPointEnergy(){return nenergy;};

  I GetListType(){return ec_type;};
  I GetListLength(I type){return pairlist->GetLength(type);};//get length of type th exchange list
  Pair GetPair(I type,I ind){return pairlist->GetPair(type,ind);};//get exchange ind th pair of type th pair list
  I GetIndex(I ind){return index[ind];};            //get replica number of ind th temperature
  D GetTemperature(I ind){return temperature[ind];};//get temperature of ind th replica
  D GetPressure(I ind){return pressure[ind];};      //get pressure    of ind th replica
  bool GetIsExchanged(I ind){return isExchanged[ind];};//get isExchanged of ind th replica
  I GetSumHist(I ind){return histogram[ind]->GetSum();};
  I GetHistogram(I ind,I v,I p){return histogram[ind]->GetHist(v,p);};
  D GetVirialAve(I v,I p);
  //D GetPressAve(I v,I p);
  //D GetTempAve(I v,I p);
  Average GetAverages(I,I);

  void ReadTemperatureFromFile(char *filename);

  void BroadcastConditionAndIndex();

  void ShowAll();
};

#endif





