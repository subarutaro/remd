#pragma once
#include <iosfwd>
#include <iostream>
#include <sstream>
#include <string>
#include <cassert>
#include <cmath>

struct Parameter{
  typedef std::string string;

  string input_prefix;
  string output_prefix;

  string gro_in,gro_out;   // gro file prefix (postfix must be .gro)
  string trr_in,trr_out;   // trr file prefix (postfix must be .trr)
  string cdv_in,cdv_out;          // cdv file prefix (postfix must be .cdv)
  string chk_in,chk_out;
  string tkw_in;
  string xyz_in;
  string ene_out;          // energy output prefix (postfix must be .ene)
  string mtl_in;           // input file of types of molecule

  int nmol,natom;          // number of molecules and atoms
  int nsite;               // number of sites of water

  double  dt;              // delta time
  double  rcut;            // distance of cut off in direct space

  double  kcut;            // distance of cut off in wave space
  double  alpha;           // alpha for ewald sum

  double rswitch;          // start switcing funtion from rswitch to rcut

  int     tconstant;       // whether control T or not
  double  temperature;     // temperature
  double  tstat_mass;      // mass of thermostat

  int     pconstant;       // whether control P or not
  double  pressure;        // pressure
  double  bstat_mass;      // mass of barostat

  //int nstep_precalc;       // number of steps of precalculation with velocity scaling
  int ninterval;           // number of intervals
  int interval;            // number of steps for one interval
  int energy_interval;     // number of intervals of energy output
  int snap_interval;       // number of intervals of snapshot output
  int chk_interval;        // number of intervals of check point file output
  int mom_cancel_interval; // number of intervals of canceling mementum of system

  //------------------------------
  // below is for replica exchange molecular dynamics simulation
  //------------------------------
  int ngpus;               // number of gpus. if ngpu<1, cpu version will be run.
  int gpu_offset;          // the start index of devices to use
  int nthreads;            // number of OpenMP threads
  int nprocs;              // number of mpi processes

  int nreplica;            // number of replicas (= ntemp * npress)
  int ntemp;               // number of temperature condition
  int npress;              // number of pressure condition

  double temp_max;         // maximum temperature
  double temp_min;         // minimum temperature
  double press_max;        // maximum pressure
  double press_min;        // minimum pressure

  double energy_max;       // maximum potential energy for histogram
  double energy_min;       // minimum potential energy for histogram

  double volume_max;       // maximum volume for histogram
  double volume_min;       // minimum volume for histogram

  int gen_vel;             // generate velocity randomly or not
  int restart;             // read histogram(1) file or not(0)
  int adjust_center;       // move gravity center to center of system or not
  int init_scaling;        // move gravity center to center of system or not

  int cond_mode;           // generate or read temperature and pressure condition

  int rem_type;            // replica exchange type (0: normal,1: designedwalk)

  int    confined;         // confined (1 for 1D, 2 for 2D) or not(0)
  int    nfwall;           // # of force and potential to keep in array
  double wall_length;      // length of wall (radious for 1D, height for 2D)
  double sgm_wall;         // parameter for wall force
  double eps_wall;         //
  double rho_wall;         //

  int adjustor;            // use condition adjustor or not
  int adjust_interval;     // number of sampling for adjust conditions

  Parameter(){
    gro_in  = "null";
    gro_out = "null";
    chk_in  = "null";
    chk_out = "null";
    tkw_in  = "null";
    xyz_in  = "null";
    trr_in  = "null";
    trr_out = "null";
    ene_out = "null";
    cdv_in  = "null";
    cdv_out = "null";
    mtl_in  = "moltype.txt";

    nmol    = 0;
    natom   = 0;
    nsite   = 4;

    dt      = 2.0;
    rcut    = 10.0;

    kcut    = 10.0;
    alpha   = 0.30;

    rswitch = 8.0;

    tconstant  = 0;
    tstat_mass = -1.0;
    pconstant  = 0;
    bstat_mass = -1.0;

    temperature = 300.0;
    pressure    = 1.0;

    //nstep_precalc = 0;
    ninterval           = 0;
    interval            = 0;
    energy_interval     = 0;
    snap_interval       = 0;
    mom_cancel_interval = 0;

    ngpus      = 0;
    gpu_offset = 0;
    nprocs     = 1;
    nthreads   = 1;


    nreplica = 1;
    ntemp    = 1;
    npress   = 1;

    temp_max = 300.0;
    temp_min = 300.0;
    press_max = 1.0;
    press_min = 1.0;

    energy_max = -2000.;
    energy_min = -2000.;
    volume_max = 200.;
    volume_min = 200.;

    gen_vel   = 0;
    cond_mode = 0;
    restart   = 0;
    adjust_center = 0;
    init_scaling = 0;

    rem_type = 0;

    confined    = 0;
    nfwall      = 0;
    wall_length = 0.0;
  // default values of sgm and eps are from T Werder et al., Nano Letters,1,12,697-702,2001
    sgm_wall    = 3.19;
    eps_wall    = 0.3135;
    rho_wall    = 0.3816;

    adjustor        = 0;
    adjust_interval = 0;
  }

  void read(std::istream &is){
    for(std::string line; std::getline(is, line);){
      std::istringstream ss(line);
      std::string tag;
      ss >> tag;
      if(ss.fail()) continue;
      if(!isalpha(tag[0])) continue; // check whether first string is alphabet or not
#define READ(name)                      \
      if(tag == std::string(#name)){    \
	ss >> name;			\
	continue;			\
      }
      READ(input_prefix);
      READ(output_prefix);
      READ(gro_in);READ(gro_out);
      READ(trr_in);READ(trr_out);
      READ(chk_in);READ(chk_out);
      READ(tkw_in);
      READ(xyz_in);
      READ(cdv_in);
      READ(cdv_out);
      READ(ene_out);
      READ(mtl_in);

      READ(nmol);READ(natom);
      READ(nsite);
      READ(dt);
      READ(rcut);
      READ(kcut);READ(alpha);
      READ(rswitch);
      READ(temperature);READ(tconstant);READ(tstat_mass);
      READ(pressure);READ(pconstant);READ(bstat_mass);

      //READ(nstep_precalc);
      READ(ninterval);READ(interval);
      READ(energy_interval);
      READ(snap_interval);
      READ(chk_interval);
      READ(mom_cancel_interval);

      READ(ngpus);READ(nprocs);READ(nthreads);
      READ(gpu_offset);
      READ(nreplica);
      READ(ntemp);READ(npress);
      READ(temp_max);READ(temp_min);
      READ(press_max);READ(press_min);

      READ(energy_max);READ(energy_min);
      READ(volume_max);READ(volume_min);

      READ(gen_vel);
      READ(cond_mode);
      READ(restart);
      READ(adjust_center);
      READ(init_scaling);

      READ(rem_type);

      READ(confined);
      READ(nfwall);
      READ(wall_length);
      READ(sgm_wall);
      READ(eps_wall);
      READ(rho_wall);

      READ(adjustor);
      READ(adjust_interval);
#undef READ
    }
  }
  void print(std::ostream &os) const {
#define PRINT(x) std::cout << #x << " " << x << std::endl
      PRINT(gro_in);PRINT(gro_out);
      PRINT(trr_in);PRINT(trr_out);
      PRINT(chk_in);PRINT(chk_out);
      PRINT(tkw_in);
      PRINT(xyz_in);
      PRINT(cdv_in);
      PRINT(cdv_out);
      PRINT(ene_out);
      PRINT(mtl_in);
      PRINT(nmol);PRINT(natom);
      PRINT(nsite);
      PRINT(dt);
      PRINT(rcut);
      PRINT(kcut);PRINT(alpha);
      PRINT(rswitch);
      PRINT(tconstant);PRINT(pconstant);
      PRINT(tstat_mass);PRINT(bstat_mass);
      PRINT(temperature);PRINT(pressure);
      //PRINT(nstep_precalc);
      PRINT(ninterval);PRINT(interval);
      PRINT(energy_interval);
      PRINT(snap_interval);
      PRINT(chk_interval);
      PRINT(mom_cancel_interval);

      PRINT(ngpus);PRINT(nprocs);PRINT(nthreads);
      PRINT(gpu_offset);
      PRINT(nreplica);
      PRINT(ntemp);PRINT(npress);
      PRINT(temp_max);PRINT(temp_min);
      PRINT(press_max);PRINT(press_min);

      PRINT(energy_max);PRINT(energy_min);
      PRINT(volume_max);PRINT(volume_min);

      PRINT(gen_vel);
      PRINT(cond_mode);
      PRINT(restart);
      PRINT(adjust_center);
      PRINT(init_scaling);

      PRINT(rem_type);

      PRINT(confined);
      PRINT(nfwall);
      PRINT(wall_length);
      PRINT(sgm_wall);
      PRINT(eps_wall);
      PRINT(rho_wall);

      PRINT(adjustor);
      PRINT(adjust_interval);
#undef PRINT
  }
};

//#define TEST_MAIN
#ifdef TEST_MAIN

int main(){
  Parameter param;
  param.read(std::cin);
  std::cout << std::endl;
  param.print(std::cout);
  return 0;
}
#endif
