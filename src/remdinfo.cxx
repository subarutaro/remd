#include "remdinfo.h"

#include <cassert>

#ifdef ENABLE_MPI
#include <mpi.h>
#endif
//===================== class Histogram ========================
Histogram::Histogram(int _vnum,int _pnum){
  vnum = _vnum;
  pnum = _pnum;
  histogram = new int*[vnum];
  for(int i=0;i<vnum;++i) histogram[i] = new int[pnum];
  for(int i=0;i<vnum;++i){
    for(int j=0;j<pnum;++j){
      histogram[i][j] = 0;
  }}

  ave = new Average*[vnum];
  for(int i=0;i<vnum;++i) ave[i] = new Average[pnum];
  for(int i=0;i<vnum;++i){
    for(int j=0;j<pnum;++j){
      ave[i][j].flush();
  }}

  sum = 0;
};

Histogram::~Histogram(){
  printf("===destruct Histogram===\n");
  for(int i=0;i<vnum;++i) delete[] histogram[i];
  delete[] histogram;

  for(int i=0;i<vnum;++i) delete[] ave[i];
  delete[] ave;
};

void Histogram::Increment(int v,int p){
  if(sum==0){
    vmax = vmin = v;
    pmax = pmin = p;
  }

  if(v < vmin) vmin = v;
  if(v > vmax) vmax = v;
  if(p < pmin) pmin = p;
  if(p > pmax) pmax = p;

  histogram[v][p]++;
  sum++;
};

//===================== class ExchangeList =======================
ExchangeList::ExchangeList(int in_dim_temp,int in_dim_press){
  printf("====construct ExchangeList====\n");
  dim_temp  = in_dim_temp;
  dim_press = in_dim_press;
  printf("dim_temp=%d,dim_press=%d\n",dim_temp,dim_press);
  list = new Pair*[4];
  for(int i=0;i<4;++i) list[i] = new Pair[dim_temp*dim_press/2];
  rate = new AcceptRate*[4];
  for(int i=0;i<4;++i) rate[i] = new AcceptRate[dim_temp*dim_press/2];

  //setting pairlist 0
  int count0 = 0;
  for(int j=0;j<dim_press;++j){
    for(int i=0;i<dim_temp/2;++i){
      list[0][count0].a = 2*i   + j*dim_temp;
      list[0][count0].b = 2*i+1 + j*dim_temp;
      printf("list[0][%d]= (%d,%d)\n",count0,list[0][count0].a,list[0][count0].b);
      count0++;
    }
  }
  length[0] = count0;
  //setting pairlist 1
  int count1 = 0;
  for(int j=0;j<dim_press/2;++j){
    for(int i=0;i<dim_temp;++i){
      list[1][count1].a = i + 2*j*dim_temp;
      list[1][count1].b = i + (2*j+1)*dim_temp;
      printf("list[1][%d]= (%d,%d)\n",count1,list[1][count1].a,list[1][count1].b);
      count1++;
    }
  }
  length[1] = count1;
  //setting pairlist 2
  int count2 = 0;
  for(int j=0;j<dim_press;++j){
    for(int i=0;i<dim_temp/2-(dim_temp+1)%2;++i){
      list[2][count2].a = 2*i   + j*dim_temp +1;
      list[2][count2].b = 2*i+1 + j*dim_temp +1;
      printf("list[2][%d]= (%d,%d)\n",count2,list[2][count2].a,list[2][count2].b);
      count2++;
    }
  }
  length[2] = count2;
  //setting pairlist 3
  int count3 = 0;
  for(int j=0;j<dim_press/2-(dim_press+1)%2;++j){
    for(int i=0;i<dim_temp;++i){
      list[3][count3].a = i + (2*j+1)*dim_temp;
      list[3][count3].b = i + (2*j+2)*dim_temp;
      printf("list[3][%d]= (%d,%d)\n",count3,list[3][count3].a,list[3][count3].b);
      count3++;
    }
  }
  length[3] = count3;
};

ExchangeList::~ExchangeList(){
  printf("===destruct ExchangeList===\n");
  for(int i=0;i<4;++i) delete[] list[i];
  delete[] list;
  for(int i=0;i<4;++i) delete[] rate[i];
  delete[] rate;
};

//===================== class REMDInfo =========================

static void read(int *target,FILE *fp){
  char unko[256];
  fscanf(fp,"%s %i",unko,target);
};

static void read(unsigned long *target,FILE *fp){
  char unko[256];
  fscanf(fp,"%s %ld",unko,target);
};

static void read(double *target,FILE *fp){
  char unko[256];
  fscanf(fp,"%s %lf",unko,target);
};

static void read(char *target,FILE *fp){
  char unko[256];
  fscanf(fp,"%s %s",unko,target);
};

REMDInfo::REMDInfo(const Parameter _p){
  printf("=== REMDInfo::REMDInfo ===\n");

  nmol = _p.nmol;
  step_max = _p.ninterval*_p.interval;
  interval = _p.interval;
#ifdef ENABLE_MPI
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#else
  nproc = 1;
  rank = 0;
#endif
  ngpu  = _p.ngpus;

  sprintf(input_dir,"%s",_p.input_prefix.c_str());
  sprintf(output_dir,"%s",_p.output_prefix.c_str());

  nreplica_global = _p.nreplica;
  assert(nreplica_global % nproc == 0);
  nreplica = _p.nreplica / nproc;

  dim_temp = _p.ntemp;
  dim_press = _p.npress;
  if(nreplica_global != (dim_temp*dim_press)){
    fprintf(stderr,"error: # of replicas must be equal to DimTemp*DimPress\n");
    exit(EXIT_FAILURE);
  }
  temperature_max = _p.temp_max  / unit_temp; // [K]
  temperature_min = _p.temp_min  / unit_temp;
  pressure_max    = _p.press_max / unit_press;// [Pa]
  pressure_min    = _p.press_min / unit_press;
#ifdef LOG_SAMPLING
  energy_max = _p.energy_max; // [kJ/mol]
  energy_min = _p.energy_min;
  delta_energy = log(energy_max - energy_min)/(D)NUM_POT;
  volume_max = _p.volume_max; // [A^3/mol]
  volume_min = _p.volume_min;
  delta_volume = log(volume_max - volume_min)/(D)NUM_VOL;
#else
  energy_max = _p.energy_max; // [kJ/mol]
  energy_min = _p.energy_min;
  delta_energy = (energy_max - energy_min)/(D)NUM_POT;
  volume_max = _p.volume_max; // [A^3/mol]
  volume_min = _p.volume_min;
  delta_volume = (volume_max - volume_min)/(D)NUM_VOL;
#endif
  ensemble = (_p.confined<<CSHIFT) + (_p.pconstant<<PSHIFT) + (_p.tconstant<<TSHIFT);
  bkup_ninterval = _p.energy_interval;

  restart = _p.restart;

  ec_type = 0;
  if(dim_temp==1) ec_type = 1;
#ifdef ENABLE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  if(rank == 0){
    int cond_mode = _p.cond_mode;
    histogram = new Histogram*[nreplica_global];
    tc = new TunnelCount*[nreplica_global];
    for(int i=0;i<nreplica_global;++i){
      histogram[i] = new Histogram(NUM_VOL,NUM_POT);
      tc[i] = new TunnelCount();
    }

    pairlist = new ExchangeList(dim_temp,dim_press);

    index       = new int[nreplica_global];
    temperature = new D[nreplica_global];
    pressure    = new D[nreplica_global];
    isExchanged = new bool[nreplica_global];
    for(int i=0;i<nreplica_global;++i) index[i]=i;
    for(int i=0;i<nreplica_global;++i) isExchanged[i]=false;

    assert(dim_temp * dim_press == nreplica_global);
    if(cond_mode == 0){
      SetConditionGeometrical(dim_temp,dim_press);
    }else if(cond_mode == 1){
      char cond_file[256];
      sprintf(cond_file,"%s/phys_value.dat",input_dir);
      SetConditionGeometrical(dim_temp,dim_press);
      SetConditionFromHeatCapacity(cond_file,dim_temp,dim_press,true);
    }else if(cond_mode == 2){
      char cond_file[256];
      sprintf(cond_file,"%s/condition.dat",input_dir);
      SetConditionGeometrical(dim_temp,dim_press);
      SetConditionFromFile(cond_file,dim_temp,dim_press);
    }else if(cond_mode == 3){
      char cond_file[256];
      SetConditionArithmetic(dim_temp,dim_press);
    }else if(cond_mode == 4){
      char cond_file[256];
      sprintf(cond_file,"%s/phys_value.dat",input_dir);
      SetConditionGeometrical(dim_temp,dim_press);
      SetConditionFromHeatCapacity(cond_file,dim_temp,dim_press,false);
    }else{
      std::cerr << "error: set cond_mode 0, 1 or 2" << std::endl;
      exit(EXIT_FAILURE);
    }
  }else{
    printf("%d: receive index, temperature and pressure\n",rank);
    index       = new int[nreplica];
    temperature = new D[nreplica];
    pressure    = new D[nreplica];
    isExchanged = new bool[nreplica];
  }
  BroadcastConditionAndIndex();

  if(rank == 0){
    printf("Conditions at master rank:\n");
    for(int i=0;i<nreplica_global;i++){
      printf("%d %d %lf %lf\n",i,index[i],temperature[i],pressure[i]);
    }
  }
#ifdef ENABLE_MPI
  for(int i=0;i<nproc;i++){
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == i){
      if(i==0) printf("Conditions at each rank:\n");
      for(int i=0;i<nreplica;i++){
	printf("%d: %d %d %lf %lf\n",rank,i,index[i],temperature[i],pressure[i]);
      }
    }
  }
#endif
}

void REMDInfo::BroadcastConditionAndIndex(){
#ifdef ENABLE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank == 0){
    // broadcast index, temperature, pressure
    //printf("0: broadcast index, temperature and pressure\n");
    MPI_Request req[nproc];
    MPI_Status  stat[nproc];
    for(int i=1; i<nproc;i++) MPI_Isend(index      +i*nreplica,nreplica,MPI_INT,   i,i+0*nproc,MPI_COMM_WORLD,&req[i-1]);
    MPI_Waitall(nproc-1,req,stat);
    for(int i=1; i<nproc;i++) MPI_Isend(temperature+i*nreplica,nreplica,MPI_DOUBLE,i,i+1*nproc,MPI_COMM_WORLD,&req[i-1]);
    MPI_Waitall(nproc-1,req,stat);
    for(int i=1; i<nproc;i++) MPI_Isend(pressure   +i*nreplica,nreplica,MPI_DOUBLE,i,i+2*nproc,MPI_COMM_WORLD,&req[i-1]);
    MPI_Waitall(nproc-1,req,stat);
    for(int i=1; i<nproc;i++) MPI_Isend(isExchanged+i*nreplica,nreplica,MPI_C_BOOL,i,i+3*nproc,MPI_COMM_WORLD,&req[i-1]);
    MPI_Waitall(nproc-1,req,stat);
 }else{
    //printf("%d: receive index, temperature and pressure\n",rank);
    // recv index, temperature, pressure
    MPI_Recv(index,      nreplica,MPI_INT,   0,rank+0*nproc,MPI_COMM_WORLD,NULL);
    MPI_Recv(temperature,nreplica,MPI_DOUBLE,0,rank+1*nproc,MPI_COMM_WORLD,NULL);
    MPI_Recv(pressure,   nreplica,MPI_DOUBLE,0,rank+2*nproc,MPI_COMM_WORLD,NULL);
    MPI_Recv(isExchanged,nreplica,MPI_C_BOOL,0,rank+3*nproc,MPI_COMM_WORLD,NULL);
  }
#endif
}

REMDInfo::~REMDInfo(){
  printf("===destruct REMDInfo===\n");
  delete[] pressure;
  delete[] temperature;
  delete[] index;

  delete pairlist;
  if(rank == 0){
    for(int i=0;i<nreplica;++i) delete[] histogram[i];
    delete[] histogram;
  }

};

void REMDInfo::ReadTemperatureFromFile(char *filename){
  FILE *fp=fopen(filename,"r");
  int i;
  double d;
  while(fscanf(fp,"%d %lf",&i,&d)!=EOF){
    temperature[i]=d;
    printf("temperature[%d]= %lf\n",i,temperature[i]);
  }
};

void REMDInfo::ShowAll(){
  printf("====REMDInfo::ShowAll()====\n");
  printf("step_max= %d\n",step_max);
  printf("interval= %d\n",interval);

  printf("nreplica= %d\n",nreplica);
  printf("temperature_max= %lf\n",temperature_max);
  printf("temperature_min= %lf\n",temperature_min);
  printf("pressure_max= %lf\n",pressure_max);
  printf("pressure_min= %lf\n",pressure_min);
  printf("delta_energy= %lf\n",delta_energy);
  printf("energy_max= %lf\n",energy_max);
  printf("energy_min= %lf\n",energy_min);
  printf("delta_volume= %lf\n",delta_volume);
  printf("volume_max= %lf\n",volume_max);
  printf("volume_min= %lf\n",volume_min);

  printf("output_dir= %s\n",output_dir);

  printf("nproc= %d\n",nproc);
  printf("ngpu= %d\n",ngpu);

  printf("mode = %d\n",mode);

  for(int rep=0;rep<nreplica;++rep){
    printf("rep%d T=%lf,P=%lf\n",rep,temperature[rep],pressure[rep]);
  }
};

void REMDInfo::SetConditionGeometrical(int dim_temp,int dim_press){
  D t = (dim_temp==1)  ? 1.0 : powf(temperature_max/temperature_min,1./(D)(dim_temp-1));
  D p = (dim_press==1) ? 0.0 : ((pressure_max - pressure_min)/(D)(dim_press-1));
  for(int j=0;j<dim_press;++j){
    for(int i=0;i<dim_temp;++i){
      int ind = i+j*dim_temp;
      temperature[ind] = temperature_min*powf(t,i);
      pressure[ind] = pressure_min + p * j;
      //printf(" md[%d]: T= %lf,P= %lf\n",ind,temperature[i+j*dim_temp],pressure[i+j*dim_temp]);
    }
  }
};

void REMDInfo::SetConditionArithmetic(int dim_temp,int dim_press){
  D t = (dim_temp==1)? 0.0 : ((temperature_max - temperature_min)/(D)(dim_temp-1));
  D p = (dim_press==1)?0.0 : ((pressure_max - pressure_min)/(D)(dim_press-1));
  for(int j=0;j<dim_press;++j){
    for(int i=0;i<dim_temp;++i){
      int ind = i+j*dim_temp;
      temperature[ind] = temperature_min + t*i;
      pressure[ind] = pressure_min + p*j;
      printf(" md[%d]: T= %lf,P= %lf\n",ind,temperature[i+j*dim_temp],pressure[i+j*dim_temp]);
    }
  }
};

void REMDInfo::SetConditionFromFile(const char* filename,int dim_temp,int dim_press){
  std::ifstream ifs(filename);
  if(ifs.fail()){
    std::cerr << "error: open " << filename << " failed" << std::endl;
    exit(EXIT_FAILURE);
  }
  int idx,count=0,pcount=0,tcount=0;
  for(std::string line; std::getline(ifs,line);){
    std::stringstream ss(line);
    ss >> idx >> pressure[count] >> temperature[count];
    std::cout << idx << " " << pressure[count] << " " << temperature[count] << std::endl;
    if(idx!=0){
      if(pressure[idx] != pressure[idx-1]){
	tcount = 0;
	pcount++;
      }
    }
    tcount++;
    count++;
  }
  if(count != nreplica){
    std::cerr << "error: # of condition != that in input file" << std::endl;
    exit(EXIT_FAILURE);
  }
}

/*
static D calc_hctmp(double *p, double *t, int pcount, int tcount, int pdim){
      if(t[tcount+1]<ttmp && tcount < ntemp-1) tcount++;
      double tmp1 = hc[pcount][tcount]   + (hc[pcount][tcount+1]   - hc[pcount][tcount])   / (t[tcount+1]-t[tcount]) * (ttmp - t[tcount]);
      if(pdim==1) return tmp1;

      double tmp2 = hc[pcount+1][tcount] + (hc[pcount+1][tcount+1] - hc[pcount+1][tcount]) / (t[tcount+1]-t[tcount]) * (ttmp - t[tcount]);
      return tmp1 + (tmp2 - tmp1) / (p[pcount+1]-p[pcount]) * (ptmp - p[pcount]);
}
//*/

void ReadHeatCapacityFromFile(const std::string filename, double *p,double *t,double **hc,int &pcount,int &tcount,const int nmol,const bool type){
  // read heat capacity from file
  std::ifstream ifs(filename.c_str());
  if(ifs.fail()){
    std::cerr << "error: open " << filename << " failed" << std::endl;
    exit(EXIT_FAILURE);
  }
  //*
  std::string comment;
  std::getline(ifs, comment);

  double ptmp0,ttmp,hctmp,dummy;
  double ptmp1;
  pcount=0,tcount=0;
  for(std::string line; std::getline(ifs, line);){
    if(line.empty()) continue;
    if(line[0]=='#') continue;
    std::istringstream ss(line);
    if(ss.fail()) continue;
    if(type) ss >> ptmp0 >> ttmp >> dummy >> hctmp;
    else     ss >> dummy >> ptmp0 >> ttmp >> dummy >> dummy >> dummy >> hctmp;
    if(pcount==0 && tcount==0) ptmp1 = ptmp0;
    if(ptmp0!=ptmp1){
      tcount = 0;
      pcount++;
    }
    p[pcount] = ptmp0;
    t[tcount] = ttmp;
    hc[pcount][tcount] = hctmp * nmol;
    //std::cout << pcount << " " << tcount << " " << ttmp << " " << ptmp0 << " " << hctmp << std::endl;
    ptmp1 = ptmp0;
    tcount++;
  }
  pcount++;
  /*
  std::cout << pcount << " " << tcount << std::endl;
  for(int j=0;j<pcount;j++)
    for(int i=0;i<tcount;i++)
      std::cout << p[j] << " " << t[i] << " " << hc[j][i] << std::endl;
  //*/
}

void REMDInfo::SetConditionFromHeatCapacity(std::string filename,int dim_temp,int dim_press,const bool type){
  std::cout << "generaging temperature condition from " << filename << std::endl;
  const int npmax = 10000, ntmax = 10000;
  double *p,*t,**hc;
  const int psize = sizeof(double)*npmax;
  const int tsize = sizeof(double)*ntmax;
  SAFE_MALLOC(p,psize,double);
  SAFE_MALLOC(t,tsize,double);
  SAFE_MALLOC(hc,psize,double*);
  for(int i=0;i<npmax;i++)  SAFE_MALLOC(hc[i],tsize,double);

  int ntemp, npress;
  ReadHeatCapacityFromFile(filename,p,t,hc,npress,ntemp,nmol,type);
  assert(ntemp  < ntmax);
  assert(npress < npmax);

  unsigned int ndiv = 100000;
#if 0
  //pressure is determined geometrically
  double pconst = pow(pressure_max/pressure_min, 1.0/(double)(dim_press - 1));
  for(int j=0;j<dim_press;j++){
    for(int i=0;i<dim_temp;i++){
      pressure[j*dim_temp + i] = pressure_min * pow(pconst,j);
  }}
#else
  //pressure is determined arithmetrically
  double pconst = (pressure_max - pressure_min) / (double)(dim_press - 1);
  for(int j=0;j<dim_press;j++){
     for(int i=0;i<dim_temp;i++){
       if(dim_press > 1) pressure[j*dim_temp + i] = pressure_min + pconst*j;
       else              pressure[j*dim_temp + i] = pressure_min;
  }}
#endif
  double dt = (temperature_max - temperature_min)/(double)(ndiv - 1);

  int pcount = 0;
  for(int j=0;j<dim_press;j++){
    double ptmp = pressure[j*dim_temp];

    double S = 0.0;
    int tcount = 0;
    double ttmp = temperature_min;
    while(p[pcount+1]< ptmp && pcount < npress-1) pcount++;
    //std::cout << "pcount " << pcount << " " << "p= " << p[pcount] << std::endl;
    while(ttmp<temperature_max){
      if(t[tcount+1]<ttmp && tcount < ntemp-1) tcount++;
      double tmp1 = hc[pcount][tcount]   + (hc[pcount][tcount+1]   - hc[pcount][tcount])   / (t[tcount+1]-t[tcount]) * (ttmp - t[tcount]);
      double tmp2 = hc[pcount+1][tcount] + (hc[pcount+1][tcount+1] - hc[pcount+1][tcount]) / (t[tcount+1]-t[tcount]) * (ttmp - t[tcount]);
      double hctmp = tmp1 + (tmp2 - tmp1) / (p[pcount+1]-p[pcount]) * (ptmp - p[pcount]);
      //std::cout << ttmp << " " << hctmp << " " << tmp1 << " " << tmp2 << " " << p[pcount] << " " << p[pcount+1] << " " << ptmp << std::endl;
      S += dt*sqrt(hctmp);
      //S += dt*hctmp;
      ttmp += dt;
    }
    double c  = pow(temperature_max/temperature_min,1.0/(double)(dim_temp-1));
    double S0 = (1.0 - c) / (1.0 - pow(c,(double)(dim_temp-1))) * S;
    double sum = 0.0;
    //std::cout << "S= " << S << " S0= " << S0 << std::endl;

    temperature[j*dim_temp] = ttmp  = temperature_min;
    tcount = 0;
    for(int i=1;i<dim_temp-1;i++){
      double Stmp = (1.0 - pow(c,i))/(1.0 - c) * S0;
      //std::cout << "Stmp= " << Stmp << std::endl;
      while(sum < Stmp){
	if(t[tcount+1]<ttmp && tcount < ntemp) tcount++;
	double tmp1 = hc[pcount][tcount] + (hc[pcount][tcount+1] - hc[pcount][tcount]) / (t[tcount+1]-t[tcount]) * (ttmp - t[tcount]);
	double tmp2 = hc[pcount+1][tcount] + (hc[pcount+1][tcount+1] - hc[pcount+1][tcount]) / (t[tcount+1]-t[tcount]) * (ttmp - t[tcount]);
	double hctmp = tmp1 + (tmp2 - tmp1) / (p[pcount+1]-p[pcount]) * (ptmp - p[pcount]);
	if(hctmp<0.0){
	  std::cout << "hctmp<0:" << hctmp << " " << tmp1 << " " << tmp2 << " " << p[pcount+1] << " " << p[pcount] << " " << ptmp << std::endl;
	  continue;
	}
	//sum  += dt * hctmp;
	sum  += dt * sqrt(hctmp);
	ttmp += dt;
      }
      //std::cout << "sum= " << sum << std::endl;
      temperature[j*dim_temp+i] = ttmp;
    }
    temperature[(j+1)*dim_temp-1] = temperature_max;
  }

  for(int i=0;i<nreplica_global;i++){
   std::cout << "rep" << i << ": " << temperature[i] << " " << pressure[i] << std::endl;
  }

  free(p);
  free(t);
  for(int i=0;i<npress;i++) free(hc[i]);
  free(hc);
}

void REMDInfo::SetConditionFromHeatCapacity2(std::string filename,int dim_temp,int dim_press){
  std::cout << "generaging temperature condition from " << filename << std::endl;
  const double ar = log(0.6);

  const int npmax = 1000, ntmax = 1000;
  double *p,*t,**hc;
  const int psize = sizeof(double)*npmax;
  const int tsize = sizeof(double)*ntmax;
  SAFE_MALLOC(p,psize,double);
  SAFE_MALLOC(t,tsize,double);
  SAFE_MALLOC(hc,psize,double*);
  for(int i=0;i<npmax;i++)  SAFE_MALLOC(hc[i],tsize,double);

  int ntemp, npress;
  ReadHeatCapacityFromFile(filename,p,t,hc,npress,ntemp,nmol,true);

  unsigned int ndiv = 100000;
  //pressure is determined geometrically
  double pconst = pow(pressure_max/pressure_min, 1.0/(double)(dim_press - 1));
  for(int j=0;j<dim_press;j++){
    for(int i=0;i<dim_temp;i++){
      pressure[j*dim_temp + i] = pressure_min * pow(pconst,j);
  }}

  double dt = (temperature_max - temperature_min)/(double)(ndiv - 1);

  int pcount = 0;
  for(int j=0;j<dim_press;j++){
    double ptmp = pressure[j*dim_temp];

    double S = 0.;
    int tcount = 0;
    double ttmp = temperature_min;
    double delta = 0;
    while(p[pcount+1]< ptmp && pcount < npress-1) pcount++;
    //std::cout << "pcount " << pcount << " " << "p= " << p[pcount] << std::endl;
    for(int i=1;i<dim_temp;i++){
      if(t[tcount+1]<ttmp && tcount < ntemp-1) tcount++;
      while(S < - ar / (1./(ttmp+delta) - 1./ttmp)){
	double tmp1 = hc[pcount][tcount]   + (hc[pcount][tcount+1]   - hc[pcount][tcount])   / (t[tcount+1]-t[tcount]) * (ttmp + delta - t[tcount]);
	double tmp2 = hc[pcount+1][tcount] + (hc[pcount+1][tcount+1] - hc[pcount+1][tcount]) / (t[tcount+1]-t[tcount]) * (ttmp + delta - t[tcount]);
	double hctmp = tmp1 + (tmp2 - tmp1) / (p[pcount+1]-p[pcount]) * (ptmp - p[pcount]);
	//S += dt*sqrt(hctmp);
	S += dt*hctmp;
	delta += dt;
      }
      ttmp += delta;
      temperature[j*dim_temp + i] = ttmp;
    }
  }

  for(int i=0;i<nreplica_global;i++){
   std::cout << "rep" << i << ": " << temperature[i] << " " << pressure[i] << std::endl;
  }

  free(p);
  free(t);
  for(int i=0;i<npress;i++) free(hc[i]);
  free(hc);
}

bool REMDInfo::ConditionAdjustor(){
  std::cout << "warning: coudition adjustor only supports 1D REM for different temperatures. pressures won't be adjusted" << std::endl;
  int npair = 0;
  for(int type=0;type<4;++type) npair += pairlist->GetLength(type);

  double dT[npair],rate[npair];
  double ave = 0.0;
  double delta;
  int count = 0;
  for(int type=0;type<4;++type){
    for(int i=0;i<pairlist->GetLength(type);++i){
      const Pair pair = GetPair(type,i);
      const int ind1 = GetIndex(pair.a);
      const int ind2 = GetIndex(pair.b);
      rate[count] = pairlist->GetRate(type,i);
      ave += rate[count];
      dT[count] = temperature[ind2] - temperature[ind1];

      std::cout << temperature[ind1] << " " << temperature[ind2] << " : " << rate[count] << std::endl;

      if(rate[count]==0.){
	std::cerr << "temperature cant be adusted if accpet rate is 0" << std::endl;
	return true;
      }
      if(count==0) delta = dT[count];
      if(delta>dT[count]) delta = dT[count];
      count++;
    }
  }

  for(int rep=0;rep<nreplica;rep++){
    histogram[rep]->flush();
  }
  pairlist->FlushRate();

  ave /= (double)(npair);
  std::cout << "average rate is " << ave << std::endl;
  // check adjusting is enough or not
  bool isEnough = true;
  for(int i=0;i<npair;i++){
    if(fabs(rate[i] - ave)>0.05) isEnough = false;
    std::cout << "|rate - ave| = " << fabs(rate[i] - ave) << std::endl;
  }
  if(isEnough) return false;

  delta *= 0.1; // delta is half of minimum of dT
  for(int rep=0;rep<npair;rep++)
    dT[rep] += (rate[rep]-ave)/ave * delta;

  double max_len = 0;
  for(int type=0;type<4;++type)
    if(max_len<pairlist->GetLength(type)) max_len = pairlist->GetLength(type);

  count = 0;
  for(int i=0;i<max_len;++i){
    for(int type=0;type<4;++type){
      if(i < pairlist->GetLength(type)){
	const Pair pair = GetPair(type,i);
	const int ind1 = GetIndex(pair.a);
	const int ind2 = GetIndex(pair.b);
	temperature[ind2] = temperature[ind1] + dT[count++];
      }
    }
  }

  std::cout << "temperatures are adjusted." << std::endl;
  for(int rep=0;rep<nreplica;rep++){
    int ind = index[rep];
    std::cout <<  rep << " " << temperature[ind] << std::endl;
  }

  return true;
}

void REMDInfo::ReadHistogramFromBackUp(){
  for(int rep=0;rep<nreplica;++rep){
    char filename[256];
    //remdinfo has conditions of rep th replica in order (T_i can be less than T_rep+1 in remdinfo)
    //you have to use index to get conditions of rep th condition in order
    sprintf(filename,"%s/histo_P%lf_T%lf.dat",input_dir,pressure[index[rep]],temperature[index[rep]]);
    std::istream *is;
    is = new std::ifstream(filename,std::ios_base::out);
    if(is->fail()){
      std::cerr << "error: input histogram file " << filename << " is not found" << std::endl;
      exit(EXIT_FAILURE);
    }

    std::string dummy;
    int vmax,vmin,pmax,pmin,hist;
    *is >> dummy >> dummy >> dummy >> dummy;
    *is >> dummy >> vmax >> vmin >> pmax >> pmin;
    std::cout << "reading " << filename << std::endl;
    std::cout << vmax << " " << vmin << " " << pmax << " " << pmin << std::endl;

    histogram[rep]->SetVolMax(vmax);
    histogram[rep]->SetVolMin(vmin);
    histogram[rep]->SetPotMax(pmax);
    histogram[rep]->SetPotMin(pmin);

    int sum = 0;
    for(int v=vmin;v<=vmax;v++){
      for(int p=pmin;p<=pmax;p++){
	*is >> dummy >> dummy >> hist;
	histogram[rep]->SetHist(v,p,hist);
	sum += hist;
	Average tmp;tmp.flush();
	for(int i=0;i<Average::NUM_ELEM;i++){
	  *is >> tmp.sum[i];
	}
	histogram[rep]->SetAve(v,p,tmp);
      }
    }
    histogram[rep]->SetSum(sum);
  }
}

void REMDInfo::OutputHistogram(){
  for(int rep=0;rep<nreplica;++rep){
    char filename[256];
    //remdinfo has conditions of rep th replica in order (T_i can be less than T_rep+1 in remdinfo)
    //you have to use index to get conditions of rep th condition in order
    sprintf(filename,"%s/histo_P%lf_T%lf.dat",output_dir,pressure[index[rep]],temperature[index[rep]]);
    FILE *fp;
    SAFE_FILEOPEN(fp,filename,"w");
    const int vmax = histogram[rep]->GetVolMax();
    const int vmin = histogram[rep]->GetVolMin();
    const int pmax = histogram[rep]->GetPotMax();
    const int pmin = histogram[rep]->GetPotMin();
    fprintf(fp,"#vmax vmin pmax pmin\n");
    fprintf(fp,"# %d %d %d %d\n",vmax,vmin,pmax,pmin);
#ifdef LOG_SAMPLING
    for(int v=vmin;v<=vmax;v++){
      D vol = volume_min + exp(v*delta_volume);
      for(int p=pmin;p<=pmax;p++){
	D pot = energy_min + exp(p*delta_energy);
#else
    for(int v=vmin;v<=vmax;v++){
      D vol = volume_min + v*delta_volume;
      for(int p=pmin;p<=pmax;p++){
	D pot = energy_min + p*delta_energy;
#endif
	fprintf(fp," %10.6e %10.6e %10d",
		vol,pot,histogram[rep]->GetHist(v,p));
	Average tmp = histogram[rep]->GetAverages(v,p);
	for(int i=0;i<Average::NUM_ELEM;i++){
	  fprintf(fp," %10.6e",tmp.sum[i]);
	}
	fprintf(fp,"\n");
      }
      fprintf(fp,"\n");
    }
    fclose(fp);
  }
};

void REMDInfo::OutputAcceptRate(){
  FILE *fp;
  char filename[256];
  sprintf(filename,"%s/accept_rate.dat",output_dir);
  //printf("==== making %s ====\n",filename);
  SAFE_FILEOPEN(fp,filename,"w");
  fprintf(fp,"#Ti,Pi,Tj,Pj,AR\n");
#if 0
  for(int type=0;type<4;++type){
    for(int i=0;i<pairlist->GetLength(type);++i){
      const Pair pair = GetPair(type,i);
      const int ind1 = GetIndex(pair.a);
      const int ind2 = GetIndex(pair.b);
      fprintf(fp,"%lf %lf %lf %lf %lf\n",temperature[ind1],pressure[ind1],temperature[ind2],pressure[ind2],pairlist->GetRate(type,i));
    }
  }
#else
  const int ll0 = pairlist->GetLength(0);
  const int ll1 = pairlist->GetLength(1);
  const int ll2 = pairlist->GetLength(2);
  const int ll3 = pairlist->GetLength(3);

  for(int p=0;p<dim_press;p++){
    const int i0 = p * ll0 / dim_press;
    const int i2 = p * ll2 / dim_press;
    for(int i = 0;i<ll0/dim_press;i++){
      Pair pair = GetPair(0,i+i0);
      int ind0 = GetIndex(pair.a);
      int ind1 = GetIndex(pair.b);
      fprintf(fp,"%lf %lf %lf %lf %lf\n",temperature[ind0],pressure[ind0],temperature[ind1],pressure[ind1],pairlist->GetRate(0,i+i0));
      if(i >= ll2/dim_press) break;
      pair = GetPair(2,i+i2);
      ind0 = GetIndex(pair.a);
      ind1 = GetIndex(pair.b);
      fprintf(fp,"%lf %lf %lf %lf %lf\n",temperature[ind0],pressure[ind0],temperature[ind1],pressure[ind1],pairlist->GetRate(2,i+i2));
    }
    fprintf(fp,"\n");
  }

 for(int t=0;t<dim_temp;t++){
   const int i1 = t * ll1 / dim_temp;
   const int i3 = t * ll3 / dim_temp;
   for(int i = 0;i<ll1/dim_temp;i++){
     Pair pair = GetPair(1,i+i1);
     int ind0 = GetIndex(pair.a);
     int ind1 = GetIndex(pair.b);
     fprintf(fp,"%lf %lf %lf %lf %lf\n",temperature[ind0],pressure[ind0],temperature[ind1],pressure[ind1],pairlist->GetRate(1,i+i1));
     if(i >= ll3/dim_temp) break;
     pair = GetPair(3,i+i3);
     ind0 = GetIndex(pair.a);
     ind1 = GetIndex(pair.b);
     fprintf(fp,"%lf %lf %lf %lf %lf\n",temperature[ind0],pressure[ind0],temperature[ind1],pressure[ind1],pairlist->GetRate(3,i+i3));
   }
   fprintf(fp,"\n");
  }

#endif
  fclose(fp);
};

void REMDInfo::OutputTunnelCount(){
  FILE *fp;
  char filename[256];
  sprintf(filename,"%s/tunnel_count.dat",output_dir);
  SAFE_FILEOPEN(fp,filename,"w");
  for(int rep=0;rep<nreplica;rep++){
    fprintf(fp," %d",tc[rep]->GetCount());
  }
  fprintf(fp,"\n");
  fclose(fp);
}

Average REMDInfo::GetAverages(int v,int p){
  int sum_hist=0;
  Average sum_ave;
  sum_ave.flush();
  for(int rep=0;rep<nreplica;++rep){
    const int hist = histogram[rep]->GetHist(v,p);
    Average tmp = histogram[rep]->GetAverages(v,p);
    if(hist != 0){
      sum_hist += hist;
      sum_ave  += tmp*hist;
    }
  }
  if(sum_hist>0){
    return sum_ave/(D)sum_hist;
  }else{
    sum_ave.flush();
    return sum_ave;
  }
};
