#include "wham.h"

WHAM::WHAM(REMDInfo *remdinfo){
  SAFE_MALLOC(dens_state,sizeof(D*)*NUM_VOL,D*);
  for(int i=0;i<NUM_VOL;++i){
    SAFE_MALLOC(dens_state[i],sizeof(D)*NUM_POT,D);
  }
  for(int i=0;i<NUM_VOL;++i){
    for(int j=0;j<NUM_POT;++j){
      dens_state[i][j] = 0.;
  }}

  const I nreplica = remdinfo->GetNumReplica();
  SAFE_MALLOC(g,sizeof(D)*nreplica,D);
  for(int i=0;i<nreplica;++i){ g[i]=0.;}

  dos_min = 0.;
};

static inline double sum_log(double log_a,double log_b){
  if((log_a==0.)||(log_b==0.)) return log_a+log_b;
  if(log_a>log_b)
    return log_a + log(1.e0 + exp(log_b-log_a) );
  return log_b + log(1.e0 + exp(log_a-log_b) );
};

void WHAM::CalcDensState(REMDInfo *remdinfo){
  //printf("====WHAM::CalcDensState====\n");
  const I nmol     = remdinfo->GetNmol();
  const I nreplica = remdinfo->GetNumReplica();
  const I nproc    = remdinfo->GetNumProc();
  //printf("nproc=%d\n",nproc);

#pragma omp parallel for
  for(int i=0;i<NUM_VOL;++i){
    const D v = remdinfo->GetAbsoluteVolume(i);
    for(int j=0;j<NUM_POT;++j){
      const D p = remdinfo->GetAbsoluteEnergy(j);
      D denominator=0.,numerator=0.;
      for(int rep=0;rep<nreplica;++rep){
	const I hst = remdinfo->GetHistogram(rep,i,j);
	if(hst>0){
	  const I ind   = remdinfo->GetIndex(rep);
	  const D press = remdinfo->GetPressure(ind);
	  const D temp  = remdinfo->GetTemperature(ind);
	  const D h     = (p+press*v)*nmol;
	  const D bf    = h/(temp);
	  const I sum   = remdinfo->GetSumHist(rep);
	  numerator     = sum_log(log((D)hst),numerator);
	  denominator   = sum_log(log((D)sum)-bf+g[rep],denominator);
	}
      }//roop of rep
      //calc log of density state
      if(numerator!=0){
	if(pot_max < j || pot_max == 0) pot_max = j;
	if(pot_min > j || pot_min == 0) pot_min = j;
	if(vol_max < i || vol_max == 0) vol_max = i;
	if(vol_min > i || vol_min == 0) vol_min = i;
	dens_state[i][j] = numerator - denominator;
      }
    }//roop of j
  }//roop of i

  //set minimum density of state
#pragma omp parallel for
  for(int i=0;i<NUM_VOL;++i){
    for(int j=0;j<NUM_POT;++j){
      if((dos_min==0.)||(dos_min>dens_state[i][j])){
	dos_min = dens_state[i][j];
      }
  }}
#if 0
  //* // if min of dens_state < zero, calculation will be broken ?
#pragma omp parallel for
  for(int i=0;i<NUM_VOL;++i){
    for(int j=0;j<NUM_POT;++j){
      if(dens_state[i][j]!=0.){
	dens_state[i][j] -= ds_min;
      }
  }}
#endif
};

void WHAM::CalcG(REMDInfo *remdinfo){
  //printf("====WHAM::CalcG====\n");
  const I nmol = remdinfo->GetNmol();
  const I nreplica = remdinfo->GetNumReplica();
  const I nproc= remdinfo->GetNumProc();

#pragma omp parallel for
  for(int rep=0;rep<nreplica;++rep){
    const I ind   = remdinfo->GetIndex(rep);
    const D press = remdinfo->GetPressure(ind);
    const D temp  = remdinfo->GetTemperature(ind);
    D tmp=0.;
    for(int i=vol_min;i<=vol_max;++i){
      for(int j=pot_min;j<=pot_max;++j){
	const D p     = remdinfo->GetAbsoluteEnergy(j);
	const D v     = remdinfo->GetAbsoluteVolume(i);
	const D h     = (p + press*v)*nmol;
	const D bf    = h/temp;
	if(dens_state[i][j]!=0){
	  tmp = sum_log(tmp,dens_state[i][j]-bf);
	}
      }//roop of j
    }//roop of i
    g[rep] = -tmp;
  }//roop of rep
};

void WHAM::CalcPhysValue(REMDInfo *remdinfo){
  //printf("====WHAM::CalcPhysValue====\n");
  const I nmol = remdinfo->GetNmol();
  const I nproc= remdinfo->GetNumProc();
  const D Tmax = remdinfo->GetTemperatureMax();
  const D Tmin = remdinfo->GetTemperatureMin();
  const D DeltaT = (Tmax -Tmin)/(D)(NumTemp-1);
  //printf("Tmax=%lf,Tmin=%lf,DeltaT=%lf\n",Tmax,Tmin,DeltaT);

  const D Pmax = remdinfo->GetPressureMax();
  const D Pmin = remdinfo->GetPressureMin();
  const D DeltaP = (Pmax -Pmin)/(D)(NumPress-1);
  //printf("Pmax=%lf,Pmin=%lf,DeltaP=%lf\n",Pmax,Pmin,DeltaP);

#pragma omp parallel for
  for(int p=0;p<NumPress;++p){
    const D press = Pmin + DeltaP*p;
    const D p_min = remdinfo->GetEnergyMin()*nmol;
    const D v_min = remdinfo->GetVolumeMin()*nmol;
    const D h_min = p_min + press*v_min;
    //printf("min of p = %lf,v =%lf,h = %lf\n",p_min,v_min,h_min);
    for(int t=0;t<NumTemp;++t){
      const D temp = (Tmin + DeltaT*t);
      const D beta = 1./(temp);
      /*
      // set pressure min
      D prs_min = 0.0;
      for(int i=0;i<NUM_VOL;++i){
	for(int j=0;j<NUM_POT;++j){
	  if(dens_state[i][j]!=0){
	    const D prs_tmp = remdinfo->GetPressAve(i,j);
	    if(prs_tmp<prs_min){ prs_min = prs_tmp;}
	  }
      }}
      //*/
      // set minimum of averages
      Average ave_min;
      ave_min.flush();
      for(int i=0;i<NUM_VOL;++i){
	for(int j=0;j<NUM_POT;++j){
	  if(dens_state[i][j]!=0.){
	    for(int a=0;a<Average::NUM_ELEM;a++){
	      Average ave_tmp = remdinfo->GetAverages(i,j);
	      if(ave_tmp.sum[a] < ave_min.sum[a]){ ave_min.sum[a] = ave_tmp.sum[a];}
	    }
	  }
	}
      }

      //start calculation physical value
      D log_ent=0.0,log_ent_sq=0.0,pf=0.0;
      D log_pot=0.0,log_vol=0.0;
      //D log_press=0.0,log_tmp=0.0;
      Average log_ave; log_ave.flush();
      for(int i=vol_min;i<=vol_max;++i){
	for(int j=pot_min;j<=pot_max;++j){
	  if(dens_state[i][j]!=0.){
	    const D p   = remdinfo->GetRelativeEnergy(j)*nmol;
	    const D v   = remdinfo->GetRelativeVolume(i)*nmol;
	    const D h   = p + press*v;
	    const D bf  = (h+h_min)*beta;
	    //const D prs = remdinfo->GetPressAve(i,j) - prs_min;
	    //const D tmp = remdinfo->GetTempAve(i,j);
	    const Average ave = remdinfo->GetAverages(i,j);
	    //if(prs<0) printf("prs= %lf,prs_min=%lf\n",prs,prs_min);
	    log_ent    = sum_log(log_ent,   log(h)                   + dens_state[i][j] - bf);
	    log_ent_sq = sum_log(log_ent_sq,log((h+h_min)*(h+h_min)) + dens_state[i][j] - bf);
	    pf         = sum_log(pf,                                   dens_state[i][j] - bf);
	    log_pot    = sum_log(log_pot,   log(p)                   + dens_state[i][j] - bf);
	    log_vol    = sum_log(log_vol,   log(v)                   + dens_state[i][j] - bf);
	    //log_press  = sum_log(log_press, log(prs)                 + dens_state[i][j] - bf);
	    //log_tmp    = sum_log(log_tmp,   log(tmp)                 + dens_state[i][j] - bf);
	    for(int a=0;a<Average::NUM_ELEM;a++){
	      log_ave.sum[a] = sum_log(log_ave.sum[a],   log(ave.sum[a]-ave_min.sum[a]) + dens_state[i][j] - bf);
	    }
	  }
	}// roop of j
      }// roop of i
      //printf("log of ent= %lf,pf= %lf,sub = %lf\n",log_ent,pf,log_ent-pf);
      log_ent    -= pf;
      log_ent_sq -= pf;
      log_pot    -= pf;
      log_vol    -= pf;
      //log_press  -= pf;
      //log_tmp    -= pf;
      log_ave    -= pf;

      ent[p][t]       = (exp(log_ent) + h_min);
      hc[p][t]        = (exp(log_ent_sq) - ent[p][t]*ent[p][t]) / (temp*temp);
      fe[p][t]          = -temp*pf;
      potential[p][t]   = exp(log_pot)   + p_min;
      volume[p][t]      = exp(log_vol)   + v_min;
      //pressure[p][t]    = exp(log_press) + prs_min;
      //temperature[p][t] = exp(log_tmp);

      for(int a=0;a<Average::NUM_ELEM;a++){
	average[p][t].sum[a] = (exp(log_ave.sum[a]) + ave_min.sum[a]) / (D)nmol;
      }
      average[p][t].sum[Average::T]  *= (D)nmol;
      average[p][t].sum[Average::Px] *= (D)nmol;
      average[p][t].sum[Average::Py] *= (D)nmol;
      average[p][t].sum[Average::Pz] *= (D)nmol;

      ent[p][t]       /= (D)nmol;
      hc[p][t]        /= (D)nmol;
      fe[p][t]        /= (D)nmol;
      potential[p][t] /= (D)nmol;
      volume[p][t]    /= (D)nmol;
    }// roop of t
  }// roop of p
};

void WHAM::Output(REMDInfo *remdinfo){
  //printf("====WHAM::Output====\n");
  const D Tmax = remdinfo->GetTemperatureMax();
  const D Tmin = remdinfo->GetTemperatureMin();
  const D DeltaT = (Tmax -Tmin)/(D)NumTemp;

  const D Pmax = remdinfo->GetPressureMax();
  const D Pmin = remdinfo->GetPressureMin();
  const D DeltaP = (Pmax -Pmin)/(D)NumPress;

  const char *output_dir = remdinfo->GetOutputDir();
  char filename[256];

  //output physical value
  sprintf(filename,"%s/phys_value.dat",output_dir);
  FILE *fp;
  SAFE_FILEOPEN(fp,filename,"w");
  fprintf(fp,"#1pressure,2tempetature,3enthalpy,4heat capacity,5free energy,6potential,7volume,8average_pressure,9average_temperature");
  for(int i=0;i<Average::NUM_ELEM;i++){
    fprintf(fp,"%d%s",i,Average::name(i));
  }
  fprintf(fp,"\n");
  for(int p=0;p<NumPress;++p){
    const D press = Pmin + DeltaP*p;
    for(int t=0;t<NumTemp;++t){
      const D temp = Tmin + DeltaT*t;
      fprintf(fp," %e %e %e %e %e %e %e %e %e",press,temp,ent[p][t],hc[p][t],fe[p][t],potential[p][t],volume[p][t],pressure[p][t],temperature[p][t]);
      for(int a=0;a<Average::NUM_ELEM;a++){
	fprintf(fp," %e",average[p][t].sum[a]);
      }
      fprintf(fp,"\n");
    }
    fprintf(fp,"\n\n");
    if(DeltaP == 0.) break;
  }
  fclose(fp);

  //output density of state
  sprintf(filename,"%s/dens_state.dat",output_dir);
  SAFE_FILEOPEN(fp,filename,"w");
  for(int i=vol_min;i<=vol_max;++i){
    const D v = remdinfo->GetAbsoluteVolume(i);
    for(int j=pot_min;j<=pot_max;++j){
      const D p = remdinfo->GetAbsoluteEnergy(j);
      if(dens_state[i][j]!=0){
	fprintf(fp,"%lf %lf %e\n",v,p,dens_state[i][j]);
      }else{
	fprintf(fp,"%lf %lf %e\n",v,p,dos_min);
      }
    }
    fprintf(fp,"\n");
  }
  fclose(fp);

  //output fe_surface of each condition
  const int nreplica = remdinfo->GetNumReplica();
  const int nmol = remdinfo->GetNmol();
  for(int rep=0;rep<nreplica;rep++){
    const I ind   = remdinfo->GetIndex(rep);
    const D press = remdinfo->GetPressure(ind);
    const D temp  = remdinfo->GetTemperature(ind);
    const char *output_dir = remdinfo->GetOutputDir();
    sprintf(filename,"%s/fe_surface_T%lfP%lf.dat",output_dir,temp,press);
    SAFE_FILEOPEN(fp,filename,"w");
    for(int i=vol_min;i<vol_max;++i){
      const D v  = remdinfo->GetAbsoluteVolume(i);
      for(int j=pot_min;j<pot_max;++j){
	const D p  = remdinfo->GetAbsoluteEnergy(j);
	if(dens_state[i][j]!=0.){
	  const D bf = (p + press*v)/temp*(D)nmol;
	  fprintf(fp,"%lf %lf %lf\n",v,p,(-dens_state[i][j]+bf)/(D)nmol);
	}else{
	  fprintf(fp,"%lf %lf %lf\n",v,p,0.0);
	}
      }
      fprintf(fp,"\n");
    }
    fclose(fp);
  }
};

