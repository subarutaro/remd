#include "bond.h"

double distance(const CDV c1, const CDV c2){
  const double dx = c1.r[0] - c2.r[0];
  const double dy = c1.r[1] - c2.r[1];
  const double dz = c1.r[2] - c2.r[2];

  return sqrt(dx*dx + dy*dy + dz*dz);
}

void make_intra_bond(const CDV c1, const CDV c2, map< pair<int,int>, int> &bond){
  if(c1.i > c2.i) return;
  if(c1.t>2 || c2.t>2) return;

  if(c1.m == c2.m){  // intra bond
    if(c1.t != c2.t){
      pair<int,int> key(c1.i,c2.i);
      bond.insert( pair< pair<int,int> ,int>(key,0) );
    }
  }
}

void make_hydrogen_bond(const vector<CDV>::iterator c1, const vector<CDV>::iterator c2, map<pair<int,int>,int> &bond){
  /*
  double pot = 0.0;

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      const double r1[3] = {(c1+i)->r[0],(c1+i)->r[1],(c1+i)->r[2]};
      const double r2[3] = {(c2+j)->r[0],(c2+j)->r[1],(c2+j)->r[2]};

      const double dx = r1[0] - r2[0];
      const double dy = r1[1] - r2[1];
      const double dz = r1[2] - r2[2];

      const double ri = rsqrt(dx*dx + dy*dy + dz*dz);
      pot += q1 * q2 * r;
    }
  }
  if(pot < 3.0){
      pair<int,int> key(c1.i,c2.i);
      bond.insert( pair< pair<int,int> ,int>(key,1) );
  }
  //*/
  if(c1->m == c2->m) return;
  if(c1->t == 0 && c2->t == 0){
    const CDV O1 = *c1;
    const CDV O2 = *c2;

    const double dOOx = O1.r[0] - O2.r[0];
    const double dOOy = O1.r[1] - O2.r[1];
    const double dOOz = O1.r[2] - O2.r[2];
    const double dOO  = sqrt(dOOx*dOOx + dOOy*dOOy + dOOz*dOOz);
    //if(dOO > 3.0) return;

    for(int i=1;i<=2;i++){
      const CDV H = *(c1+i);
      
      const double dHO1x = H.r[0] - O1.r[0];
      const double dHO1y = H.r[1] - O1.r[1];
      const double dHO1z = H.r[2] - O1.r[2];
      const double dHO1  = sqrt(dHO1x*dHO1x + dHO1y*dHO1y + dHO1z*dHO1z);
      
      const double dHO2x = H.r[0] - O2.r[0];
      const double dHO2y = H.r[1] - O2.r[1];
      const double dHO2z = H.r[2] - O2.r[2];
      const double dHO2  = sqrt(dHO2x*dHO2x + dHO2y*dHO2y + dHO2z*dHO2z);

      if(dHO2 > 2.0) continue;
      const double angle = acos( (dHO1*dHO1 + dHO2*dHO2 - dOO*dOO) / (2.0 * dHO1 * dHO2) ) * 180 / M_PI;
      if(angle > 155.0){
	pair<int,int> key(H.i,O2.i);
	bond.insert( pair< pair<int,int> ,int>(key,1) );
	return;
      }
    }
  }
}

void make_bonds(vector<CDV> cdv, map< pair<int,int>, int> &bond){
  for(vector<CDV>::iterator it0 = cdv.begin();it0 != cdv.end(); ++it0){
    for(vector<CDV>::iterator it1 = cdv.begin();it1 != cdv.end(); ++it1){
      make_intra_bond(*it0,*it1,bond);
      make_hydrogen_bond(it0,it1,bond);
    }
  }
}
