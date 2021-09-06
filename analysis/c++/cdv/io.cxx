#include "io.h"

bool ks_strval(char *buf, char *key, double &v)
{
  int i;
  char c0[256];
  char *cp;

  cp = strstr(buf,key);
  if(cp == NULL){
    return false;
  }
  cp += strlen(key);
  //skip space space if there is no line break
  for(; *cp && *cp == ' ' && *cp != 0x0d && *cp != 0x0a;cp++);
  //read chars until next space or line break
  for(i = 0; *cp && *cp != ' ' && *cp != 0x0d && *cp != 0x0a;i++,cp++){
    c0[i] = *cp;
  }
  c0[i] = '\0';
  v = atof(c0);
  return true;
}

bool ks_strval3(char *buf, char *key, double &v0,double &v1,double &v2)
{
  int i;
  char c0[256];
  char *cp;

  cp = strstr(buf,key);
  if(cp == NULL){
    return false;
  }
  cp += strlen(key);
  //skip left bracket
  if(*cp == '(') cp++;else return false;
  //read chars until next comma or line break
  for(i = 0; *cp && *cp != ',' && *cp != 0x0d && *cp != 0x0a;i++,cp++){
    c0[i] = *cp;
  }
  c0[i] = '\0';
  v0 = atof(c0);
  //skip comma
  if(*cp == ',') cp++;else return false;
  //read chars until next comma or line break
  for(i = 0; *cp && *cp != ',' && *cp != 0x0d && *cp != 0x0a;i++,cp++){
    c0[i] = *cp;
  }
  c0[i] = '\0';
  v1 = atof(c0);
  //skip comma
  if(*cp == ',') cp++;else return false;
  //read chars until next right bracket or line break
  for(i = 0; *cp && *cp != ')' && *cp != 0x0d && *cp != 0x0a;i++,cp++){
    c0[i] = *cp;
  }
  c0[i] = '\0';
  v2 = atof(c0);
  return true;
}

void read_header(string line,Header &header){
  char c[256], key[256];
  sprintf(c,"%s",line.c_str());

  sprintf(key,"box_sx=");
  ks_strval(c,key,header.x.s);
  sprintf(key,"box_ex=");
  ks_strval(c,key,header.x.e);
  sprintf(key,"box_sy=");
  ks_strval(c,key,header.y.s);
  sprintf(key,"box_ey=");
  ks_strval(c,key,header.y.e);
  sprintf(key,"box_sz=");
  ks_strval(c,key,header.z.s);
  sprintf(key,"box_ez=");
  ks_strval(c,key,header.z.e);

  for(int i=0;i<MAX_ATOMTYPE;i++){
    sprintf(key,"r%d=",i);
    ks_strval(c,key, header.r[i]);
    sprintf(key,"c%d=",i);
    ks_strval3(c,key, header.c[i].r,header.c[i].g,header.c[i].b);
  }
  for(int i=0;i<MAX_BONDTYPE;i++){
    sprintf(key,"bond%d_c=",i);
    ks_strval3(c,key,header.b[i].r,header.b[i].g,header.b[i].b);
    sprintf(key,"bond%d_wt=",i);
    ks_strval(c,key,header.w[i]);
  }
}

void read_cdv(std::istream &is, Header &header, vector<CDV> &cdv, map< pair<int,int>, int> &bond){
  string dummy;
  cdv.clear();
  bond.clear();

  int count = 0;
  for(std::string line; std::getline(is, line);){
    std::istringstream ss(line);
    if(ss.fail()) continue;

    if(line[0]=='\'' || line[0]=='#'){
      read_header(line,header);
      continue;
    }
    if(isalpha(line[0])){
      int key1,key2,val;
      ss >> dummy >> key1 >> key2 >> val;
      pair<int,int> key(key1,key2);
      bond.insert( pair< pair<int,int> ,int>(key,val) );
      continue;
    }
    CDV tmp;
    ss >> tmp.m >> tmp.t >> tmp.r[0] >> tmp.r[1] >> tmp.r[2];
    if(tmp.t > 2)  continue;
    tmp.i = count++;
    if(tmp.t == 2) tmp.t = 1;
    cdv.push_back(tmp);
  }
}

void write_cdv(std::ostream &os, Header header,vector<CDV> cdv, map< pair<int,int>, int> bond){
  header.Print(os);
  for(auto it = cdv.begin();it != cdv.end(); ++it){
    it->Print(os);
  }
  for(auto it = bond.begin();it != bond.end(); ++it){
    pair<int,int> index = it->first;
    int type = it->second;
    os << "CDVIEW_BOND " << index.first << " " << index.second << " " << type << " " << endl;
  }
}
