#include "gro.h"

void read_gro_natom(char *buffer,int &natom){
  char c[10];
  for(int i = 0;i<5;++i){
    c[i] = buffer[i];
  }
  c[5] = '\0';
  natom = atoi(c);
};

static inline std::string trim(const std::string& string, const char* trimCharacterList = " \t\v\r\n")
{
  std::string result;
  std::string::size_type left = string.find_first_not_of(trimCharacterList);
  if (left != std::string::npos){
    std::string::size_type right = string.find_last_not_of(trimCharacterList);
    result = string.substr(left, right - left + 1);
  }
  return result;
}

void read_gro_atom(char *buffer,AtomGro *atoms,int n){

  char value[9];
  char str[6];
  int count = 0;
  //read ID of molecule
  for(int i = 0;i<5;++i){
    str[i] = buffer[count++];
  }
  str[5] = '\0';
  atoms[n].MolID = atoi(str);
  //read name of molecule
  for(int i = 0;i<5;++i){
    str[i] = buffer[count++];
  }
  str[5] = '\0';
  //atoms[n].MolName = std::string(c);
  //atoms[n].MolName = trim(atoms[n].MolName);
  sprintf(atoms[n].MolName,"%s",str);

  //read name of atom
  for(int i = 0;i<5;++i){
    str[i] = buffer[count++];
  }
  str[5] = '\0';
  //atoms[n].AtomName = std::string(c);
  sprintf(atoms[n].AtomName,"%s",str);

  //read id of atom
  for(int i = 0;i<5;++i){
    str[i] = buffer[count++];
  }
  str[5] = '\0';
  atoms[n].AtomID = atoi(str);
  //read coordinates of atom
  for(int d=0;d<3;d++){
    for(int i = 0;i<8;++i){
      value[i] = buffer[count++];
    }
    value[8] = '\0';
    atoms[n].coor[d] = atof(value);
  }
  //velocities of atom are set to zero;
  for(int d=0;d<3;d++) atoms[n].vel[d] = 0.0;
};

void read_gro_cellsize(char *buffer,dvec3 &cellsize){
  char c[11];
  int count = 0;
  for(int d = 0;d<3;++d){
    for(int i = 0;i<10;++i){
      c[i] = buffer[count++];
    }
    c[10] = '\0';
    cellsize[d] = atof(c);
  }
};

/* gro file format
0        1         2         3         4         5         6         7
1234567890123456789012345678901234567890123456789012345678901234567890
MD of 2 waters, t= 0.0
    6
    1WATER  OW1    1   0.126   1.624   1.679  0.1227 -0.0580  0.0434
    1WATER  HW2    2   0.190   1.661   1.747  0.8085  0.3191 -0.7791
    1WATER  HW3    3   0.177   1.568   1.613 -0.9045 -2.6469  1.3180
    2WATER  OW1    4   1.275   0.053   0.622  0.2519  0.3140 -0.1734
    2WATER  HW2    5   1.337   0.002   0.680 -1.0641 -1.1349  0.0257
    2WATER  HW3    6   1.326   0.120   0.568  1.9427 -0.8216 -0.0244
   1.82060   1.82060   1.82060
//*/

void read_gro_file(char *filename,int &natom,AtomGro *atom,dvec3 &cellsize,int &nmol){
  char buffer[256];

  FILE *fp;
  if((fp=fopen(filename,"r")) == NULL){
    fprintf(stderr,"error: fopen %s failed\n",filename);
    exit(EXIT_FAILURE);
  }
  fgets(buffer,70,fp);//skip comment
  fgets(buffer,70,fp);
  read_gro_natom(buffer,natom);
  for(int i = 0;i<natom;++i){
    fgets(buffer,70,fp);
    read_gro_atom(buffer,atom,i);
    std::cout << atom[i] << std::endl;
  }
  //read_gro_natom(buffer,nmol);
  fgets(buffer,70,fp);
  read_gro_cellsize(buffer,cellsize);

  fclose(fp);
};

void write_gro_file(const char *filename,const int natom,const AtomGro *atom,const dvec3 L,const char *comment){
  FILE *fp;
  if((fp=fopen(filename,"w")) == NULL){
    fprintf(stderr,"error: fopen %s failed\n",filename);
    exit(EXIT_FAILURE);
  }
  fprintf(fp,"%s\n",comment);
  fprintf(fp,"%5d\n",natom);
  for(int i=0;i<natom;i++){
    AtomGro a = atom[i];
    fprintf(fp,"%5d%s%s%5d%8.3lf%8.3lf%8.3lf%8.4lf%8.4lf%8.4lf\n",
	    a.MolID, a.MolName, a.AtomName, a.AtomID,
	    a.coor[0], a.coor[1], a.coor[2],
	    a.vel[0],  a.vel[1],  a.vel[2]);
  }
  fprintf(fp,"%10.5lf%10.5lf%10.5lf\n",L[0],L[1],L[2]);
  fclose(fp);
}
