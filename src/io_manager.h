#ifndef H_INPUT
#define H_INPUT

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <assert.h>

#include <iostream>
#include <fstream>

#include <molfile_plugin.h>
#include "vector.h"
#include "gro.h"
#include "molecules.h"

#include "parameter.h"

class InputParam{
 private:
 public:
  int  nmol;
  int  natom;
  char igro[256];  // input gro file name
  char itrr[256];  // input trr file name
  char ichk[256];  // input check point file
  char itkw[256];  // input takaiwa file
  char ixyz[256];  // input xyz file
  char icdv[256];  // input cdv file

  char ogro[256];  // output gro file name
  char otrr[256];  // output trr file name
  char ochk[256];

  char oene[256];  // output energy file name
  char ocdv[256];

  InputParam(int,char**);
  InputParam(const Parameter);
};

class IOManager{
 private:
  //void *rg,*wg;
  FILE *rg;
  char buffer[256];

  void *rt,*wt; //trr read and write pointer

  int      natom;   // number of atoms
  AtomGro  *atom;   // information of atoms in .gro format
  int      nmol;    // number of molecules
  //dvec3    L;       // system size
  //molfile_timestep_t timestep;

  MolType  SearchMolType(const std::string,MolTypeList);
  //int      ReadInputFile();

  dvec4    CalcAngleFromAtoms(const dvec3,AtomGro*);
  dvec4    CalcAngVelFromAtoms(const Molecule,AtomGro*,const MolType);

  void     ReadGroInput(dvec3&);
  void     ReadTrrInput();
  void     ReadCDVInput(dvec3&);
  void     ReadXYZInput(dvec3&);
  void     ReadCheckPoint(Molecule*,Thermostat*,Barostat*,dvec3&,const MolTypeList,const int,const int,Property&);
  void     ReadTakaiwaFormat(Molecule*,dvec3&,const int);

  void     UpdateTrrOutput(const dvec3);
  void     UpdateGroOutput(const dvec3);

  Molecule ConvertAtomsToMolecule(AtomGro*, const MolType,const dvec3);
  void     ConvertAtomsToMolecules(Molecule*,MolTypeList,const dvec3);
  void     ConvertMoleculesToAtoms(const Molecule*,MolTypeList,const dvec3);

 public:
  InputParam *iprm;   // input parameter

  IOManager(char *filename,dvec3 &L,int nmol);
  IOManager(const Parameter);
  IOManager(int,char**);
  ~IOManager();

  void     OpenFiles();
  void     ReadInputs(Molecule*,Thermostat*,Barostat*,MolTypeList,const int,const int,dvec3&,Property&);
  void     ReadInputs(Molecules* m){ReadInputs(m->mlcl,m->tst,m->bst,m->mlist,m->nmol,m->natom,m->L,m->prop);};

  void     UpdateOutputs(const Molecule*,const Thermostat*,const Barostat*,const Property,const dvec3,MolTypeList);
  void     UpdateOutputs(const Molecules* m){UpdateOutputs(m->mlcl,m->tst,m->bst,m->prop,m->L,m->mlist);};
  void     WriteCDV(const Molecule*,const MolTypeList,const dvec3,const long);
  void     WriteCDV(const Molecules* m,const long s){WriteCDV(m->mlcl,m->mlist,m->L,s);};
  void     WriteCheckPoint(const Molecule*,const Thermostat*,const Barostat*,const dvec3,const Property);

  int      GetNumberOfMolecules(){return nmol;};
  //dvec3    GetSytemSize(){return L;};

  void     PrintOptions();
};
#endif

