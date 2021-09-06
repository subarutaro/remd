#include "vector.h"
#include "molecules.h"
#include "io_manager.h"
#include "parameter.h"
#include "gro.h"

#ifndef GROMACSPLUGIN
#define GROMACSPLUGIN
#include <gromacsplugin.C>
#endif

#define NMOLMAX  2048
#define NATOMMAX (5*NMOLMAX)

static
void print_help(char **argv){
	fprintf(stderr,"%s is a MD simulation program of rigid-body rotational molecules\n",argv[0]);
	fprintf(stderr,"Options  Default     Type        Description\n");
	fprintf(stderr,"------------------------------------------------------------------------------\n");
	fprintf(stderr,"-nm      0           Input       number of molecules\n");
	fprintf(stderr,"-na      0           Input       number of atoms\n");
	fprintf(stderr,"-c       input.gro   Input       Input gro file (use *.gro), you can input \"fcc\" to use face center cubic structure to this option\n");
	fprintf(stderr,"-t       input.trr   Input Opt.  Input trr file (use *.trr)\n");
	fprintf(stderr,"-oc      output.gro  Output      Output gro file\n");
	fprintf(stderr,"-ot      output.trr  Output Opt. Output trr file\n");
	fprintf(stderr,"-oe      output.ene  Output Opt. Output energy file\n");
};

InputParam::InputParam(int argc,char **argv){
  //input default value
  nmol = 0;
  natom = 0;
  sprintf(igro, "input.gro");
  sprintf(itrr, "null");
  sprintf(ogro, "output.gro");
  sprintf(otrr, "output.trr");
  sprintf(oene, "output.ene");
  sprintf(ocdv, "output");
  //read options from input
  for(int i=1;i<argc;i++){
    if(*argv[i]=='-'){
      switch(*(argv[i]+1)){
      case 'n':
	switch(*(argv[i]+2)){
	case 'm':
	  nmol = atoi(argv[++i]);break;
	case 'a':
	  natom = atoi(argv[++i]);break;
	}
      case 'c':
	sprintf(igro,"%s",argv[++i]);break;
      case 't':
	sprintf(itrr,"%s",argv[++i]);break;
      case 'h':
	print_help(argv);exit(0);break;
      case 'o':
	printf("output file input\n");
	switch(*(argv[i]+2)){
	case 'c':
	  sprintf(ogro,"%s",argv[++i]);break;
	case 't':
	  sprintf(otrr,"%s",argv[++i]);break;
	case 'e':
	  sprintf(oene,"%s",argv[++i]);break;
	default:
	  fprintf(stderr,"error: undifined option %s\n",argv[i]);
	  exit(EXIT_FAILURE);
	  break;
	}
	break;
      default:
	fprintf(stderr,"error: undefined option \"%s\"\n",argv[i]);
	exit(EXIT_FAILURE);
	break;
      }
    }else{
      fprintf(stderr,"error: undifined option! %s\n",argv[i]);
      exit(EXIT_FAILURE);
    }
  }
  if(nmol == 0){
    fprintf(stderr,"error: nmol must be defined with \"-n\"!\n");
    exit(EXIT_FAILURE);
  }
}

InputParam::InputParam(const Parameter param){
  //input default value
  nmol = param.nmol;
  natom = param.natom;

  if(param.gro_in!="null"){
    sprintf(igro,"%s",(param.input_prefix+param.gro_in+".gro").c_str());
  }else{
    sprintf(igro,"%s",param.gro_in.c_str());
  }
  if(param.trr_in!="null"){
    sprintf(itrr,"%s",(param.input_prefix+param.trr_in+".trr").c_str());
  }else{
    sprintf(itrr,"%s",param.trr_in.c_str());
  }
  if(param.chk_in!="null"){
    sprintf(ichk,"%s",(param.input_prefix+param.chk_in+".chk").c_str());
  }else{
    sprintf(ichk,"%s",param.chk_in.c_str());
  }

  if(param.tkw_in!="null"){
    sprintf(itkw,"%s",(param.input_prefix+param.tkw_in+".tkw").c_str());
  }else{
    sprintf(itkw,"%s",param.tkw_in.c_str());
  }

  if(param.xyz_in!="null"){
    sprintf(ixyz,"%s",(param.input_prefix+param.xyz_in+".xyz").c_str());
  }else{
    sprintf(ixyz,"%s",param.xyz_in.c_str());
  }

  if(param.cdv_in!="null"){
    sprintf(icdv,"%s",(param.input_prefix+param.cdv_in+".cdv").c_str());
  }else{
    sprintf(icdv,"%s",param.cdv_in.c_str());
  }

  if(param.gro_out!="null"){
    sprintf(ogro,"%s",(param.output_prefix+param.gro_out+".gro").c_str());
  }else{
    sprintf(ogro,"%s",param.gro_out.c_str());
  }
  if(param.trr_out!="null"){
    sprintf(otrr,"%s",(param.output_prefix+param.trr_out+".trr").c_str());
  }else{
    sprintf(otrr,"%s",param.trr_out.c_str());
  }
  if(param.chk_out!="null"){
    sprintf(ochk,"%s",(param.output_prefix+param.chk_out+".chk").c_str());
  }else{
    sprintf(ochk,"%s",param.chk_out.c_str());
  }
  if(param.ene_out!="null"){
    sprintf(oene,"%s",(param.output_prefix+param.ene_out+".ene").c_str());
  }else{
    sprintf(oene,"%s",param.ene_out.c_str());
  }
  if(param.cdv_out!="null"){
    sprintf(ocdv,"%s",(param.output_prefix+param.cdv_out).c_str());
  }else{
    sprintf(ocdv,"%s",param.cdv_out.c_str());
  }
}

static inline
dvec4 euler_to_quartanion(const dvec3 e){
  dvec4 q;
  q[0] = cos(0.5*e[1])*cos(0.5*(e[0]+e[2]));
  q[1] = sin(0.5*e[1])*cos(0.5*(e[0]-e[2]));
  q[2] = sin(0.5*e[1])*sin(0.5*(e[0]-e[2]));
  q[3] = cos(0.5*e[1])*sin(0.5*(e[0]+e[2]));
  return q;
}

/* --------------------------------------------
// only for water
// all atom OW, HW1 and HW2 must be on xy plane
// vector GO must be on x axis of body space
// vector GH1.y must > 0
// vector GH2.y must < 0
//   y
//   ^
//  H|
//   \
// --|O--> x
//   /
//  H|
///////////////////////////////
 ------------------------------------------- */

IOManager::IOManager(const Parameter param){
  iprm = new InputParam(param);
  OpenFiles();
}

IOManager::IOManager(int argc,char **argv){
  iprm = new InputParam(argc,argv);
  OpenFiles();
};

IOManager::~IOManager(){
  fclose(rg);
  //close_gro_read(rg);
  //close_gro_write(wg);
  if(strcmp(iprm->itrr,"null")){
    close_trr_read(rt);
  }
  close_trr_write(wt);
};

void IOManager::OpenFiles(){
  if(strcmp(iprm->ichk,"null")){
    FILE *fchk;
    if((fchk=fopen(iprm->ichk,"r")) == NULL){
      fprintf(stderr,"error: fopen %s failed\n",iprm->ichk);
      exit(EXIT_FAILURE);
    }
    fscanf(fchk,"%d %d",&nmol,&natom);
    fclose(fchk);
  }else if(strcmp(iprm->igro,"null")){
    if((rg=fopen(iprm->igro,"r")) == NULL){
      fprintf(stderr,"error: fopen %s failed\n",iprm->igro);
      exit(EXIT_FAILURE);
    }
    fgets(buffer,70,rg);//skip comment
    fgets(buffer,70,rg);
    read_gro_natom(buffer,natom);

    if(strcmp(iprm->itrr,"null")){
      rt = open_trr_read(iprm->itrr,"trr",&natom);
      if (!rt) {
	fprintf(stderr, "%s: line %d\n", __FILE__,__LINE__);
	exit(EXIT_FAILURE);
      }
      std::cout << "open_trr_read scceeded" << std::endl;
    }
  }

  if(strcmp(iprm->otrr,"null")){
    wt = open_trr_write(iprm->otrr,"trr",natom);
    if (!wt) {
      fprintf(stderr, "%s: line %d\n", __FILE__,__LINE__);
      exit(EXIT_FAILURE);
    }
  }
}

void IOManager::ReadGroInput(dvec3 &L){
  const double nm_to_m = 1e-9;
  for(int i = 0;i<natom;++i){
    fgets(buffer,70,rg);
    read_gro_atom(buffer,atom,i);
  }
  nmol = atom[natom-1].MolID;

  fgets(buffer,70,rg);
  read_gro_cellsize(buffer,L);

  for(int i=0;i<natom;i++){
    atom[i].coor *= nm_to_m / unit_length;
    atom[i].vel  *= unit_time / 1e3*unit_length;
  }
  L *= nm_to_m / unit_length;
};

//*
//timestep with class AtomGro. the file format has a header with each record
static int read_trr_timestep(void *v, int natoms, AtomGro *atom, dvec3 &L) {
  gmxdata *gmx = (gmxdata *)v;
  md_ts mdts;
  memset(&mdts, 0, sizeof(md_ts));
  mdts.natoms = natoms;

  if (mdio_timestep(gmx->mf, &mdts) < 0) {
    if (mdio_errno() == MDIO_EOF || mdio_errno() == MDIO_IOERROR) {
      // XXX Lame, why does mdio treat IOERROR like EOF?  
      return MOLFILE_ERROR;
    }
    fprintf(stderr, "gromacsplugin) Error reading timestep, %s\n", 
            mdio_errmsg(mdio_errno()));
    return MOLFILE_ERROR;
  }
  if (mdts.natoms != natoms) {
    fprintf(stderr, "gromacsplugin) Timestep in file contains wrong number of atoms\n");
    fprintf(stderr, "gromacsplugin) Found %d, expected %d\n", mdts.natoms, natoms);
    mdio_tsfree(&mdts);
    return MOLFILE_ERROR;
  }

  if (atom) {
    for(int i=0;i<natoms;i++){
      for(int d=0;d<3;d++){
	std::cout << mdts.pos[3*i+d] << " ";
	atom[i].coor[d] = (double)mdts.pos[3*i+d];
      }
      std::cout << std::endl;
    }
    if (mdts.box) {
      L[0] = mdts.box->A;
      L[1] = mdts.box->B;
      L[2] = mdts.box->C;
      /*
      ts->alpha = mdts.box->alpha;
      ts->beta = mdts.box->beta;
      ts->gamma = mdts.box->gamma;
      //*/
    }
  }
  mdio_tsfree(&mdts);
  return MOLFILE_SUCCESS;
}
//*/

void IOManager::ReadTrrInput(){
  dvec3 dummyL;
  read_trr_timestep(rt,natom,atom,dummyL);
};

void IOManager::ReadTakaiwaFormat
(Molecule   *m,
 dvec3      &L,
 const int  _nmol
){
  std::cout << "ReadTakaiwaFormat" << std::endl;
  std::istream *is;
  is = new std::ifstream(iprm->itkw,std::ios_base::out);

  std::string dummy;
  *is >> L[0];
  L[1] = L[2] = L[0];
  std::cout << L[0] << std::endl;

  for(int i=0;i<nmol;i++){
    *is >> dummy;
    *is >> m[i].r[0] >> m[i].r[1] >> m[i].r[2];
    *is >> m[i].q[0] >> m[i].q[1] >> m[i].q[2] >> m[i].q[3];
    std::cout << i;
    std::cout << " " << m[i].r[0] << " " << m[i].r[1] << " " << m[i].r[2];
    std::cout << " " << m[i].q[0] << " " << m[i].q[1] << " " << m[i].q[2] << " " << m[i].q[3];
    std::cout << std::endl;
  }

  delete is;

  //*
  for(int i=0;i<nmol;i++){
    m[i].r /= L;
  }
}

void IOManager::ReadXYZInput(dvec3 &L){
  std::cout << "ReadXYZInput" << std::endl;
  std::cout << "opening " << iprm->ixyz << std::endl;
  std::istream *is;
  is = new std::ifstream(iprm->ixyz,std::ios_base::out);

  int na;
  std::string dummy;
  *is >> na;
  *is >> L[0] >> dummy;
  L[1] = L[2] = L[0];
  std::cout << na << std::endl;
  std::cout << L[0] << std::endl;

  int count=0;
  for(int i=0;i<nmol;i++){
    for(int j=0;j<na/nmol;j++){
      *is >> dummy;
      std::cout << " " << dummy;
      for(int d=0;d<3;d++){
	*is >> atom[count].coor[d];
	atom[count].vel[d] = 0.0;
	std::cout << " " << atom[count].coor[d];
      }
      std::cout << std::endl;
      count++;
    }
  }

  delete is;
}

void IOManager::ReadCDVInput(dvec3 &L){
  std::cout << "Reading CDV file " << iprm->icdv << std::endl;
  std::string line;
  std::ifstream is(iprm->icdv);

  std::getline(is,line);
  size_t beg,end;
  beg = line.find(" box_sz=") + 8;
  end = line.find(" box_ex=");
  std::string sz = line.substr(beg,end-beg);
  //std::cout << beg << ' ' << end << ' ' << sz << std::endl;
  L[2] = -2.0 * std::stod(sz);
  for(int count=0;count < natom;count++){
    std::getline(is,line);
    std::stringstream strs(line);
    int dummy;
    double x,y,z;
    strs >> dummy >> dummy >> x >> y >> z;
    atom[count].coor[0] = x + 0.5 * L[0];
    atom[count].coor[1] = y + 0.5 * L[1];
    atom[count].coor[2] = z + 0.5 * L[2];
  }
}


void IOManager::ReadInputs(Molecule *mlcl,Thermostat *tst,Barostat *bst,MolTypeList mtype,const int _nmol,const int _natom,dvec3 &L,Property &prop){
  if(strcmp(iprm->ichk,"null")){
    std::cout << iprm->ichk << std::endl;
    ReadCheckPoint(mlcl,tst,bst,L,mtype,_nmol,_natom,prop);
    if(strcmp(iprm->itkw,"null"))
      ReadTakaiwaFormat(mlcl,L,_nmol);
  }else if(strcmp(iprm->igro,"null")){
    std::cout << iprm->igro << std::endl;
    assert(natom==_natom);
    atom = (AtomGro*)malloc(natom*sizeof(AtomGro));
    ReadGroInput(L);
    assert(nmol==_nmol);
    if(strcmp(iprm->itrr,"null")) ReadTrrInput();
    if(strcmp(iprm->ixyz,"null")) ReadXYZInput(L);
    if(strcmp(iprm->icdv,"null")) ReadCDVInput(L);
    ConvertAtomsToMolecules(mlcl,mtype,L);
    free(atom);
  }else{
    fprintf(stderr,"error: no input .gro or .chk file\n");
    exit(EXIT_FAILURE);
  }
};

void IOManager::UpdateGroOutput(const dvec3 L){
  char comment[256];
  sprintf(comment,"test");
  write_gro_file(iprm->ogro,natom,atom,L,comment);
};

static int write_trr_timestep(void *mydata, const AtomGro *ts, const dvec3 L,const dvec3 angle=90.0) 
{
  const float nm=0.1;
  gmxdata *gmx = (gmxdata *)mydata;
  // determine and write header from structure info.
  // write trr header. XXX: move this to Gromacs.h ??
  if (gmx->mf->fmt == MDFMT_TRR) {
    int i,d;
    if ( put_trx_int(gmx->mf, TRX_MAGIC)            // ID
         || put_trx_string(gmx->mf, "GMX_trn_file") // version
         || put_trx_int(gmx->mf, 0)                 // ir_size (ignored)
         || put_trx_int(gmx->mf, 0)                 // e_size (ignored)
         || put_trx_int(gmx->mf, 9*sizeof(float))   // box
         || put_trx_int(gmx->mf, 0)                 // vir_size (ignored)
         || put_trx_int(gmx->mf, 0)                 // pres_size (ignored)
         || put_trx_int(gmx->mf, 0)                 // top_size (ignored)
         || put_trx_int(gmx->mf, 0)                 // sym_size (ignored)
         || put_trx_int(gmx->mf, 3*sizeof(float)*gmx->natoms) // coordinates 
         || put_trx_int(gmx->mf, 0)                 // no velocities
         || put_trx_int(gmx->mf, 0)                 // no forces
         || put_trx_int(gmx->mf, gmx->natoms)       // number of atoms
         || put_trx_int(gmx->mf, gmx->step)         // current step number
         || put_trx_int(gmx->mf, 0)                 // nre (ignored)
         || put_trx_real(gmx->mf, 0.1*gmx->step)    // current time. (dummy value: 0.1)
         || put_trx_real(gmx->mf, 0.0))             // current lambda
      return MOLFILE_ERROR;

    // set up box according to the VMD unitcell conventions.
    // the a-vector is collinear with the x-axis and
    // the b-vector is in the xy-plane. 
    const float sa = sin((double)angle[0]/180.0*M_PI);
    const float ca = cos((double)angle[0]/180.0*M_PI);
    const float cb = cos((double)angle[1]/180.0*M_PI);
    const float cg = cos((double)angle[2]/180.0*M_PI);
    const float sg = sin((double)angle[2]/180.0*M_PI);
    float box[9];
    box[0] = L[0];    box[1] = 0.0;     box[2] = 0.0;
    box[3] = L[1]*ca; box[4] = L[1]*sa; box[5] = 0.0;
    box[6] = L[2]*cb; box[7] = L[2]*(ca - cb*cg)/sg;
    box[8] = L[2]*sqrt((double)(1.0 + 2.0*ca*cb*cg 
                                 - ca*ca - cb*cb - cg*cg)/(1.0 - cg*cg));

    for (i=0; i<9; ++i) {
      if (put_trx_real(gmx->mf, box[i]*nm))
        return MOLFILE_ERROR;
    }

    // write coordinates
    for (i=0; i<gmx->natoms; ++i){
      for (d=0; d<3; ++d){
	if (put_trx_real(gmx->mf, (float)ts[i].coor[d]*nm))
	  return MOLFILE_ERROR;
      }
    }
  } else {
    fprintf(stderr, "gromacsplugin) only .trr is supported for writing\n");
    return MOLFILE_ERROR;
  }

  ++ gmx->step;
  return MOLFILE_SUCCESS;
}

void IOManager::UpdateTrrOutput(const dvec3 L){
  write_trr_timestep(wt,atom,L);
};

void IOManager::WriteCDV
(const Molecule *mlcl,const MolTypeList mtype,const dvec3 L,const long step)
{
  if(strcmp(iprm->ocdv,"null")){
    const double A_to_nm = 0.1;
    atom = (AtomGro*)malloc(natom*sizeof(AtomGro));
    ConvertMoleculesToAtoms(mlcl,mtype,L);

    std::stringstream strs;
    strs << iprm->ocdv << "_" << std::setw(8) << std::setfill('0') << step << ".cdv";
    std::string filename = strs.str();

    std::ostream *s;
    //s = new std::iostream(std::cout.rdbuf());
    s = new std::ofstream(filename.c_str(),std::ios_base::out);

    //output comments
    *s << "' box_sx=0.0 box_sy=0.0 box_sz=0.0"; // start of box for all axises
    *s << " box_ex=" << L[0]*A_to_nm; //end of box for x axis
    *s << " box_ey=" << L[1]*A_to_nm; //end of box for y axis
    *s << " box_ez=" << L[2]*A_to_nm; //end of box for z axis
    *s << std::endl;
    //output
    for(int i=0;i<natom;i++){
      *s << std::setw(10) << i << " " << atom[i].AtomID%3 << " ";
      *s << atom[i].coor << std::endl;
    }
    delete s;
    free(atom);
  }
}

void IOManager::ReadCheckPoint
(Molecule   *m,
 Thermostat *tst,
 Barostat   *bst,
 dvec3      &L,
 const MolTypeList ml,
 const int  _nmol,
 const int  _natom,
 Property   &prop
)
{
  std::istream *is;
  is = new std::ifstream(iprm->ichk,std::ios_base::out);

  *is >> nmol >> natom;
  //std::cout << nmol << " " << natom << std::endl;
  //std::cout << _nmol << " " << _natom << std::endl;

  //assert(nmol==_nmol && natom==_natom);
  if(nmol != _nmol || natom != _natom){
    fprintf(stderr,"warning: input # of moleucles(=%d) and atoms(=%d) are not equal to those of input file(=%d, %d)",
	    nmol,natom,_nmol,_natom);
    fprintf(stderr,"warning: nmol and natom is fixed as value in input file\n");
    nmol = _nmol;
    natom = _natom;
  }

  int atom_count = 0;
  for(int i=0;i<nmol;i++){
    *is >> m[i].r[0] >> m[i].r[1] >> m[i].r[2];
    *is >> m[i].v[0] >> m[i].v[1] >> m[i].v[2];
    *is >> m[i].q[0] >> m[i].q[1] >> m[i].q[2] >> m[i].q[3];
    *is >> m[i].p[0] >> m[i].p[1] >> m[i].p[2] >> m[i].p[3];
    *is >> m[i].m;
    *is >> m[i].i[0] >> m[i].i[1] >> m[i].i[2];
    *is >> m[i].type;
    *is >> m[i].id;
    if(m[i].id != atom_count){
      fprintf(stderr,"warning: input atom id(=%d) is not equal to that of on time counting=(%d).\n",m[i].id,atom_count);
      m[i].id = atom_count;
    }
    atom_count += ml[m[i].type].a.size();
  }
  *is >> tst->s >> tst->Ps >> tst->Q;;
  *is >> bst->Pv[0] >> bst->Pv[1] >> bst->Pv[2] >> bst->W;
  *is >> L[0] >> L[1] >> L[2];

  delete is;

  for(int i=0;i<nmol;i++){
    m[i].r /= L;
    m[i].v *= L*tst->s;
    m[i].p *= tst->s;
  }
  /*
  for(int i=0;i<nmol;i++){
    std::cout << m[i].r;
    std::cout << m[i].v;
    std::cout << m[i].q;
    std::cout << m[i].p;
    std::cout << std::endl;
  }
  //*/
}

void IOManager::WriteCheckPoint
(const Molecule   *m,
 const Thermostat *tst,
 const Barostat   *bst,
 const dvec3      L,
 const Property   prop
)
{
  //std::ostream *is;
  //s = new std::iostream(std::cout.rdbuf());
  //s = new std::ofstream(filename.c_str(),std::ios_base::out);

  std::ostream *ofs;
  ofs = new std::ofstream(iprm->ochk,std::ios_base::out);

  *ofs << nmol << " " << natom << std::endl;
  for(int i=0;i<nmol;i++){
    if(m[i].r[0]>1.0 || m[i].r[0]<0.0 || m[i].r[1]>1.0 || m[i].r[1]<0.0 || m[i].r[2]>1.0 || m[i].r[2]<0.0){
      fprintf(stderr, "error: %d th molecule is out of box\n",i);
      exit(EXIT_FAILURE);
    }
    *ofs << m[i].r * L;
    *ofs << m[i].v / (L * tst->s);
    *ofs << m[i].q;
    *ofs << m[i].p / tst->s;
    *ofs << " " << m[i].m;
    *ofs << " " << m[i].i;
    *ofs << " " << m[i].type;
    *ofs << " " << m[i].id;
    *ofs << std::endl;
  }
  *ofs << tst->s << " " << tst->Ps << " " << tst->Q << std::endl;
  *ofs << bst->Pv << " " << bst->W << std::endl;
  *ofs << L << std::endl;
  *ofs << prop << std::endl;

  delete ofs;
}

void IOManager::UpdateOutputs
(
 const Molecule   *mlcl,
 const Thermostat *tst,
 const Barostat   *bst,
 const Property   prop,
 const dvec3      L,
 MolTypeList      mtype
 ){
  const double A_to_nm = 0.1;
  atom = (AtomGro*)malloc(natom*sizeof(AtomGro));
  ConvertMoleculesToAtoms(mlcl,mtype,L);
  if(strcmp(iprm->ogro,"null"))
    UpdateGroOutput(L*A_to_nm);
  if(strcmp(iprm->otrr,"null"))
    UpdateTrrOutput(L*A_to_nm);
  if(strcmp(iprm->ochk,"null"))
    WriteCheckPoint(mlcl,tst,bst,L,prop);

  free(atom);
};

dvec4 IOManager::CalcAngleFromAtoms(const dvec3 R, AtomGro *atom){
  dvec3 rx = atom[0].coor - R;
  //dvec3 ry = r[1] - R;
  dvec3 rz = atom[2].coor - R;

  dvec3 x = rx / norm(rx);
  dvec3 z = vector_prod(rz,rx);
  z /= norm(z);
  dvec3 y = vector_prod(z,x);

  dvec3 euler;
  const double e = 1e-10;

  dvec3 xaxis,yaxis,zaxis;
  xaxis[0]=1.0;xaxis[1]=0.0;xaxis[2]=0.0;
  yaxis[0]=0.0;yaxis[1]=1.0;yaxis[2]=0.0;
  zaxis[0]=0.0;zaxis[1]=0.0;zaxis[2]=1.0;

  dvec3 test;
  for(int i=0;i<4;i++){
  for(int j=0;j<3;j++){
  for(int k=0;k<4;k++){
	test[0] =i*0.5*M_PI;
	test[1] =j*0.5*M_PI;
	test[2] =k*0.5*M_PI;
	dvec4 tmp = euler_to_quartanion(test);
	if(norm(x-space_to_body(tmp,xaxis))<e){// operator< is not suitable for this code
	  if(norm(y-space_to_body(tmp,yaxis)),e){
	    if(norm(z-space_to_body(tmp,zaxis)),e){
	      euler = test;
	      return euler_to_quartanion(euler);
	}}}
  }}}

  euler[0] = atan2(z[0],-z[1]); // 0 < phi < 2pi
  euler[1] = acos(z[2]);
  euler[2] = atan2(x[2],y[2]);  // 0 < phi < 2pi
  return euler_to_quartanion(euler);
};

dvec4 IOManager::CalcAngVelFromAtoms
(const Molecule m,AtomGro *atoms,const MolType mt)
{
  dvec3 w = 0.0;
  for(unsigned int i=0;i<mt.a.size();i++){
    w += vector_prod(mt.a[i].r, space_to_body(m.q,atoms[i].vel))*mt.a[i].m;
  }
  w /= m.i;

  dvec4 p,w4;
  w4[0] = 0.0;
  w4[1] = w[0];
  w4[2] = w[1];
  w4[3] = w[2];
  p = scalar_prod(S(m.q),w4)*2.0;
  p[1] *= m.i[0];
  p[2] *= m.i[1];
  p[3] *= m.i[2];
};
//*/
MolType IOManager::SearchMolType
(const std::string mname,MolTypeList mtype)
{
  for(MolTypeList::iterator it = mtype.begin(); it != mtype.end(); it++){
    if(it->name == mname){
      //std::cout << it->name << std::endl;
      return (*it);
    }
  }
  fprintf(stderr,"error: type name \"%s\" is not available!\n",mname.c_str());
  exit(EXIT_FAILURE);
};

Molecule IOManager::ConvertAtomsToMolecule
(AtomGro *a, const MolType mt,const dvec3 L)
{
  const int n = mt.a.size();
  Molecule m;
  for(int i=0;i<n;i++){
    if(!strcmp(a[i].AtomName,mt.a[i].name.c_str())){
      if(i==0) m.id = a[i].AtomID - 1; // index starts from 1 in *.gro
      double mass = mt.a[i].m;
      m.r += a[i].coor * mass;
      m.v += a[i].vel  * mass;
      m.m += mass;
    }else{
      fprintf(stderr,"error: %s read from file is not correspond to %s\n",a[i].AtomName,mt.a[i].name.c_str());
      exit(EXIT_FAILURE);
    }
  }
  m.r /= m.m;
  m.v /= m.m;

  m.q = CalcAngleFromAtoms(m.r, a);
  //m.p = CalcAngVelFromAtoms(m, a, mt);

  m.i = mt.i;

  // convert to normalized coord and vel
  m.r /= L;
  m.v *= L;

  //std::cout << *a << std::endl;
  //std::cout << m << std::endl;
  return m;
};

void IOManager::ConvertAtomsToMolecules
(Molecule *mlcl,MolTypeList mtype,const dvec3 L)
{
  AtomGro *a = atom;
  int mid = 0;
  while(mid < nmol){
    MolType mt = SearchMolType(a->MolName,mtype);
    for(unsigned int i=0;i<mtype.size();i++)
      if(mt.name==mtype[i].name)
	 mlcl[mid].type = i;
    mlcl[mid++] = ConvertAtomsToMolecule(a, mt, L);
    a += mt.a.size();
    //std::cout << "r= " << mlcl[mid-1].r << std::endl;
    //std::cout << "q= " << mlcl[mid-1].q << " " << norm(mlcl[mid-1].q) << std::endl;
  }
}

static dvec3 CalcVelFromMolecule(const Molecule m,const AtomType at){
  dvec3 v_tra = m.v;
  dvec4 w4 = scalar_prod(S_T(m.q),m.p)*0.5;
  dvec3 w;
  w[0] = w4[1]*m.i[0];
  w[1] = w4[2]*m.i[1];
  w[2] = w4[3]*m.i[2];
  double r = sum(at.r * at.r);
  dvec3 v_rot = body_to_space(m.q,w*r);

  return v_tra+v_rot;
}

void IOManager::ConvertMoleculesToAtoms
(const Molecule *mlcl,MolTypeList mtype,const dvec3 L)
{
  const double A_to_nm = unit_length / 1e-9;
  int count = 0;
  for(int i=0;i<nmol;i++){
    Molecule m = mlcl[i];
    MolType mt = mtype[mlcl[i].type];
    for(unsigned int j=0;j<mt.a.size();j++){
      sprintf(atom[count].MolName,"%s",mt.name.c_str());
      atom[count].MolID = i+1;
      sprintf(atom[count].AtomName,"%s",mt.a[j].name.c_str());
      atom[count].AtomID = count+1;
      atom[count].coor = (m.r*L + body_to_space(mlcl[i].q,mt.a[j].r))*A_to_nm;// convert angstrom to nm
      //atom[count].vel = CalcVelFromMolecule(m,mt.a[j]);
      atom[count].vel = 0.f;
      count++;
    }
  }
}

void IOManager::PrintOptions()
{
  std::cout << "-c  " << iprm->igro  << std::endl;
  std::cout << "-t  " << iprm->itrr  << std::endl;
  std::cout << "-oc " << iprm->ogro  << std::endl;
  std::cout << "-ot " << iprm->otrr  << std::endl;
};

IOManager::IOManager
(char *filename,dvec3 &L,int nmol)
{
  //int nmol_tmp;
  //atom = (AtomGro*)malloc(5*nmol*isizeof(AtomGro)); //natom in a molecule must be < 6

  //ReadInputFile(filename,natom,atom,L,nmol_tmp);
  //assert(nmol != nmol_tmp);
  //assert(natom > 5*nmol);
};
/*
int IOManager::ReadInputFile(){
  printf("---reading %s---\n",iprm->igro);
  read_gro_file(iprm->igro,natom,atom,L,nmol);
  if(strcmp(iprm->itrr,"null")){
    printf("---reading %s---\n",iprm->itrr);
    timestep = read_trr_file(iprm->itrr,natom);
    for(int i=0;i<natom;i++){
      for(int d=0;d<3;d++){
	atom[i].coor[d] = timestep.coords[3*i+d] * 1e-1;
      }
      std::cout << atom[i] << std::endl;
    }
    printf("---writing %s---\n",iprm->otrr);
    write_trr_file(iprm->otrr,timestep,natom);
  }
  printf("---writing %s---\n",iprm->ogro);
  write_gro_file(iprm->ogro,natom,atom,L);
  std::cout << "end read input file " << std::endl;

  return 0;
}
//*/
