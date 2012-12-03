#include <fstream.h>
#include "phonlib.h"
#include "getvalue.h"
#include "parse.h"
#include "version.h"
#include "plugin.h"

void read_spring_file(LinkedList<SpringSvsL> *springlist, int *pmaxlpow, int *pdirdep, int *pconcpow, istream &file, const Array<AutoString> &label) {
  *pmaxlpow=-1;
  *pdirdep=0;
  *pconcpow=0;
  while (skip_delim(file)) {
    AutoString atom;
    SpringSvsL s;
    get_string(&atom,file);
    if (strcmp(atom,"maxlpow")==0) {
      file >> *pmaxlpow;
    }
    else if (strcmp(atom,"dirdep")==0) {
      file >> *pdirdep;
    }
    else if (strcmp(atom,"concpow")==0) {
      file >> *pconcpow;
    }
    else {
      s.atom_type[0]=index_in_array(label,atom);
      skip_delim(file);
      get_string(&atom,file);
      s.atom_type[1]=index_in_array(label,atom);
      s.forcek.resize(2);
      file >> s.forcek(0);
      file >> s.forcek(1);
      (*springlist) << new SpringSvsL(s);
    }
  }
  if (*pmaxlpow==-1) {
    LinkedListIterator<SpringSvsL> i(*springlist);
    *pmaxlpow=i->forcek(0).get_size()-1;
  }
}

void read_spring_file(LinkedList<ExtendableSpringSvsL> *springlist, istream &file, const Array<AutoString> &label, char *potplugin) {
  while (skip_delim(file)) {
    AutoString atom;
    ExtendableSpringSvsL *pspring=new ExtendableSpringSvsL;
    get_string(&atom,file);
    pspring->atom_type[0]=index_in_array(label,atom);
    skip_delim(file);
    get_string(&atom,file);
    pspring->atom_type[1]=index_in_array(label,atom);
    pspring->pfunc=GenericPlugIn<ArrayFunctionArray<Real> >::create(potplugin);
    if (!pspring->pfunc->read(file)) {
      ERRORQUIT("Error reading spring file.");
    }
    (*springlist) << pspring;
  }
}

void read_mass_file(Array<Real> *masses, istream &file, const Array<AutoString> &label) {
  masses->resize(label.get_size());
  for (int i=0; i<masses->get_size(); i++) {(*masses)(i)=-1.;}
  while (skip_delim(file)) {
    AutoString curlab;
    get_string(&curlab,file);
    Real m;
    file >> m;
    int i=index_in_array(label,curlab);
    if (i!=-1) {(*masses)(i)=m;}
  }
  for (int i=0; i<masses->get_size(); i++) {
    if ((*masses)(i)==-1.) {
      cerr << "Mass of " << label(i) << " is not known";
      ERRORQUIT("Aborting");
    }
  }
}

Real find_local_minimum(const Array<Real> &a,Real *pos=NULL) {
  int i=1;
  Real curmin=a(1);
  for (int j=1; j<a.get_size()-1; j++) {
    if (a(j+1)>a(j) && a(j-1)>a(j)) {
      if (a(j)<curmin) {
	i=j;
	curmin=a(j);
      }
    }
  }
  if (i==1 && a(1)>a(2)) {
    if (a(0)<a(1)) {
      i=0;
    }
    else {
      i=a.get_size()-1;
    }
    if (pos) (*pos)=(Real)i;
    return a(i);
  }
  Real l=(a(i+1)-a(i-1))/2.;
  Real q=(a(i+1)+a(i-1)-2.*a(i))/2.;
  Real x=-l/(2.*q);
  Real m=a(i)-sqr(l)/(4.*q);
  /*
  {
    ofstream file("tmp.tmp");
    for (int i=0; i<a.get_size(); i++) {
      file << a(i) << endl;
    }
  }
  cerr << ((Real)i)+x << endl;
  cerr << m << endl;
  cerr << i << endl << l << endl << q << endl << x << endl << m << endl;
  */
  if (pos) (*pos)=((Real)i)+x;
  return m;
}

extern char *helpstring;

GenericPlugIn<ArrayFunctionArray<Real> > *GenericPlugIn<ArrayFunctionArray<Real> >::list=NULL;

int main(int argc, char *argv[]) {
  char *latfilename="";
  char *strfilename="str.out";
  char *rel_strfilename="str_relax.out";
  char *sprfilename="slspring.out";
  char *potplugin="";
  char *massfilename="";
  Real bulkmod=0.;
  char *bulkmodfile="";
  Real maxs=0.05;
  int ns=1;
  Real T0=0.;
  Real T1=2000.;
  Real dT=100.;
  Real kppra=1000;
  rVector3d kshift(0.5,0.5,0.5);
  Real hplanck=6.6260755e-34/1.60217733e-19;  // h in eV s (not h bar);
  Real kboltzman=1.380658e-23/1.60217733e-19; // k_B in eV/K;
  Real convfk=1.60217733e-19*1e20; // converts eV/A^2 into J/m^2;
  Real mass_unit=1.6605402e-27;    // convert a.u. mass into kg;
  Real robust_len=0.;
  Real perfewer=1.;
  int peratom=0;
  int forceneg=0;
  char *kdispfile="";
  char *strainfilename="";
  int sigdig=5;
  int dummy=0;
  int dohelp=0;
  Real maxfklen=MAXFLOAT;
  if (file_exists("../Trange.in")) {
    ifstream tfile("../Trange.in");
    cerr << "Reading Trange.in file." << endl;
    T0=0;
    int nT=0;
    tfile >> T1 >> nT;
    dT=(T1-T0)/MAX(1,nT-1);
  }

  AskStruct options[]={
    {"","Vibrational free energy calculator using the Stiffness VS Length method " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-l","Input file defining the lattice (defaults: lat.in, ../lat.in, str.out)",STRINGVAL,&latfilename},
    {"-us","Input file defining the unrelaxed structure (default: str.out)",STRINGVAL,&strfilename},
    {"-rs","Input file defining the   relaxed structure (default: str_relax.out)",STRINGVAL,&rel_strfilename},
    {"-sp","Input file defining the springs (default: slspring.out)",STRINGVAL,&sprfilename},
    {"-m","Input file defining the atomic masses (default: ${atatdir}/data/masses.in)",STRINGVAL,&massfilename},
    {"-b","Bulk modulus",REALVAL,&bulkmod},
    {"-bf","Bulk modulus file",STRINGVAL,&bulkmodfile},
    {"-ms","Maximum stretching of the lattice parameter (default: 0.05)",REALVAL,&maxs},
    {"-ns","Number of lattice parameter stretching step (default: 1 => harmonic approximation)",INTVAL,&ns},
    {"-T0","Minimum temperature (default: 0) ",REALVAL,&T0},
    {"-T1","Maximum temperature (default: 2000)",REALVAL,&T1},
    {"-dT","Temperature step (default: 100)",REALVAL,&dT},
    {"-kp","Number of k-points per reciprocal atom (default: 1000)",REALVAL,&kppra},
    {"-sx","k-point shift (along 1st recip lat. vect.)",REALVAL,&kshift(0)},
    {"-sy","k-point shift (along 2nd recip lat. vect.)",REALVAL,&kshift(1)},
    {"-sz","k-point shift (along 3rd recip lat. vect.)",REALVAL,&kshift(2)},
    {"-hp","Planck's constant (default in (eV s))",REALVAL,&hplanck},
    {"-kb","Boltzman's constant (default in eV/K)",REALVAL,&kboltzman},
    {"-cfk","Conversion factor for force constants into energy/dist^2 (default: converts eV/A^2 into J/m^2)",REALVAL,&convfk},
    {"-mu","Mass units (default: converts a.u. mass into kg)",REALVAL,&mass_unit},
    {"-rl","Robust Length algorithm parameter for soft modes (beta)",REALVAL,&robust_len},
    {"-pa","Output free energy per atom instead of per unit cell",BOOLVAL,&peratom},
    {"-sc","Correction factor if spectator ion are present (default: 1)",REALVAL,&perfewer},
    {"-fn","Force continuation of calculations even if unstable",BOOLVAL,&forceneg},
    {"-df","Phonon dispersion curve calculation input file.",STRINGVAL,&kdispfile},
    {"-sf","Extra strain file",STRINGVAL,&strainfilename},
    {"-msl","Maximum spring length",REALVAL,&maxfklen},
    {"-sig","Number of significant digits to print in output files",INTVAL,&sigdig},
    {"-up","User-specified Potential Plug-in (optional)",STRINGVAL,&potplugin},
    {"-d","Use all default values",BOOLVAL,&dummy},
    {"-h","Display more help",BOOLVAL,&dohelp}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    return 1;
  }
  if (dohelp) {
    cout << helpstring;
    return 1;
  }

  int maxlpow=1;
  int dirdep=0;
  int concpow=0;

  Structure lat;
  Array<Arrayint> labellookup;
  Array<AutoString> label;
  rMatrix3d axes;
  {
    if (strlen(latfilename)==0) {
      if (file_exists("lat.in")) {latfilename="lat.in";}
      else if (file_exists("../lat.in")) {latfilename="../lat.in";}
      else if (file_exists("str.out")) {latfilename="str.out";}
    }
    ifstream latfile(latfilename);
    if (!latfile) ERRORQUIT("Unable to open lattice file");
    parse_lattice_file(&lat.cell, &lat.atom_pos, &lat.atom_type, &labellookup, &label, latfile, &axes);
    wrap_inside_cell(&lat.atom_pos,lat.atom_pos,lat.cell);
  }

  Structure str;
  {
    ifstream strfile(strfilename);
    if (!strfile) ERRORQUIT("Unable to open structure file");
    parse_structure_file(&str.cell, &str.atom_pos, &str.atom_type, label, strfile, NULL);
    wrap_inside_cell(&str.atom_pos,str.atom_pos,str.cell);
    rMatrix3d extrastrain;
    extrastrain.identity();
    if (strlen(strainfilename)>0) {
      ifstream strainfile(strainfilename);
      strainfile >> extrastrain;
    }
    str.cell=extrastrain*(str.cell);
    for (int at=0; at<str.atom_pos(); at++) {
      str.atom_pos(at)=extrastrain*(str.atom_pos(at));
    }
    lat.cell=extrastrain*(lat.cell);
    for (int at=0; at<lat.atom_pos(); at++) {
      lat.atom_pos(at)=extrastrain*(lat.atom_pos(at));
    }
  }

  LinkedList<SpringSvsL> springlist;
  LinkedList<ExtendableSpringSvsL> extspringlist;

  {
    ifstream sprfile(sprfilename);
    if (!sprfile) {
      AutoString filename("../");
      filename+=sprfilename;
      sprfile.clear();
      sprfile.open(filename);
      if (!sprfile) {
	ERRORQUIT("Unable to open spring file.");
      }
    }
    if (strlen(potplugin)>0) {
      if (!check_plug_in(ArrayFunctionArray<Real>(),potplugin)) {
	ERRORQUIT("Aborting");
      }
      read_spring_file(&extspringlist,sprfile,label,potplugin);
      dirdep=1;
    }
    else {
      read_spring_file(&springlist,&maxlpow,&dirdep,&concpow,sprfile,label);
    }
  }

  LinkedList<rMatrix3d> dirdep_mat;
  if (!dirdep) {
    rMatrix3d *pid=new rMatrix3d;
    pid->identity();
    dirdep_mat << pid;
  }
  else {
    SpaceGroup spacegroup;
    find_spacegroup(&spacegroup.point_op,&spacegroup.trans,lat.cell,lat.atom_pos,lat.atom_type);
    Array<rMatrix3d> pointgroup;
    pointgroup_from_spacegroup(&pointgroup, spacegroup.point_op);
    find_sym_const_tensor(&dirdep_mat, pointgroup,1);
  }

  Array2d<Real> redtonored;
  Array<Real> conc;
  Array<Real> red_conc;
  {
    Array2d<Real> m;
    Array<Real> v;
    calc_noredtored(&m,&v,lat,labellookup);
    calc_redtonored(&redtonored, m,v);
    //    if (!calc_concentration(&red_conc,lat,labellookup,str)) {
    //      ERRORQUIT("Error calculating concentration: cannot map an atom in structure onto lattice site.");
    //    }
    calc_concentration(&red_conc, str.atom_type,label.get_size());
    product(&conc,redtonored,red_conc);
  }

  MultiDimPoly multipoly;
  {
    int ddms=dirdep_mat.get_size();
    Array<Real> maxpows(1+ddms+redtonored.get_size()(0));
    Array<Real> dlineq(maxpows.get_size()),duineq(maxpows.get_size());
    Array<Real> cineq(maxpows.get_size());
    zero_array(&dlineq);
    zero_array(&duineq);
    zero_array(&cineq);
    maxpows(0)=maxlpow;
    for (int i=0; i<dirdep_mat.get_size(); i++) {
      maxpows(1+i)=1;
      dlineq(1+i)=-1;
      duineq(1+i)=1;
    }
    for (int i=0; i<redtonored.get_size()(0); i++) {
      maxpows(1+ddms+i)=concpow;
      cineq(1+ddms+i)=1;
    }
    LinkedList<LinearInequality> linineq;
    linineq << new LinearInequality(dlineq,-1) << new LinearInequality(duineq,1) << new LinearInequality(cineq,concpow);
    multipoly.init(maxpows,linineq);
  }

  Array<Real> atomic_masses;
  {
    ifstream massfile;
    if (strlen(massfilename)>0) {
      massfile.open(massfilename);
      if (!massfile) {
	ERRORQUIT("Unable to open atomic masses file");
      }
    }
    else {
      AutoString configfilename(getenv("HOME"));
      configfilename+="/.atat.rc";
      ifstream configfile(configfilename);
      if (!configfile) {
	ERRORQUIT("$HOME/.atat.rc was not found.");
      }
      while (configfile.get()!='=') {};
      skip_delim(configfile," \t");
      AutoString massfilename2;
      get_string(&massfilename2,configfile);
      massfilename2+="/data/masses.in";
      massfile.open(massfilename2);
      if (!massfile) {
	ERRORQUIT("Unable to open atomic masses file");
      }
    }
    read_mass_file(&atomic_masses,massfile,label);
  }

  remove_vacancies(&str,label);

  Structure rel_str;
  {
    ifstream rel_strfile(rel_strfilename);
    if (!rel_strfile) ERRORQUIT("Unable to open relaxed structure file");
    parse_structure_file(&rel_str.cell,&rel_str.atom_pos,&rel_str.atom_type,label,rel_strfile);
    wrap_inside_cell(&rel_str.atom_pos,rel_str.atom_pos,rel_str.cell);
    reorder_atoms(&rel_str,str);
  }

  Array<Real> masses(rel_str.atom_type.get_size());
  for (int at=0; at<rel_str.atom_type.get_size(); at++) {
    masses(at)=atomic_masses(rel_str.atom_type(at));
  }

  if (strlen(kdispfile)>0) {
    ifstream infile(kdispfile);
    if (!infile) {
      ERRORQUIT("Unable to open dispersion curve calculation input file.");
    }
    ofstream outffile("eigenfreq.out");
    ofstream outvfile("eigenvect.out");
    rMatrix3d recip=(!(~((rel_str.cell)*(!str.cell)*axes)));
    while (skip_delim(infile)) {
      int nbpts;
      rVector3d k1,k2;
      infile >> nbpts >> k1 >> k2;
      k1=recip*k1;
      k2=recip*k2;
      LinkedList<PairSpring> fklist;
      if (strlen(potplugin)>0) {
	make_nn_forcek(&fklist, str, rel_str, extspringlist, dirdep_mat, red_conc);
      }
      else {
	make_nn_forcek(&fklist, str, rel_str, springlist, dirdep_mat, conc, multipoly);
      }
      {
	LinkedListIterator<PairSpring> ifk(fklist);
	while (ifk) {
	  if (norm(ifk->dr)>maxfklen) {
	    delete fklist.detach(ifk);
	  }
	  else {
	    ifk++;
	  }
	}
      }
      LinkedList<Array<Real> > freq;
      LinkedList<Array2d<Complex> > eigvect;
      calc_dispersion_curve(&freq,&eigvect, k1,k2,nbpts,fklist,masses,convfk/mass_unit);
      LinkedListIterator<Array<Real> > i_freq(freq);
      LinkedListIterator<Array2d<Complex> > i_eigvect(eigvect);
      for (; i_freq; i_freq++, i_eigvect++) {
	for (int i=0; i<i_freq->get_size(); i++) {
	  outffile << (*i_freq)(i) << '\t';
	}
	outffile << endl;
	outvfile << (*i_eigvect);
      }
    }
  }
  else {
    iVector3d kmesh;
    calc_k_mesh(&kmesh,str.cell,kppra/(Real)(str.atom_pos.get_size()));
    int nT=(int)((T1-T0)/dT)+1;
    Array<Array<Real> > F;
    resize(&F,nT,ns);
    ofstream vdosfile("vdos.out");
    ofstream logfile("svsl.log");
    logfile.setf(ios::fixed);
    logfile.precision(sigdig);
    logfile << "kmesh= " << kmesh << endl;
    for (int is=0; is<ns; is++) {
      cerr << is+1 << "/" << ns << endl;
      Real iso_stretch=(1.+maxs*(Real)is/(Real)MAX(1,ns-1));
      Structure stretched_str;
      stretched_str.atom_type=rel_str.atom_type;
      stretched_str.cell=iso_stretch*rel_str.cell;
      stretched_str.atom_pos.resize(rel_str.atom_pos.get_size());
      for (int at=0; at<rel_str.atom_pos.get_size(); at++) {
	stretched_str.atom_pos(at)=iso_stretch*rel_str.atom_pos(at);
      }
      LinkedList<PairSpring> fklist;
      if (strlen(potplugin)>0) {
	make_nn_forcek(&fklist, str, stretched_str, extspringlist, dirdep_mat, red_conc);
      }
      else {
	make_nn_forcek(&fklist, str, stretched_str, springlist, dirdep_mat, conc, multipoly);
      }
      {
	LinkedListIterator<PairSpring> ifk(fklist);
	while (ifk) {
	  if (norm(ifk->dr)>maxfklen) {
	    delete fklist.detach(ifk);
	  }
	  else {
	    ifk++;
	  }
	}
      }
      logfile << "Number of NN springs: " << fklist.get_size() << endl;
      LinkedListIterator<PairSpring> ifk(fklist);
      for (; ifk; ifk++) {
	logfile << "atoms: " << ifk->whichatom[0] << " " << ifk->whichatom[1] << endl;
	logfile << "atom types: " << label(stretched_str.atom_type(ifk->whichatom[0])) << "-" << label(stretched_str.atom_type(ifk->whichatom[1])) << endl;
	logfile << "|dr|: " << norm(ifk->dr) << endl;
	logfile << "dr: " << ifk->dr << endl;
	logfile << "cell shift: " <<  ifk->cell_shift << endl;
	logfile << "fk matrix: " << endl << ifk->forcek << endl;
      }
      logfile << "end" << endl << endl;
      if (strlen(bulkmodfile)) {
	ifstream file(bulkmodfile);
	if (!file) ERRORQUIT("Unable to open bulk modulus file.");
	file >> bulkmod;
      }
      if (bulkmod==0.) bulkmod=calc_bulk_modulus(fklist,det(rel_str.cell));
      logfile << "Bulk modulus: " << bulkmod << endl;
      Real elasE=det(rel_str.cell)*bulkmod*pow(pow(iso_stretch,3.)-1.,2.)/2./(Real)(rel_str.atom_type.get_size());
      LinkedList<Real> freq,weight;
      if (!list_phonon_freq(&freq, &weight, kmesh, kshift, fklist, stretched_str.cell,masses,convfk/mass_unit)) {
	if (!forceneg && robust_len==0.) {
	  ERRORQUIT("Unstable modes found. Aborting.");
	}
      }
      
      LinkedListIterator<Real> ifreq(freq);
      logfile << "Phonon frequencies:" << endl;
      for ( ; ifreq; ifreq++) {
	logfile << *ifreq << endl;
      }
      logfile << "end" << endl << endl;
      
      Array<Real> freqa;
      LinkedList_to_Array(&freqa,freq);
      Real xmin,xmax;
      Array<Real> vdos;
      smooth_density(&xmin,&xmax,&vdos,freqa);
      Real df=(xmax-xmin)/(Real)(vdos.get_size()-1);
      LinkedListIterator<Real> iw(weight);
      Real w=(*iw)*(Real)(peratom ? 1:rel_str.atom_type.get_size())*perfewer;
      for (int i=0; i<vdos.get_size(); i++) {
	vdosfile << xmin+df*(Real)i << "\t" << w*vdos(i) << endl;
      }
      vdosfile << endl;

      if (is==0) {
	ofstream file("svib_ht");
	file << calc_vib_entropy(freq,weight)*(Real)(peratom ? 1:rel_str.atom_type.get_size())*perfewer << endl;
      }

      logfile << "Free energy" << endl;
      for (int it=0; it<nT; it++) {
	Real T=T0+(Real)it*dT;
	if (robust_len>0) {
	  Real geomass=1.;
	  for (int at=0; at<masses.get_size(); at++) {geomass*=masses(at);}
	  geomass=mass_unit*pow(geomass,1./(Real)(masses.get_size()));
	  Real L=pow(det(rel_str.cell)/(Real)(rel_str.atom_type.get_size())*exp(1.),1./3.)*sqrt(geomass/convfk)/2.;
	  F(it)(is)=calc_vib_free_energy_robust(freq,weight,kboltzman*T,hplanck,robust_len*L)*(Real)(peratom ? 1:rel_str.atom_type.get_size());
	}
	else {
	  F(it)(is)=(calc_vib_free_energy(freq,weight,kboltzman*T,hplanck)+elasE)*(Real)(peratom ? 1:rel_str.atom_type.get_size());
	}
	logfile << T << " " << iso_stretch << " " << F(it)(is) << endl;
      }
      logfile << "end" << endl << endl;
    }
    {
      ofstream outfile("svsl.out");
      ofstream fvibfile("fvib");
      for (int it=0; it<nT; it++) {
	Real pos=0;
	Real minF=(ns==1 ? F(it)(0) : find_local_minimum(F(it),&pos));
	outfile << (T0+(Real)it*dT) << " " << minF*perfewer << " " << pos*maxs/(Real)MAX(1,ns-1) << endl;
	fvibfile << minF*perfewer << endl;
      }
    }
  }
}
