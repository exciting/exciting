#include <fstream.h>
#include <strstream.h>
#include <iomanip.h>
#include <unistd.h>
#include <sys/stat.h>
#include "parse.h"
#include "mrefine.h"
#include "getvalue.h"

#include "version.h"

extern char *helpstring;

int sigdig=6;

// Creates a directory containing a file describing a new structure.;
// ARG:
//  str: (IN) structure to write out;
//  lat: (IN)               \ ;
//  site_type_list: (IN)    | describes the lattice (see parse.h);
//  atom_label: (IN)        / ;
//  axes: (IN) coordinate system in which to express coordinates;
void write_structure_file(const StructureInfo &str, const Structure &lat,
			  const Array<Arrayint> &site_type_list, const Array<AutoString> &atom_label, const rMatrix3d &axes) {
  mkdir(str.label,S_IRWXU | S_IRWXG | S_IRWXO);
  if (chdir(str.label)==0) {
    ofstream file("str.out");
    file.setf(ios::fixed);
    file.precision(sigdig);
    write_structure(str,lat,site_type_list,atom_label, axes, file);
    ofstream waitfile("wait");
    chdir("..");
  }
}

// Writes all the output files describing the current state;
// of the cluster expansion (see class CEFitInfo in raffine.h).
void write_fit_info(const CEFitInfo &fitinfo, const Structure &lattice,
		    const Array<Arrayint> &site_type_list, const Array<AutoString> &atom_label, const rMatrix3d &axes, const char *propnameext="eci.out") {
  ofstream log("maps.log");
  log << "Maps version " MAPS_VERSION << endl;
  if (fitinfo.status & CEFitInfo::fit_impossible) {
    log << "Not enough known energies to fit CE" << endl;
  }
  else {
    log << "The internal database of structures extends up to " << fitinfo.max_volume*lattice.atom_type.get_size() << " atoms/unit cell, see predstr.out" << endl;
    if (fitinfo.status & CEFitInfo::gs_problem) {
      log << "Among structures of known energy, true ground states differ from fitted ground states" << endl;
    }
    else {
      log << "Among structures of known energy, true and predicted ground states agree" << endl;
    }
    if (fitinfo.status & CEFitInfo::new_gs) {
      log << "New ground states with at most " << fitinfo.max_volume*lattice.atom_type.get_size() << " atoms/unit cell predicted , see predstr.out" << endl;
    }
    if (fitinfo.status==CEFitInfo::fit_ok) {
      log << "No other ground states of " << fitinfo.max_volume*lattice.atom_type.get_size() << " atoms/unit cell or less exist." << endl;
    }
    // insert range todo;
    log << "Crossvalidation score: " << fitinfo.cv << endl;
    ofstream fit_file("fit.out");
    fit_file.setf(ios::fixed);
    fit_file.precision(sigdig);
    for (int i=0; i<fitinfo.energy.get_size(); i++) {
      for (int j=0; j<fitinfo.concentration(i).get_size(); j++) { 
	fit_file << fitinfo.concentration(i)(j) << " ";
      }
      fit_file << fitinfo.energy(i)<< " " 
               << fitinfo.fitted_energy(i)<< " " 
               << fitinfo.energy(i)-fitinfo.fitted_energy(i)<< " " 
	       << fitinfo.weight(i) << " " 
	       << fitinfo.pstr(i)->label << endl;
    }
    ofstream gs_file("gs.out");
    gs_file.setf(ios::fixed);
    gs_file.precision(sigdig);
    for (int i=0; i<fitinfo.gs_list.get_size(); i++) {
      int gs=fitinfo.gs_list(i);
      for (int j=0; j<fitinfo.concentration(gs).get_size(); j++) {
	gs_file << fitinfo.concentration(gs)(j) << " ";
      }
      gs_file << fitinfo.energy(gs) << " " 
	      << fitinfo.fitted_energy(gs) << " " 
	      << fitinfo.energy(gs)-fitinfo.fitted_energy(gs) << " " 
	      << fitinfo.pstr(gs)->label << endl;
    }

    {
      ofstream gs_connect("gs_connect.out");
      LinkedListIterator<PolytopeFace> i(fitinfo.hull);
      for (; i; i++) {
	for (int j=0; j<i->pts.get_size(); j++) {
	  gs_connect << fitinfo.pstr(i->pts(j))->label << " ";
	}
	gs_connect << endl;
      }
    }

    {
      // print chemical potentials;
      int dimlong=fitinfo.conc_to_fullconc.get_size()(0);
      int dimshort=fitinfo.conc_to_fullconc.get_size()(1);
      Array2d<Real> CC,CCi,A,CCiAt;
      Array2d<Real> mu_to_fullmu;
      inner_product(&CC,fitinfo.corr_to_fullconc,fitinfo.corr_to_fullconc);
      invert_matrix(&CCi,CC);
      A.resize(fitinfo.corr_to_conc.get_size()(0),fitinfo.corr_to_conc.get_size()(1)+1);
      zero_array(&A);
      for (int i=0; i<fitinfo.corr_to_conc.get_size()(0); i++) {
	for (int j=0; j<fitinfo.corr_to_conc.get_size()(1); j++) {
	  A(i,1+j)=fitinfo.corr_to_conc(i,j);
	}
      }
      outer_product(&CCiAt,CCi,A);
      product(&mu_to_fullmu,fitinfo.corr_to_fullconc,CCiAt);
      ofstream mufile("chempot.out");
      Array<Array<Real> > gs_mu(fitinfo.gs_list.get_size());
      for (int j=0; j<fitinfo.gs_list.get_size(); j++) {
        gs_mu(j).resize(dimlong);
	zero_array(&gs_mu(j));
      }
      Array<int> gs_nb_mu(fitinfo.gs_list.get_size());
      zero_array(&gs_nb_mu);
      mufile << "chemical potentials stabilizing multiphase equilibria" << endl;
      mufile << "ground state indices , structure labels : chemical potentials" << endl;
      LinkedListIterator<PolytopeFace> i(fitinfo.hull);
      for (; i; i++) {
        Array<Real> full_mu;
	Array<Real> short_mu(i->up.get_size()-1);
	Real e=i->up(short_mu.get_size());
	for (int j=0; j<short_mu.get_size(); j++) {
	  short_mu(j)=-(i->up(j))/e;
	}
	product(&full_mu,mu_to_fullmu,short_mu);
	for (int j=0; j<full_mu.get_size(); j++) {
	  full_mu(j)+=(i->c)/e;
	}
	/*
	for (int j=0; j<i->pts.get_size(); j++) {
	  cerr << j << " " << fitinfo.fitted_energy(i->pts(j))-inner_product(full_mu,fitinfo.concentration(i->pts(j))) << " ";
	}
	cerr << endl;
	*/
	for (int j=0; j<i->pts.get_size(); j++) {
          int gs=index_in_array(fitinfo.gs_list,i->pts(j));
          sum(&gs_mu(gs),gs_mu(gs),full_mu);
	  gs_nb_mu(gs)++;
	  mufile << gs << " ";
	}
	mufile << " , ";
	for (int j=0; j<i->pts.get_size(); j++) {
	  mufile << fitinfo.pstr(i->pts(j))->label << " ";
	}
        mufile << ": ";
        for (int j=0; j<full_mu.get_size(); j++) {
          mufile << full_mu(j) << " ";
        }
        mufile << endl;
      }
      mufile << "chemical potentials stabilizing each ground state" << endl;
      mufile << "ground state index , structure label : chemical potentials" << endl;
      for (int j=0; j<gs_mu.get_size(); j++) {
        mufile << j << " , " << fitinfo.pstr(fitinfo.gs_list(j))->label << " : ";
	for (int i=0; i<gs_mu(j).get_size(); i++) {
	  mufile << gs_mu(j)(i)/(Real)(gs_nb_mu(j)) << " ";
	}
        mufile << endl;
      }
    }

    ofstream eci_file(propnameext);
    eci_file.setf(ios::fixed);
    eci_file.precision(sigdig);
    ofstream cluster_file("clusters.out");
    cluster_file.setf(ios::fixed);
    cluster_file.precision(sigdig);
    rMatrix3d iaxes=!axes;
    for (int i=0; i<fitinfo.eci.get_size(); i++) {
      eci_file << fitinfo.eci(i) << endl;
      cluster_file << (int)(fitinfo.multiplicity(i)) << endl;
      cluster_file << get_cluster_length(fitinfo.cluster(i).clus) << endl;
      cluster_file << fitinfo.cluster(i).clus.get_size() << endl;
      for (int j=0; j<fitinfo.cluster(i).clus.get_size(); j++) {
	cluster_file << (iaxes*fitinfo.cluster(i).clus(j)) << " " << fitinfo.cluster(i).site_type(j) << " " << fitinfo.cluster(i).func(j) << endl;
      }
      cluster_file << endl;
    }
    ofstream gsstr_file("gs_str.out");
    gsstr_file.setf(ios::fixed);
    gsstr_file.precision(sigdig);
    for (int i=0; i<fitinfo.gs_list.get_size(); i++) {
      write_structure(*fitinfo.pstr(fitinfo.gs_list(i)),lattice,site_type_list,atom_label, axes, gsstr_file);
      gsstr_file << "end" << endl << endl;
    }
  }
}

int update_structure(StructureInfo *pstr, Real atom_factor, const char *propname="energy") {
  int changed=0;
  StructureInfo::Status new_status;
  Real new_energy;
  if (MyMPIobj.file_exists("error")) {
    new_status=StructureInfo::error;
  }
  else if (MyMPIobj.file_exists(propname)) {
    new_status=StructureInfo::calculated;
    new_energy=MAXFLOAT;
    if (MyMPIobj.is_root()) {
      ifstream file(propname);
      if (file) {
	file >> new_energy;
      }
    }
    MyMPI_BcastStream(&new_energy);
    if (new_energy!=MAXFLOAT) {
      // if the file energy contains a valid number use it as energy;
      if (atom_factor>0) {
	new_energy*=atom_factor/(Real)(pstr->atom_pos.get_size());
      }
    }
    else {
      // if anythings goes wrong: error;
      new_status=StructureInfo::error;	  }
  }
  else {
    // if neither error nor calculated => is being calculated;
    new_status=StructureInfo::busy;
  }
  // a structure has changed if its status or energy has changed;
  if (new_status!=pstr->status) changed=1;
  if (new_status==StructureInfo::calculated && pstr->status==StructureInfo::calculated) {
    if (new_energy!=pstr->energy) changed=1;
  }
  // update values in memory;
  pstr->status=new_status;
  pstr->energy=new_energy;
  // indicates if the structure has changed;
  return changed;
}

int update_all_structures(StructureBank<StructureInfo> *pstr_bank, const Structure &lattice, const Array<Arrayint> &site_type_list, const Array<AutoString> &atom_label, char *propname="energy", int intensive=0) {
  Real atom_factor=(intensive ? -1 : lattice.atom_pos.get_size());
  int changed=0;
  // list structures on disk;
  LinkedList<AutoString> label_on_disk;
  if (MyMPIobj.is_root()) {
    system("ls */str.out 2> /dev/null | sed 's+/str.out++g' > strlist.out");
    ifstream strfile("strlist.out");
    while (strfile && !strfile.eof()) {
      char buf[MAX_LINE_LEN];
      buf[0]=0;
      strfile.get(buf,MAX_LINE_LEN-1);
      char tmp;
      strfile.get(tmp);
      if (strlen(buf)==0) break;
      label_on_disk << new AutoString(buf);
    }
    unlink("strlist.out"); // cleanup structure list file;
  }
  MyMPI_BcastStream(&label_on_disk);

  LinkedListIterator<StructureInfo> i(pstr_bank->get_structure_list());
  for (; i; i++) {
    if (!(i->status & StructureInfo::unknown)) {
      LinkedListIterator<AutoString> i_disk(label_on_disk);
      for ( ; i_disk; i_disk++) {
	if (i->label==*i_disk) break;
      }
      if (!i_disk) {
	i->status=StructureInfo::unknown;
	changed=1;
      }
      else {
	if (chdir(i->label)==0) {
	  if (update_structure(&(*i),atom_factor,propname)) changed=1;
	  delete label_on_disk.detach(i_disk);
	  chdir("..");
	}
      }
    }
  }

  LinkedListIterator<AutoString> i_disk(label_on_disk);
  for (; i_disk; i_disk++) {
    int str_is_ok=0;
    StructureInfo str;
    if (MyMPIobj.is_root()) {
      if (chdir(*i_disk)==0) {
	ifstream strfile("str.out");
	if (strfile) {
	  str_is_ok=parse_structure_file(&str.cell,&str.atom_pos,&str.atom_type,atom_label,strfile);
	}
	chdir("..");
      }
      else {
	cerr << "Unable to cd to " << *i_disk << endl;
      }
      if (!str_is_ok) cerr << "Error while reading structure " << *i_disk << endl;
    }
    MyMPI_BcastStream(&str_is_ok);
    MyMPI_BcastStream(&str);

    if (str_is_ok) {
      if (str.atom_pos.get_size()==0) {
	cerr << "Problem reading structure " << *i_disk << endl;
	ERRORQUIT("Aborting.");
      }
      if (fix_atom_type(&str,lattice,site_type_list,0)) {
	StructureInfo *pstr;
	if (pstr_bank->add_structure(str,&pstr)) {
	  pstr->label.set(*i_disk);
	  if (update_structure(pstr,atom_factor,propname)) changed=1;
	}
	else {
	  if (pstr->status & StructureInfo::unknown) {
	    pstr->label.set(*i_disk);
	    if (update_structure(pstr,atom_factor,propname)) changed=1;
	  }
	  else {
	    if (MyMPIobj.is_root()) cerr << "Structure " << *i_disk <<" not loaded since it is the same as structure " << pstr->label << endl;
	  }
	}
      }
      else {
	if (MyMPIobj.is_root()) cerr << "Unable to read structure " << *i_disk << " due to geometric incompatibilities with the lattice." << endl;
      }
    }
  }
  return changed;
}

#define MAPS_IS_RUNNING "maps_is_running"

int main(int argc, char *argv[]) {
  MyMPIobj.init(argc,argv);
  // parse command line arguments or display help (see getvalue.h);
  int dohelp=0;
  char *latticefilename="lat.in";
  int polltime=10;
  int max_multiplet=4;
  int dummy=0;
  int quiet=0;
  int gsnbatom=0;
  char *conc_file="crange.in";
  int do2D=0;
  Real complexity_exp=3.;
  Real high_energy=MAXFLOAT;
  char *predictor_labels="";
  char *algo_label="std";
  char *corrfunc_label="trigo";
  char *propname="energy";
  int ignore_gs=0;
  int intensive=0;
  AskStruct options[]={
    {"","MIT Ab initio Phase Stability (MAPS) code " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-h","Display more help",BOOLVAL,&dohelp},
    {"-l","Input file defining the lattice (default: lat.in)",STRINGVAL,&latticefilename},
    {"-z","Tolerance for finding symmetry operations (default: 1e-3)",REALVAL,&zero_tolerance},
    {"-c","Exponent of the order of complexity (default: 3)",REALVAL,&complexity_exp},
    {"-t","Time between disk reads in sec (default: 10 sec)",INTVAL,&polltime},
    {"-m","Maximum number of points in cluster (default 4)",INTVAL,&max_multiplet},
    {"-g","Extend ground state search up to structures having at least that many atoms.",INTVAL,&gsnbatom},
    {"-cr","Concentration region input file",STRINGVAL,&conc_file},
    {"-he","Highest predicted energy, above ground state hull, allowed when generating structures (default: no limit)",REALVAL,&high_energy},
    {"-2d","Find supercells along a and b axes only",BOOLVAL,&do2D},
    {"-p","Predictor plugins to use (examples: -p=es or -p=es_cs)",STRINGVAL,&predictor_labels},
    {"-ks","same as -p",STRINGVAL,&predictor_labels},
    {"-fa","Select fitting algorithm (default: built-in)",STRINGVAL,&algo_label},
    {"-crf","Select correlation functions (default: trigo)",STRINGVAL,&corrfunc_label},
    {"-pn","Property to cluster expand (default: energy)",STRINGVAL,&propname},
    {"-ig","Ignore whether cluster expansion predicts correct ground states",BOOLVAL,&ignore_gs},
    {"-pa","Quantity to expand is already per atom",BOOLVAL,&intensive},
    {"-q","Quiet mode (do not print status to stderr)",BOOLVAL,&quiet},
    {"-sig","Number of significant digits to print in output files",INTVAL,&sigdig},
    {"-d","Use all default values",BOOLVAL,&dummy}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    return 1;
  }
  if (dohelp) {
    cout << helpstring;
    return 1;
  }

  if (MyMPIobj.file_exists(MAPS_IS_RUNNING)) {
    ERRORQUIT("Maps is already running in this directory. To override this error, type rm " MAPS_IS_RUNNING);
  }
  if (MyMPIobj.is_root()) {
    ofstream tmp(MAPS_IS_RUNNING);
  }
  
  // read in lattice (see parse.h);
  Structure lat;
  Array<Arrayint> site_type_list;
  Array<AutoString> atom_label;
  rMatrix3d axes;
  if (MyMPIobj.is_root()) {
    ifstream file(latticefilename);
    if (!file) ERRORQUIT("Unable to open lattice file.");
    parse_lattice_file(&lat.cell, &lat.atom_pos, &lat.atom_type, &site_type_list, &atom_label, file, &axes);
    if (fabs(det(lat.cell))<zero_tolerance) ERRORQUIT("Lattice vectors are coplanar.");
    wrap_inside_cell(&lat.atom_pos,lat.atom_pos,lat.cell);
    {
      ofstream labelfile("atoms.out");
      for (int i=0; i<atom_label.get_size(); i++) {
	labelfile << atom_label(i) << endl;
      }
    }
  }
  MyMPI_BcastStream(&lat);
  MyMPI_BcastStream(&site_type_list);
  MyMPI_BcastStream(&atom_label);
  MyMPI_BcastStream(&axes);

  // check number of components;
  int maxt=0;
  for (int i=0; i<site_type_list.get_size(); i++) {
    maxt=MAX(site_type_list(i).get_size(),maxt);
  }
  if (maxt==1) {
    ERRORQUIT("Nothing to cluster expand!");
  }
  // find space group (see findsym.h);
  SpaceGroup spacegroup;
  spacegroup.cell=lat.cell;
  find_spacegroup(&spacegroup.point_op,&spacegroup.trans,lat.cell,lat.atom_pos,lat.atom_type);
  if (contains_pure_translations(spacegroup.point_op,spacegroup.trans)) {
    if (MyMPIobj.is_root()) cerr << "Warning: unit cell is not primitive." << endl;
  }

  // initialize a table of correlation functions;
  if (!check_plug_in(CorrFuncTable(),corrfunc_label)) {
    ERRORQUIT("Aborting");
  }
  CorrFuncTable *pcorrfunc=GenericPlugIn<CorrFuncTable>::create(corrfunc_label);

  // initialize a cluster expansion object (with user-specified algorithm, if requested);
  if (!check_plug_in(ClusterExpansionCreator(),algo_label)) {
    ERRORQUIT("Aborting");
  }
  ClusterExpansionCreator *pcec=GenericPlugIn<ClusterExpansionCreator>::create(algo_label);
  ClusterExpansion *pce=pcec->create(lat,site_type_list,atom_label,spacegroup,pcorrfunc);
  delete pcec;

  // initialize cluster expansion object parameters;
  pce->max_multiplet=max_multiplet; // sets the maximum number of point per cluster;
  pce->complexity_exp=complexity_exp;
  pce->set_highest_energy_allowed(high_energy);
  LinkedList<LinearInequality> ineq;
  if (strlen(conc_file)>0) {
    ifstream crangefile(conc_file);
    if (crangefile) {
      read_inequalities(&ineq,atom_label,crangefile);
    }
  }
  pce->set_concentration_range(ineq,ignore_gs);
  if (!check_plug_in(EnergyPredictor(),predictor_labels)) {
    ERRORQUIT("Aborting");
  }
  pce->set_predictor_labels(predictor_labels);
  pce->access_structure_bank().set_2D_mode(do2D);

  if (gsnbatom) {
    if (!quiet && MyMPIobj.is_root()) cerr << "Generating structures for ground state search." << endl;
    int maxvol=(gsnbatom+lat.atom_pos.get_size()-1)/lat.atom_pos.get_size();
    for (int i=1; i<=maxvol; i++) {
      if (!quiet && MyMPIobj.is_root()) cerr << "Volume= " << i*lat.atom_pos.get_size() << "/" << gsnbatom << endl;
      pce->access_structure_bank().find_new_structures(i);
    }
    if (!quiet && MyMPIobj.is_root()) cerr << "done!" << endl;
  }

  // clear up fit status file;
  CEFitInfo fitinfo;
  fitinfo.status=CEFitInfo::fit_impossible;
  ostrstream propnameext;
  if (strcmp(propname,"energy")==0) {
    propnameext << "eci.out" << '\0';
  }
  else {
    propnameext << propname << ".eci" << '\0';
  }
  if (MyMPIobj.is_root()) write_fit_info(fitinfo,lat,site_type_list,atom_label,axes,propnameext.str());

  // main loop: wait for "events";
  while (1) {
    MyMPIobj.barrier();
    // look for updated structure status/energy;
    while (update_all_structures(&(pce->access_structure_bank()),lat,site_type_list,atom_label,propname,intensive) || MyMPIobj.file_exists("refresh")) {
      while (MyMPIobj.file_exists("refresh") && MyMPIobj.is_root()) unlink("refresh");
      // if there are changes, fit a new CE;
      CEFitInfo fitinfo;
      if (!quiet && MyMPIobj.is_root()) cerr << "Finding best cluster expansion..." << endl;
      pce->find_best_cluster_choice(&fitinfo); //fit;
      if (MyMPIobj.is_root()) write_fit_info(fitinfo,lat,site_type_list,atom_label,axes,propnameext.str()); //print results;
      if (!quiet && MyMPIobj.is_root()) cerr << "done!" << endl;
    }
    // check if user wants to abort;
    if (MyMPIobj.file_exists("stop")) break;
    // check if user wants a new structure;
    if (MyMPIobj.file_exists("ready")) {
      StructureInfo *pstr;
      if (!quiet && MyMPIobj.is_root()) cerr << "Finding best structure..." << endl;
      // first try to find the best structure;
      pstr=pce->find_best_structure();
      // if we don't have enough info yet, find the next structure needed to fit a minimal CE;
      if (!pstr) pstr=pce->find_initial_structures();
      // if we don't have enough info yet, just find the cheapest structure to compute;
      if (!pstr) pstr=pce->find_first_unknown_structure();
      if (MyMPIobj.is_root()) {
	write_structure_file(*pstr,lat,site_type_list,atom_label,axes); // write structure;
	while (file_exists("ready") && MyMPIobj.is_root()) unlink("ready"); // tell user we're done;
      }
      pstr->status=StructureInfo::busy; // mark structure as busy;
      if (!quiet && MyMPIobj.is_root()) cerr << "done!" << endl;
    }
    for (int t=0; t<polltime; t++) {
      if (MyMPIobj.file_exists("refresh")) break;
      sleep(1); // wait a little;
    }
  }
  if (MyMPIobj.is_root()) {
    unlink("stop");
    unlink(MAPS_IS_RUNNING);
  }
  delete pce;
  delete pcorrfunc;
  MyMPIobj.barrier();
}
