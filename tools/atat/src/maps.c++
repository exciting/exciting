#include <fstream.h>
#include <strstream.h>
#include <iomanip.h>
#include <unistd.h>
#include <sys/stat.h>
#include "parse.h"
#include "refine.h"
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
		    const Array<Arrayint> &site_type_list, const Array<AutoString> &atom_label, const rMatrix3d &axes) {
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
    log <<  "Concentration range used for ground state checking: [" << (1.+fitinfo.minc_gs_ok)/2. << "," << (1.+fitinfo.maxc_gs_ok)/2. << "]." << endl;
    log << "Crossvalidation score: " << fitinfo.cv << endl;
    ofstream fit_file("fit.out");
    fit_file.setf(ios::fixed);
    fit_file.precision(sigdig);
    for (int i=0; i<fitinfo.energy.get_size(); i++) {
      fit_file << (1.+fitinfo.concentration(i))/2.
               << " " << fitinfo.energy(i)
               << " " << fitinfo.fitted_energy(i)
               << " " << fitinfo.energy(i)-fitinfo.fitted_energy(i)
	       << " " << fitinfo.weight(i)
	       << " " << fitinfo.pstr(i)->label << endl;
    }
    ofstream gs_file("gs.out");
    gs_file.setf(ios::fixed);
    gs_file.precision(sigdig);
    for (int i=0; i<fitinfo.gs_list.get_size(); i++) {
      gs_file << (1.+fitinfo.concentration(fitinfo.gs_list(i)))/2.
	      << " " << fitinfo.energy(fitinfo.gs_list(i))
	      << " " << fitinfo.fitted_energy(fitinfo.gs_list(i))
	      << " " << fitinfo.pstr(fitinfo.gs_list(i))->label << endl;
    }
    ofstream eci_file("eci.out");
    eci_file.setf(ios::fixed);
    eci_file.precision(sigdig);
    ofstream cluster_file("clusters.out");
    cluster_file.setf(ios::fixed);
    cluster_file.precision(sigdig);
    rMatrix3d iaxes=!axes;
    for (int i=0; i<fitinfo.eci.get_size(); i++) {
      eci_file << fitinfo.eci(i) << endl;
      cluster_file << (int)(fitinfo.multiplicity(i)) << endl;
      cluster_file << get_cluster_length(fitinfo.cluster(i)) << endl;
      cluster_file << fitinfo.cluster(i).get_size() << endl;
      for (int j=0; j<fitinfo.cluster(i).get_size(); j++) {
	cluster_file << (iaxes*fitinfo.cluster(i)(j)) << endl;
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

int update_structure(StructureInfo *pstr, Real atom_factor) {
  int changed=0;
  StructureInfo::Status new_status;
  Real new_energy;
  if (file_exists("energy")) {
    new_energy=MAXFLOAT;
    ifstream file("energy");
    if (file) {
      file >> new_energy;
    }
    if (new_energy!=MAXFLOAT) {
      // if the file energy contains a valid number use it as energy;
      new_energy*=atom_factor/(Real)(pstr->atom_pos.get_size());
      new_status=StructureInfo::calculated;
    }
    else {
      // if anythings goes wrong: error;
      new_energy=0.;
      new_status=StructureInfo::error;
    }
  } else {
    // if no energy file => is being calculated;
    new_energy=0.;
    new_status=StructureInfo::busy;
  }
  if (file_exists("error")) {
    new_status=StructureInfo::error;
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

int update_all_structures(StructureBank<StructureInfo> *pstr_bank, const Structure &lattice, const Array<Arrayint> &site_type_list, const Array<AutoString> &atom_label) {
  int changed=0;
  // list structures on disk;
  system("ls */str.out 2> /dev/null | sed 's+/str.out++g' > strlist.out");
  ifstream strfile("strlist.out");
  LinkedList<AutoString> label_on_disk;
  while (strfile && !strfile.eof()) {
    char buf[MAX_LINE_LEN];
    buf[0]=0;
    strfile.get(buf,MAX_LINE_LEN-1);
    char tmp;
    strfile.get(tmp);
    if (strlen(buf)==0) break;
    label_on_disk << new AutoString(buf);
  }

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
	  if (update_structure(&(*i),lattice.atom_pos.get_size())) changed=1;
	  delete label_on_disk.detach(i_disk);
	  chdir("..");
	}
      }
    }
  }

  LinkedListIterator<AutoString> i_disk(label_on_disk);
  for (; i_disk; i_disk++) {
    if (chdir(*i_disk)==0) {
      StructureInfo str;
      ifstream strfile("str.out");
      if (strfile) {
	if (parse_structure_file(&str.cell,&str.atom_pos,&str.atom_type,atom_label,strfile)) {
	  if (str.atom_pos.get_size()==0) {
	    cerr << "Problem reading structure " << *i_disk << endl;
	    ERRORQUIT("Aborting.");
	  }
	  if (fix_atom_type(&str,lattice,site_type_list,0)) {
	    StructureInfo *pstr;
	    if (pstr_bank->add_structure(str,&pstr)) {
	      pstr->label.set(*i_disk);
	      if (update_structure(pstr,lattice.atom_pos.get_size())) changed=1;
	    }
	    else {
	      if (pstr->status & StructureInfo::unknown) {
		pstr->label.set(*i_disk);
		if (update_structure(pstr,lattice.atom_pos.get_size())) changed=1;
	      }
	      else {
		cerr << "Structure " << *i_disk <<" not loaded since it is the same as structure " << pstr->label << endl;
	      }
	    }
	  }
	  else {
	    cerr << "Error while reading structure " << *i_disk << endl;
	  }
	}
	else {
	  cerr << "Error while reading structure " << *i_disk << endl;
	}
      }
      chdir("..");
    }
    else {
      cerr << "Unable to cd to " << *i_disk << endl;
    }
  }
  unlink("strlist.out"); // cleanup structure list file;
  return changed;
}

#define MAPS_IS_RUNNING "maps_is_running"

int main(int argc, char *argv[]) {
  // parse command line arguments or display help (see getvalue.h);
  int dohelp=0;
  char *latticefilename="lat.in";
  int polltime=10;
  int max_multiplet=4;
  int dummy=0;
  int quiet=0;
  int gsnbatom=0;
  Real minc_gs_ok=0.;
  Real maxc_gs_ok=1.;
  int do2D=0;
  Real complexity_exp=3.;
  char *predictor_labels="";
  char *algo_label="std";
  AskStruct options[]={
    {"","MIT Ab initio Phase Stability (MAPS) code " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-h","Display more help",BOOLVAL,&dohelp},
    {"-l","Input file defining the lattice (default: lat.in)",STRINGVAL,&latticefilename},
    {"-z","Tolerance for finding symmetry operations (default: 1e-3)",REALVAL,&zero_tolerance},
    {"-c","Exponent of the order of complexity (default: 3)",REALVAL,&complexity_exp},
    {"-t","Time between disk reads in sec (default: 10 sec)",INTVAL,&polltime},
    {"-m","Maximum number of points in cluster (default 4)",INTVAL,&max_multiplet},
    {"-g","Extend ground state search up to structures having at least that many atoms.",INTVAL,&gsnbatom},
    {"-c0","[c0,c1] is the concentration range where ground states must be correct.",REALVAL,&minc_gs_ok},
    {"-c1","",REALVAL,&maxc_gs_ok},
    {"-2d","Find supercells along a and b axes only",BOOLVAL,&do2D},
    {"-p","Predictor plugins to use (examples: -p=cs or -p=cs_el)",STRINGVAL,&predictor_labels},
    {"-ks","same as -p",STRINGVAL,&predictor_labels},
    {"-fa","Select fitting algorithm (default: built-in)",STRINGVAL,&algo_label},
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

  if (file_exists(MAPS_IS_RUNNING)) {
    ERRORQUIT("Maps is already running in this directory. To override this error, type rm " MAPS_IS_RUNNING);
  }
  {
    ofstream tmp(MAPS_IS_RUNNING);
  }
  
  // read in lattice (see parse.h);
  Structure lat;
  Array<Arrayint> site_type_list;
  Array<AutoString> atom_label;
  rMatrix3d axes;
  ifstream file(latticefilename);
  if (!file) ERRORQUIT("Unable to open lattice file.");
  parse_lattice_file(&lat.cell, &lat.atom_pos, &lat.atom_type, &site_type_list, &atom_label, file, &axes);
  if (fabs(det(lat.cell))<zero_tolerance) ERRORQUIT("Lattice vectors are coplanar.");
  wrap_inside_cell(&lat.atom_pos,lat.atom_pos,lat.cell);

  // is this a pseudo-binary system?;
  int nb_of_binary_site=0;
  for (int i=0; i<site_type_list.get_size(); i++) {
    if (site_type_list(i).get_size()>2)
      ERRORQUIT("For multicomponent systems use mmaps instead.");
    if (site_type_list(i).get_size()==2) nb_of_binary_site++;
  }
  if (nb_of_binary_site>1)
    ERRORQUIT("For coupled cluster expansions use mmaps instead.");
  if (nb_of_binary_site==0)
    ERRORQUIT("Nothing to cluster expand!");

  // find space group (see findsym.h);
  SpaceGroup spacegroup;
  spacegroup.cell=lat.cell;
  find_spacegroup(&spacegroup.point_op,&spacegroup.trans,lat.cell,lat.atom_pos,lat.atom_type);
  if (contains_pure_translations(spacegroup.point_op,spacegroup.trans)) {
    cerr << "Warning: unit cell is not primitive." << endl;
  }
  
  // initialize a cluster expansion object (with user-specified algorithm, if requested);
  if (!check_plug_in(ClusterExpansionCreator(),algo_label)) {
    ERRORQUIT("Aborting");
  }
  ClusterExpansionCreator *pcec=GenericPlugIn<ClusterExpansionCreator>::create(algo_label);
  ClusterExpansion *pce=pcec->create(lat,site_type_list,atom_label,spacegroup);
  delete pcec;

  // initialize cluster expansion object parameters;
  pce->max_multiplet=max_multiplet; // sets the maximum number of point per cluster;
  pce->complexity_exp=complexity_exp;
  pce->set_concentration_range_gs_ok(2.*minc_gs_ok-1.,2.*maxc_gs_ok-1.);
  if (!check_plug_in(EnergyPredictor(),predictor_labels)) {
    ERRORQUIT("Aborting");
  }
  pce->set_predictor_labels(predictor_labels);
  pce->access_structure_bank().set_2D_mode(do2D);

  if (gsnbatom) {
    if (!quiet) cerr << "Generating structures for ground state search." << endl;
    int maxvol=(gsnbatom+lat.atom_pos.get_size()-1)/lat.atom_pos.get_size();
    for (int i=1; i<=maxvol; i++) {
      if (!quiet) cerr << "Volume= " << i*lat.atom_pos.get_size() << "/" << gsnbatom << endl;
      pce->access_structure_bank().find_new_structures(i);
    }
    if (!quiet) cerr << "done!" << endl;
  }

  // clear up fit status file;
  CEFitInfo fitinfo;
  fitinfo.status=CEFitInfo::fit_impossible;
  write_fit_info(fitinfo,lat,site_type_list,atom_label,axes);

  // main loop: wait for "events";
  while (1) {
    // look for updated structure status/energy;
    while (update_all_structures(&(pce->access_structure_bank()),lat,site_type_list,atom_label) || file_exists("refresh")) {
      while (file_exists("refresh")) unlink("refresh");
      // if there are changes, fit a new CE;
      CEFitInfo fitinfo;
      if (!quiet) cerr << "Finding best cluster expansion..." << endl;
      pce->find_best_cluster_choice(&fitinfo); //fit;
      write_fit_info(fitinfo,lat,site_type_list,atom_label,axes); //print results;
      if (!quiet) cerr << "done!" << endl;
    }
    // check if user wants to abort;
    if (file_exists("stop")) break;
    // check if user wants a new structure;
    if (file_exists("ready")) {
      StructureInfo *pstr;
      if (!quiet) cerr << "Finding best structure..." << endl;
      // first try to find the best structure;
      pstr=pce->find_best_structure();
      // if we don't have enough info yet, find the next structure needed to fit a minimal CE;
      if (!pstr) pstr=pce->find_initial_structures();
      // if we don't have enough info yet, just find the cheapest structure to compute;
      if (!pstr) pstr=pce->find_first_unknown_structure();
      write_structure_file(*pstr,lat,site_type_list,atom_label,axes); // write structure;
      while (file_exists("ready")) unlink("ready"); // tell user we're done;
      pstr->status=StructureInfo::busy; // mark structure as busy;
      if (!quiet) cerr << "done!" << endl;
    }
    for (int t=0; t<polltime; t++) {
      if (file_exists("refresh")) break;
      sleep(1); // wait a little;
    }
  }
  unlink("stop");
  unlink(MAPS_IS_RUNNING);
  delete pce;
}
