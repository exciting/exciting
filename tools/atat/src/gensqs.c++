#include <fstream.h>
#include "clus_str.h"
#include "parse.h"
#include "getvalue.h"
#include "version.h"
#include "calccorr.h"

extern char *helpstring;

int main(int argc, char *argv[]) {
  int maxvol=0; 
  char *latticefilename="lat.in";
  char *clusterfilename="clusters.out";
  char *corrfilename="tcorr.out";
  int readcell=0;
  int printcell=0;
  int tp=1;
  int ip=0;
  Real mysqstol=zero_tolerance;
  int sigdig=6;
  char *corrfunc_label="trigo";
  int dohelp=0;
  AskStruct options[]={
    {"","GENerate Special Quasirandom Structures " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-n","nb of atom/unit cell",INTVAL,&maxvol},
    {"-sig","Number of significant digits to print in output files",INTVAL,&sigdig},
    {"-cf","Input file defining the clusters (default: clusters.out)",STRINGVAL,&clusterfilename},
    {"-tc","Input file defining the target correlations (default: tcorr.out)",STRINGVAL,&corrfilename},
    {"-tol","Tolerance for matching correlations (default: 1e-3)",REALVAL,&mysqstol},
    {"-l","Input file defining the lattice (default: lat.in)",STRINGVAL,&latticefilename},
    {"-rc","Read unit cells from file",BOOLVAL,&readcell},
    {"-pc","Print number of supercells only, and quit",BOOLVAL,&printcell},
    {"-tp","Total number of processes (for parallel operation)",INTVAL,&tp},
    {"-ip","Index of current process (0,...,tp-1) (for parallel operation)",INTVAL,&ip},
    {"-crf","Select correlation functions (default: trigo)",STRINGVAL,&corrfunc_label},
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

  // read in lattice (see parse.h);
  Structure lat;
  Array<Arrayint> site_type_list;
  Array<AutoString> atom_label;
  rMatrix3d axes;
  ifstream file(latticefilename);
  if (!file) ERRORQUIT("Unable to open lattice file.");
  parse_lattice_file(&lat.cell, &lat.atom_pos, &lat.atom_type, &site_type_list, &atom_label, file, &axes);

  Array<int> comp(lat.atom_type.get_size());
  for (int at=0; at<comp.get_size(); at++) {
    comp(at)=site_type_list(lat.atom_type(at)).get_size();
  }

  // initialize a table of correlation functions;
  if (!check_plug_in(CorrFuncTable(),corrfunc_label)) {
    ERRORQUIT("Aborting");
  }
  CorrFuncTable *pcorrfunc=GenericPlugIn<CorrFuncTable>::create(corrfunc_label);
  pcorrfunc->init(max(comp));

  LinkedList<MultiCluster> clusterlist;
  LinkedList<Real> corrlist;
  {
    ifstream clusterfile(clusterfilename);
    if (!clusterfile) ERRORQUIT("Unable to open cluster file.");
    ifstream corrfile(corrfilename);
    if (!corrfile) ERRORQUIT("Unable to open target correlation file.");
    read_clusters_and_eci(&clusterlist, &corrlist, clusterfile, corrfile, axes);
  }
  Array<Real> corrarray;
  LinkedList_to_Array(&corrarray,corrlist);

  cout.setf(ios::fixed);
  cout.precision(sigdig);

  Array<rMatrix3d> pointgroup;
  find_pointgroup(&pointgroup,lat.cell);
  SpaceGroup spacegroup;
  find_spacegroup(&spacegroup.point_op,&spacegroup.trans,lat.cell,lat.atom_pos,lat.atom_type);
  spacegroup.cell=lat.cell;

  Array<rMatrix3d> supercell;
  int v=(int)(maxvol/MAX(1,lat.atom_pos.get_size()));
  if (readcell) {
    ifstream cellfile("sqscell.out");
    int nc=0;
    cellfile >> nc;
    supercell.resize(nc);
    if (printcell) {
      cout << supercell.get_size() << endl;
      exit(0);
    }
    for (int i=0; i<supercell.get_size(); i++) {
      read_cell(&(supercell(i)),cellfile);
    }
  }
  else {
    find_supercells(&supercell, v, v, lat.cell, pointgroup);
    if (printcell) {
      cout << supercell.get_size() << endl;
      exit(0);
    }
    ofstream cellfile("sqscell.out");
    cellfile.setf(ios::fixed);
    cellfile.precision(sigdig);
    cellfile << supercell.get_size() << endl << endl;
    for (int i=0; i<supercell.get_size(); i++) {
      supercell(i)=find_symmetric_cell(supercell(i));
      write_axes(supercell(i),cellfile,0);
      cellfile << endl;
    }
  }
  if (ip<0) {
    exit(0);
  }
  ofstream outfile;
  if (tp>1) {
    ostrstream filename;
    filename << "gensqs_" << ip << "_" << tp << ".out" << '\0';
    outfile.open(filename.str());
  }

  LinkedList<Structure> strlist;
  int totcell=supercell.get_size();
  int bcell=ip*totcell/tp;
  int ecell=(ip+1)*totcell/tp;
  for (int c=bcell; c<ecell; c++) {
    {
      ostrstream filename;
      filename << "gensqs_" << ip << "_" << tp << ".stat" << '\0';
      ofstream file(filename.str());
      file << bcell << " " << c << " " << ecell << endl;
    }
    Structure blank_superstructure;
    blank_superstructure.cell=supercell(c);
    find_all_atom_in_supercell(&blank_superstructure.atom_pos,
		 &blank_superstructure.atom_type,lat.atom_pos,
		 lat.atom_type,
		 lat.cell, blank_superstructure.cell);
    for (int i=0; i<blank_superstructure.atom_type.get_size(); i++) {
      blank_superstructure.atom_type(i)=site_type_list(blank_superstructure.atom_type(i)).get_size();
    }
    Array<Array<MultiCluster> > eqclus(clusterlist.get_size());
    LinkedListIterator<MultiCluster> ic(clusterlist);
    for (int t=0; ic; t++, ic++) {
      find_equivalent_clusters(&eqclus(t),*ic,lat.cell,spacegroup.point_op,spacegroup.trans);
    }
    MultiDimIterator<Arrayint> config(blank_superstructure.atom_type);
    for (; config; config++) {
      blank_superstructure.atom_type=config;
      int t;
      for (t=0; t<corrarray.get_size(); t++) {
        if (fabs(corrarray(t)-calc_correlation(blank_superstructure,eqclus(t),lat.cell,*pcorrfunc))  > mysqstol ) break;
      }
      if (t==corrarray.get_size()) {
	if (!contains_pure_translations_or_lexico_successor(blank_superstructure,lat.cell)) {
	  LinkedListIterator<Structure> is(strlist);
	  for (; is; is++) {
	    if (spacegroup(*is,blank_superstructure)) break;
	  }
	  if (!is) {
	    strlist << new Structure(blank_superstructure);
	    if (tp>1) {
	      write_structure(blank_superstructure,lat,site_type_list,atom_label,axes,outfile);
	      outfile << "end" << endl << endl << flush;
	    }
	    else {
	      write_structure(blank_superstructure,lat,site_type_list,atom_label,axes,cout);
	      cout << "end" << endl << endl << flush;
	    }
	  }  
	}
      }
    }
  }
  if (tp>1) {
    outfile.close();
  }
}
