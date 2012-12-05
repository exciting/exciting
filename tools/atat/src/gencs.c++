#include <fstream.h>
#include "parse.h"
#include "clus_str.h"
#include "getvalue.h"
#include "ctype.h"
#include "version.h"
#include "plugin.h"
#include "tensorsym.h"

//extern char *helpstring;
char *helpstring="";

int main(int argc, char *argv[]) {
  char *delim="\t";
  int dohelp=0;
  char *latticefilename="lat.in";
  int minrank=1;
  int maxrank=1;
  int sigdig=5;
  AskStruct options[]={
    {"","Generalized Constituent Strain " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-h","Display more help",BOOLVAL,&dohelp},
    {"-l","Input file defining the lattice   (default: lat.in)",STRINGVAL,&latticefilename},
    {"-z","Tolerance for finding symmetry operations (default: 1e-3)",REALVAL,&zero_tolerance},
    {"-r","Maximum Rank of the harmonic",INTVAL,&maxrank},
    {"-mr","Minimum Rank of the harmonic (default 1)",INTVAL,&maxrank},
    {"-sig","Number of significant digits printed (default: 5)",INTVAL,&sigdig}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    return 1;
  }
  if (dohelp) {
    cout << helpstring;
    return 1;
  }
  Structure lattice;
  Array<Arrayint> labellookup;
  Array<AutoString> label;
  rMatrix3d axes;
  {
    ifstream latticefile(latticefilename);
    if (!latticefile) ERRORQUIT("Unable to open lattice file");
    parse_lattice_file(&lattice.cell, &lattice.atom_pos, &lattice.atom_type, &labellookup, &label, latticefile, &axes);
    wrap_inside_cell(&lattice.atom_pos,lattice.atom_pos,lattice.cell);
  }
  SpaceGroup spacegroup;
  spacegroup.cell=lattice.cell;
  find_spacegroup(&spacegroup.point_op,&spacegroup.trans,lattice.cell,lattice.atom_pos,lattice.atom_type);
  if (contains_pure_translations(spacegroup.point_op,spacegroup.trans)) {
    cerr << "Warning: unit cell is not primitive." << endl;
  }

  {
    ofstream symfile("sym.out");
    symfile.setf(ios::fixed);
    symfile.precision(sigdig);
    symfile << spacegroup.point_op.get_size() << endl;
    for (int i=0; i<spacegroup.point_op.get_size(); i++) {
      symfile << ((!axes)*(spacegroup.point_op(i))*axes) << endl;
      symfile << ((!axes)*(spacegroup.trans(i))) << endl;
      symfile << endl;
    }
  }

  ofstream clusterfile("harm.out");
  clusterfile.setf(ios::fixed);
  clusterfile.precision(sigdig);
  for (int rank=minrank; rank<=maxrank; rank++) {
    Array<rTensor> basis;
    calc_sym_harmonics(&basis,rank,spacegroup);
    for (int i=0; i<basis.get_size(); i++) {
      clusterfile << 1 << endl;
      clusterfile << 0 << endl;
      clusterfile << 0 << endl;
      clusterfile << "tensor" << endl;
      clusterfile << basis(i).get_size();
      clusterfile << basis(i).vectorize();
      clusterfile << endl;


      MultiDimIterator<Array<int> > j(basis(i).get_size());
      for (; j; j++) {
	Array<int> &aj=(Array<int> &)j;
	if (fabs(basis(i)(aj))>zero_tolerance) {
	  Array<int> p(3);
	  zero_array(&p);
	  int s=1;
	  for (int k=0; k<aj.get_size(); k++) {p(aj(k))++;}
	  for (int k=0; k<aj.get_size()-1; k++) {if (aj(k)>aj(k+1)) {s=0;}}
	  char c[]="xyz";
	  if (s) {
	    cerr << " + " << basis(i)(aj)*factorial(rank)/(factorial(p(0))*factorial(p(1))*factorial(p(2)));
	    for (int k=0; k<3; k++) {if (p(k)>0) {cerr << "*" << c[k] << "^" << p(k);}}
	  }
	}
      }
      cerr << endl;

    }
  }
}
