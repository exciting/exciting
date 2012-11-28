#include <fstream.h>
#include <sys/stat.h>
#include "phonlib.h"
#include "getvalue.h"
#include "parse.h"
#include "version.h"
#include "lattype.h"
#include "arraylist.h"

void find_distinct_sites(Array<int> *pwhich_atom,
		       const Structure &str,
		       const SpaceGroup &spacegroup) {
  LinkedList<int> which_atom_list;
  LinkedList<Array<rVector3d> > pointlist;
  for (int at=0; at<str.atom_pos.get_size(); at++) {
    Array<rVector3d> p(1);
    p(0)=str.atom_pos(at);
    if (add_unique(&pointlist,p,spacegroup)) {
      which_atom_list << new int(at);
    }
  }
  LinkedList_to_Array(pwhich_atom,which_atom_list);
}

// extern char *helpstring;
char *helpstring="";

int main(int argc, char *argv[]) {
  cerr.setf(ios::fixed | ios::showpos);
  cerr.precision(5);
  char *latfilename="lat.in";
  char *strfilename="str.out";
  char *strprefix="pdef";
  Real enclosed_radius=0.;
  int dohelp=0;
  int sigdig=5;
  AskStruct options[]={
    {"","Point DEFect generator " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-l","Input file defining the lattice (default: lat.in)",STRINGVAL,&latfilename},
    {"-s","Input file defining the lattice (default: str.out)",STRINGVAL,&strfilename},
    {"-p","Prefix of the directories that will contain the defected structures (Default pdef)",STRINGVAL,&strprefix},
    {"-er","Minimum distance between point defects",REALVAL,&enclosed_radius},
    {"-sig","Number of significant digits printed (default: 5)",INTVAL,&sigdig},
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

  if (enclosed_radius==0) {
    ERRORQUIT("You need to specify at least the -er parameter.");
  }

  Structure lat;
  Array<Arrayint> site_type_list;
  Array<AutoString> atom_label;
  rMatrix3d axes;
  {
    ifstream file(latfilename);
    if (!file) ERRORQUIT("Unable to open lattice file.");
    parse_lattice_file(&lat.cell, &lat.atom_pos, &lat.atom_type, &site_type_list, &atom_label, file, &axes);
    if (fabs(det(lat.cell))<zero_tolerance) ERRORQUIT("Lattice vectors are coplanar.");
    wrap_inside_cell(&lat.atom_pos,lat.atom_pos,lat.cell);
  }
  
  Structure str;
  {
    ifstream file(strfilename);
    if (!file) ERRORQUIT("Unable to open structure file.");
    parse_structure_file(&str.cell,&str.atom_pos,&str.atom_type,atom_label,file);
    wrap_inside_cell(&str.atom_pos,str.atom_pos,str.cell);
    fix_atom_type(&str, lat,site_type_list,0);
  }

  Array<int> atoms;
  SpaceGroup str_spacegroup;
  str_spacegroup.cell=str.cell;
  find_spacegroup(&str_spacegroup.point_op,&str_spacegroup.trans,str.cell,str.atom_pos,str.atom_type);
  find_distinct_sites(&atoms,str,str_spacegroup);

  Structure ideal_supercell;
  find_smallest_supercell_enclosing_sphere(&(ideal_supercell.cell),str.cell,enclosed_radius);
  find_all_atom_in_supercell(&(ideal_supercell.atom_pos), &(ideal_supercell.atom_type), str.atom_pos, str.atom_type, str.cell, ideal_supercell.cell);
  
  int pert=0;
  for (int at=0; at<atoms.get_size(); at++) {
    int i_in_lat=which_atom(lat.atom_pos,str.atom_pos(atoms(at)),!(lat.cell));
    int i_in_supstr=which_atom(ideal_supercell.atom_pos,str.atom_pos(atoms(at)),!(ideal_supercell.cell));
    for (int s=0; s<site_type_list(lat.atom_type(i_in_lat)).get_size(); s++) {
      Structure supstr(ideal_supercell);
      if (s!=ideal_supercell.atom_type(i_in_supstr)) {
	supstr.atom_type(i_in_supstr)=s;
	ostrstream line;
	line << strprefix << "_" << pert << '\0';
	mkdir(line.str(),S_IRWXU | S_IRWXG | S_IRWXO);
	chdir_robust(line.str());
	{
	  ofstream file("str.out");
	  file.setf(ios::fixed);
	  file.precision(sigdig);
	  if (!file) ERRORQUIT("Unable to open output file.");
	  write_structure(supstr,lat,site_type_list,atom_label,axes,file);
	  ofstream waitfile("wait");
	}
	chdir_robust("..");
	pert++;
      }
    }
  }
}
