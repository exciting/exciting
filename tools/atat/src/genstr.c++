#include <fstream.h>
#include "clus_str.h"
#include "parse.h"
#include "getvalue.h"
#include "version.h"
#include "mpiinterf.h"

int main(int argc, char *argv[]) {
  MyMPIobj.init(argc,argv);
  int maxvol=0; 
  char *latticefilename="lat.in";
  int sigdig=6;
  int do2D=0;
  AskStruct options[]={
    {"","GENerate STRuctures " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-n","maximum nb of atom/unit cell",INTVAL,&maxvol},
    {"-sig","Number of significant digits to print in output files",INTVAL,&sigdig},
    {"-2d","Find supercells along a and b axes only",BOOLVAL,&do2D},
    {"-l","Input file defining the lattice (default: lat.in)",STRINGVAL,&latticefilename}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
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

  int noconfig=1;
  for (int i=0; i<lat.atom_type.get_size(); i++) {
    if (site_type_list(lat.atom_type(i)).get_size()>1) noconfig=0;
  }

  cout.setf(ios::fixed);
  cout.precision(sigdig);
  if (noconfig) {
    Array<int> config(lat.atom_type.get_size());
    zero_array(&config);
    Array<rMatrix3d> pointgroup;
    if (lat.atom_pos.get_size()>0) {
      SpaceGroup spacegroup;
      spacegroup.cell=lat.cell;
      find_spacegroup(&spacegroup.point_op,&spacegroup.trans,lat.cell,lat.atom_pos,lat.atom_type);
      pointgroup_from_spacegroup(&pointgroup,spacegroup.point_op);
    }
    else {
      find_pointgroup(&pointgroup,lat.cell);
    }
    Array<rMatrix3d> supercell;
    if (do2D) {
      find_supercells_2D(&supercell, 1, maxvol/MAX(1,lat.atom_pos.get_size()), lat.cell, pointgroup);
    }
    else {
      find_supercells(&supercell, 1, maxvol/MAX(1,lat.atom_pos.get_size()), lat.cell, pointgroup);
    }
    for (int c=0; c<supercell.get_size(); c++) {
      Structure blank_superstructure;
      blank_superstructure.cell=find_symmetric_cell(supercell(c));
      find_all_atom_in_supercell(&blank_superstructure.atom_pos,
				 &blank_superstructure.atom_type,lat.atom_pos,
				 config,
				 lat.cell, blank_superstructure.cell);
      write_structure(blank_superstructure,lat,site_type_list,atom_label,axes,cout);
      cout << "end" << endl << endl;
    }
    
  }
  else {
    SpaceGroup spacegroup;
    spacegroup.cell=lat.cell;
    find_spacegroup(&spacegroup.point_op,&spacegroup.trans,lat.cell,lat.atom_pos,lat.atom_type);
    StructureBank<Structure> str_bank(lat, site_type_list, spacegroup,do2D);
    while (1) {
	if (MyMPIobj.is_root()) {
	    write_structure(str_bank.get_current_structure(),lat,site_type_list,atom_label,axes,cout);
	    cout << "end" << endl << endl;
	}
	if (str_bank.get_current_structure().atom_pos.get_size()>=maxvol && str_bank.end_of_this_volume()) break;
	str_bank++;
    }
  }
}
