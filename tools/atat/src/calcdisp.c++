#include "getvalue.h"
#include "parse.h"
#include "phonlib.h"
#include "version.h"
#include <fstream.h>

rVector3d cylindrize(const rVector3d &d, const rMatrix3d &cell) {
    rVector3d half(0.5,0.5,0.5);
    rVector3d wd=cell*(mod1((!cell)*d+half)-half);
    return wd;
}

int main(int argc, char *argv[]) {
  char *unrelfilename="str.out";
  char *relfilename="str_relax.out";
  int defval=0;
  AskStruct options[]={
    {"","CALCulate DISPlacements " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-u","Unrelaxed structure (default: str.out)",STRINGVAL,&unrelfilename},
    {"-r","Relaxed structure (default: str_relax.out)",STRINGVAL,&relfilename},
    {"-d","Use all defaults",BOOLVAL,&defval}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    return 1;
  }
  Structure ustr,rstr;
  Array<Arrayint> site_type_list;
  Array<AutoString> atom_label;
  rMatrix3d axes;
  {
    ifstream file(unrelfilename);
    if (!file) ERRORQUIT("Unable to open unrelaxed file.");
    parse_lattice_file(&ustr.cell, &ustr.atom_pos, &ustr.atom_type, &site_type_list, &atom_label, file, &axes);
  }
  {
    ifstream file(relfilename);
    if (!file) ERRORQUIT("Unable to open relaxed file.");
    parse_lattice_file(&rstr.cell, &rstr.atom_pos, &rstr.atom_type, &site_type_list, &atom_label, file, &axes);
  }
  reorder_atoms(&rstr, ustr);
  for (int i=0; i<ustr.atom_pos.get_size(); i++) {
      rVector3d d=(rstr.atom_pos(i)-ustr.atom_pos(i));
      d=cylindrize(d,ustr.cell);
      cout << atom_label(site_type_list(ustr.atom_type(i))(0)) << " " << ustr.atom_pos(i) << " " << d << endl;
  }
}
