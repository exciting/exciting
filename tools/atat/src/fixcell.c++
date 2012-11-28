#include "lattype.h"
#include "parse.h"
#include "getvalue.h"
#include "version.h"

int main(int argc, char *argv[]) {
  int dohelp=0;
  int sigdig=6;
  int noaxes=0;
  int bravaistype=0;
  int dummy=0;
  AskStruct options[]={
    {"","fixcell " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    //    {"-h","Display more help",BOOLVAL,&dohelp},
    {"-d","Use all default values",BOOLVAL,&dummy},
    {"-c","read and print cell in cartesian only (no axes specified)",BOOLVAL,&noaxes},
    {"-z","Tolerance for finding symmetry operations (default: 1e-3)",REALVAL,&zero_tolerance},
    {"-b","print bravais lattice type only and conventional (potentially not primitive) cell",BOOLVAL,&bravaistype},

    {"-sig","Number of significant digits to print in output files",INTVAL,&sigdig}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    return 1;
  }
  if (dohelp) {
    cout << "Reads a cell from stdin writes a 'symmetrized' version of it to stdout" << endl;
    return 1;
  }
  
  // read cell;
  rMatrix3d axes,cell;
  axes.identity();
  if (!noaxes) {
    read_cell(&axes,cin);
  }
  for (int i=0; i<3; i++) {
    rVector3d v;
    cin >> v;
    cell.set_column(i,axes*v);
  }
  
  if (fabs(det(cell))<zero_tolerance) ERRORQUIT("Lattice vectors are coplanar.");

  cout.setf(ios::fixed);
  cout.precision(sigdig);

  if (bravaistype) {
    AutoString bts;
    rMatrix3d conv;
    print_bravais_type(&bts,find_bravais_type(&conv,find_reduced_cell(cell)));
    cout << bts << endl;
    cout << "Conventional cell in cartesian:" << endl;
    for (int i=0; i<3; i++) {
      cout << conv.get_column(i) << endl;
    }
    conv=(!axes)*conv;
    cout << "Conventional cell in the given coordinate system:" << endl;
    for (int i=0; i<3; i++) {
      cout << conv.get_column(i) << endl;
    }
  }
  else {
    rMatrix3d sym_cell=find_symmetric_cell(cell);
    if (!noaxes) {
      for (int i=0; i<3; i++) {
	cout << axes.get_column(i) << endl;
      }
    }
    sym_cell=(!axes)*sym_cell;
    for (int i=0; i<3; i++) {
      cout << sym_cell.get_column(i) << endl;
    }
  }
}
