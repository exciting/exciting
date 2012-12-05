#include "parse.h"
#include "getvalue.h"
#include "version.h"

int main(int argc, char *argv[]) {
  int quiet=0;
  int all_mat=0;
  int dummy=0;
  AskStruct options[]={
    {"","check cell distortion " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-d","default",BOOLVAL,&dummy},
    {"-q","be quiet",BOOLVAL,&quiet},
    {"-p","print strain",BOOLVAL,&all_mat}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    return 1;
  }
  rMatrix3d before,after;
  if (!quiet) cout << "Enter original cell: (in cartesian, one lattice vector per line with coordinates separated by spaces)" << endl;
  read_cell(&before,cin);
  if (!quiet) cout << "Enter relaxed cell: (in cartesian, one lattice vector per line with coordinates seperated by spaces)" << endl;
  read_cell(&after,cin);
  if (all_mat) {
    if (!quiet) cout << "Cell distortion" << endl;
    cout << (!before)*after;
  }
  else {
    before=before/pow(det(before),1./3.);
    after=after/pow(det(after),1./3.);
    rMatrix3d diff=(!before)*after;
    rMatrix3d id;
    id.identity();
    diff=(diff+~diff)/2-id;
    cout.setf(ios::fixed);
    cout.precision(4);
    if (!quiet) cout << "Distorsion measure:" << endl;
    cout << norm(diff) << endl;
  }
}
