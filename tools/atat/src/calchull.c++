#include <fstream.h>
#include "chull.h"
#include "version.h"
#include "getvalue.h"
#include "integer.h"
#include "arraylist.h"

int main(int argc, char *argv[]) {
  char *xfilename="";
  int dummy;
  AskStruct options[]={
    {"","CALculate Convex HULL " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-f","input file",STRINGVAL,&xfilename},
    {"-z","Tolerance for finding symmetry operations (default: 1e-3)",REALVAL,&zero_tolerance},
    {"-d","Use all default values",BOOLVAL,&dummy}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    return 1;
  }

  {
    ostrstream cmd;
    cmd << "echo `cat " << xfilename << " | wc -l` "
	<< "`head -1 " << xfilename << " | wc -w` > tmp.tmp" << '\0';
    system(cmd.str());
  }
  int nb_param,nb_eqn;
  {
    ifstream size("tmp.tmp");
    size >> nb_eqn >> nb_param;
  }
  unlink("tmp.tmp");

  Array<Array<Real> > x(nb_eqn);
  ifstream xfile(xfilename);
  if (!xfile) {ERRORQUIT("Unable to open input file");}
  for (int eqn=0; eqn<nb_eqn; eqn++) {
    x(eqn).resize(nb_param);
    for (int param=0; param<nb_param; param++) {
      xfile >> x(eqn)(param);
    }
  }

  Array<Real> ground(nb_param);
  zero_array(&ground);
  ground(nb_param-1)=-1.;
  LinkedList<PolytopeFace> hull;
  calc_convex_hull(&hull, x, ground);
  {
    ofstream gs_connect("connect.out");
    LinkedListIterator<PolytopeFace> i(hull);
    for (; i; i++) {
      for (int j=0; j<i->pts.get_size(); j++) {
	gs_connect << i->pts(j) << " ";
      }
      gs_connect << endl;
    }
  }
}
