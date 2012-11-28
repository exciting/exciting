#include <fstream.h>
#include <ctype.h>
#include "parse.h"
#include "getvalue.h"
#include "version.h"
#include "stringo.h"

void to_lower(char *buf) {
  char *p;
  for (p=buf; (*p)!=0; p++) {
    *p=tolower(*p);
  }
}

int main(int argc, char *argv[]) {
  int maxvol=0; 
  char *latticefilename="lat.in";
  char *radiifilename="rad.in";
  int sigdig=6;
  int dummy=0;
  Real scale=1.;
  AskStruct options[]={
    {"","nntouch " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {""," Scales a lattice so that nearest neighbor atom just touch.",TITLEVAL,NULL},
    {"-l","Input file defining the lattice (default: lat.in)",STRINGVAL,&latticefilename},
    {"-r","Input file defining the atomic radii (default: rad.in)",STRINGVAL,&radiifilename},
    {"-s","Multiplicative scaling factor",REALVAL,&scale},
    {"-sig","Number of significant digits to print in output files",INTVAL,&sigdig},
    {"-d","Use all default values",BOOLVAL,&dummy}
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
  {
    ifstream file(latticefilename);
    if (!file) ERRORQUIT("Unable to open lattice file.");
    parse_lattice_file(&lat.cell, &lat.atom_pos, &lat.atom_type, &site_type_list, &atom_label, file, &axes);
  }

  Array<AutoString> radii_label;
  Array<Real> radii;
  {
    LinkedList<AutoString> radii_label_list;
    LinkedList<Real> radii_list;
    ifstream file(radiifilename);
    if (!file) ERRORQUIT("Unable to open atomic radii file.");
    while (!file.eof()) {
      char buf[MAX_LINE_LEN];
      buf[0]=0;
      file >> buf;
      if (strlen(buf)==0) break;
      to_lower(buf);
      radii_label_list << new AutoString(buf);
      Real r;
      file >> r;
      radii_list << new Real(r);
    }
    LinkedList_to_Array(&radii_label,radii_label_list);
    LinkedList_to_Array(&radii,radii_list);
  }

  Real max_r=0.;
  Array<Real> site_r(lat.atom_pos.get_size());
  for (int i=0; i<site_r.get_size(); i++) {
    site_r(i)=0.;
    for (int j=0; j<site_type_list(lat.atom_type(i)).get_size(); j++) {
      Real r=0.;
      int l;
      for (l=0; l<radii.get_size(); l++) {
	char tmp[MAX_LINE_LEN];
	strcpy(tmp,atom_label(site_type_list(lat.atom_type(i))(j)));
	to_lower(tmp);
	for (int i=0; i<strlen(tmp); i++) {
	  if (!isalpha(tmp[i])) {
	    tmp[i]=0;
	    break;
	  }
	}
        if (strcmp(radii_label(l),tmp)==0) {
          r=radii(l);
          break;
        }
      }
      if (l==radii.get_size()) {
        cerr << "Don't know radius of atom " << atom_label(site_type_list(lat.atom_type(i))(j)) << endl;
        ERRORQUIT("Exiting.");
      }
      site_r(i)=MAX(site_r(i),r);
      max_r=MAX(max_r,r);
    }
  }

  Real max_scale=0.;
  rMatrix3d inv_cell=!(lat.cell);
  AtomPairIterator pair(lat.cell,lat.atom_pos);
  while (max_scale*pair.length()<=2.*max_r) {
    Real sum_r=0.;
    for (int i=0; i<2; i++) {
      sum_r+=site_r(which_atom(lat.atom_pos,pair(i),inv_cell));
    }
    Real s=sum_r/pair.length();
    max_scale=MAX(max_scale,s);
    pair++;
  }
  max_scale*=scale;

  axes=axes*max_scale;
  lat.cell=lat.cell*max_scale;
  for (int i=0; i<lat.atom_pos.get_size(); i++) {
    lat.atom_pos(i)=lat.atom_pos(i)*max_scale;
  }

  cout.setf(ios::fixed);
  cout.precision(sigdig);
  write_lattice(lat,site_type_list,atom_label,axes,cout);
}
