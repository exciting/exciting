#include <fstream.h>
#include "xtalutil.h"
#include "getvalue.h"
#include "parse.h"
#include "version.h"


void list_atom(ostream &file, const Array<Arrayint> &labellookup, const Array<AutoString> &label, int thetype) {
  for (int i=0; i<labellookup(thetype).get_size(); i++) {
    if (i>0) {cout << ",";}
    cout << label(labellookup(thetype)(i));
  }
}


int main(int argc, char *argv[]) {
  char *latfilename="lat.in";
  int listbond=0;
  Real user_rad=0.0;
  int sigdig=5;
  int dummy;
  AskStruct options[]={
    {"","Nearest Neighbor Shell " MAPS_VERSION ", by Axel van de Walle and Paul Dalach",TITLEVAL,NULL},
    {"-l","Input file defining the lattice (default: lat.in)",STRINGVAL,&latfilename},
    {"-lb","List bond types in given structures",BOOLVAL,&listbond},
    {"-sig","Number of significant digits to print in output files",INTVAL,&sigdig},
    {"-r","Define a max radius for neighbor search",REALVAL,&user_rad},
    {"-d","Use all default values",BOOLVAL,&dummy}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    return 1;
  }
  Structure lat;
  Array<Arrayint> labellookup;
  Array<AutoString> label;
  rMatrix3d axes;
  {
    ifstream latfile(latfilename);
    if (!latfile) ERRORQUIT("Unable to open lattice file");
    parse_lattice_file(&lat.cell, &lat.atom_pos, &lat.atom_type, &labellookup, &label, latfile, &axes);
    wrap_inside_cell(&lat.atom_pos,lat.atom_pos,lat.cell);
  }
  cout.setf(ios::fixed);
  cout.precision(sigdig);
  if (listbond) {
    for (int at=0; at<lat.atom_type.get_size(); at++) {
      Real rnn;
      if (user_rad>0) {rnn=user_rad;}
      else {rnn=find_1nn_radius(lat,at);}
      int body=3;
      Array<Real> ang(body);
      AtomMultipletIterator clus(lat.cell, lat.atom_pos, body);
      for ( ; clus.length()<=user_rad+zero_tolerance; clus++) {
	for (int i=0; i<body; i++) {
	  rVector3d r1=clus((i+1)%body)-clus(i%body);
	  rVector3d r2=clus((i+2)%body)-clus(i%body);
	  ang(i)=acos(r1*r2/sqrt((r1*r1)*(r2*r2)));
	}
	cout << max(ang)*180/M_PI << " ";
	/*
	for (int i=0; i<body; i++) {
	  rVector3d r1=clus((i+1)%body)-clus(i%body);
	  cout << norm(r1) << " ";
	}
	*/
	int base=index_max(ang);
	Array<int> thetype(body);
	for (int i=0; i<body; i++) {
	  thetype(i)=lat.atom_type(which_atom(lat.atom_pos,clus((base+i)%body),!lat.cell));
	}
	list_atom(cout,labellookup,label,thetype(0));
	cout << " ";
	list_atom(cout,labellookup,label,MIN(thetype(1),thetype(2)));
	cout << " ";
	list_atom(cout,labellookup,label,MAX(thetype(1),thetype(2)));
	cout << endl;
      }
    }
  }
  else {
    for (int at=0; at<lat.atom_type.get_size(); at++) {
      int nbnn=0;
      cout << (!axes)*lat.atom_pos(at) << endl;
      Real rnn;
      if(user_rad>0) {rnn=user_rad;}
      else {rnn=find_1nn_radius(lat,at);}
      cout << rnn << endl;
      list_atom(cout,labellookup,label,lat.atom_type(at));
      cout << endl;
      AtomPairIterator p(lat.cell,lat.atom_pos(at),lat.atom_pos);
      while (p.length()<=rnn+zero_tolerance) {
	cout << "   " << (!axes)*(p(1)-p(0)) << " " << norm(p(1)-p(0)) << " ";
	int thetype=lat.atom_type(which_atom(lat.atom_pos,p(1),!lat.cell));
	list_atom(cout,labellookup,label,thetype);
	cout << endl;
	nbnn++;
	p++;
      }
      cout << "  " << nbnn << endl << endl;
    }
  }
}
