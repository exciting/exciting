#include <fstream.h>
#include "xtalutil.h"
#include "getvalue.h"
#include "parse.h"
#include "version.h"

Real calc_avg_spin(Real maxspin, Real x) {
  if (x>=MAXFLOAT) {
    return(maxspin);
  }
  Real num=0.;
  Real den=0.;
  for (Real spin=-maxspin; spin<=maxspin+zero_tolerance; spin+=1.) {
    num+=spin*exp(spin*x);
    den+=exp(spin*x);
  }
  return (num/den);
}

Real find_avg_spin(Real maxspin, Real x) {
  if (x>=MAXFLOAT) {
    return(maxspin);
  }
  Real below=0.;
  Real above=maxspin;
  while ((above-below)>1e-5) {
    Real avg=(above+below)/2.;
    if (calc_avg_spin(maxspin,avg*x)<avg) {
      above=avg;
    }
    else {
      below=avg;
    }
  }
  return ((above+below)/2.);
}

Real calc_mag_free_energy(Real maxspin, Real coupling, Real T) {
  Real E0=-coupling*maxspin*maxspin;
  if (T<zero_tolerance) {
    return 0.;
  }
  else {
    Real avg=find_avg_spin(maxspin,coupling/T);
    Real num=0.;
    for (Real spin=-maxspin; spin<=maxspin+zero_tolerance; spin+=1.) {
      num+=exp(coupling*spin*avg/T);
    }
    return -T*log(num)-E0;
  }
}

int main(int argc, char *argv[]) {
  char *strfilename="str.out";
  char *magatom="";
  Real coupling=0.;
  Real maxspin=0.5;
  int sigdig=5;
  int dummy;
  Real kboltzman=1.380658e-23/1.60217733e-19; // k_B in eV/K;
  Real T0=0.;
  Real T1=2000.;
  Real dT=100.;
  if (file_exists("../Trange.in")) {
    ifstream tfile("../Trange.in");
    cerr << "Reading Trange.in file." << endl;
    T0=0;
    int nT=0;
    tfile >> T1 >> nT;
    dT=(T1-T0)/MAX(1,nT-1);
  }

  AskStruct options[]={
    {"","Free energy - MAGnetic contribution " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-s","Input file defining the structure (default: str.out)",STRINGVAL,&strfilename},
    {"-ma","Magnetic atom",STRINGVAL,&magatom},
    {"-cc","Coupling constant",REALVAL,&coupling},
    {"-sp","Spin",REALVAL,&maxspin},
    {"-T0","Minimum temperature (default: 0)",REALVAL,&T0},
    {"-T1","Maximum temperature (default: 2000)",REALVAL,&T1},
    {"-dT","Temperature step (default: 100)",REALVAL,&dT},
    {"-kb","Boltzman's constant (default: in eV/K)",REALVAL,&kboltzman},
    {"-sig","Number of significant digits to print in output files",INTVAL,&sigdig},
    {"-d","Use all default values",BOOLVAL,&dummy}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    return 1;
  }
  Structure str;
  Array<Arrayint> labellookup;
  Array<AutoString> label;
  {
    ifstream strfile(strfilename);
    if (!strfile) ERRORQUIT("Unable to open structure file");
    parse_lattice_file(&str.cell, &str.atom_pos, &str.atom_type, &labellookup, &label, strfile);
    wrap_inside_cell(&str.atom_pos,str.atom_pos,str.cell);
  }
  int magindex=-1;
  for (int at=0; at<labellookup.get_size(); at++) {
    if (labellookup(at).get_size()!=1) {ERRORQUIT("Not a structure file. Aborting.");}
    if (strcmp(label(labellookup(at)(0)),magatom)==0) {
      magindex=at;
    }
  }
  Real nbmag=0.;
  for (int at=0; at<str.atom_type.get_size(); at++) {
    if (str.atom_type(at)==magindex) nbmag+=1.;
  }
  Real nbmagpair=0.;
  Real nbpair=0.;
  Real rnn=find_1nn_radius(str)+zero_tolerance;
  rMatrix3d invcell=!str.cell;
  AtomPairIterator p(str.cell,str.atom_pos);
  for ( ; p.length()<rnn; p++) {
    if (str.atom_type(which_atom(str.atom_pos,p(0),invcell))==magindex && str.atom_type(which_atom(str.atom_pos,p(1),invcell))==magindex) {nbmagpair+=1.;}
    nbpair+=1.;
  }

  cout.setf(ios::fixed);
  cout.precision(sigdig);
  for (Real T=T0; T<=T1+zero_tolerance; T+=dT) {
    Real f=nbmag*calc_mag_free_energy(maxspin,coupling*(nbmagpair/nbpair),kboltzman*T);
    Real s=nbmag*find_avg_spin(maxspin,coupling*(nbmagpair/nbpair)/(kboltzman*T));
    cout << T << '\t' << f << '\t' << s << endl;
  }
}
