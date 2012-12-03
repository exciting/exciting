#include <fstream.h>
#include "getvalue.h"
#include "integer.h"
#include "version.h"
#include "arraylist.h"
#include "parse.h"

Real calc_nelec(const Array<Real> &dos, Real mu, Real w) {
  Real ne=0.;
  for (int i=0; i<dos.get_size(); i++) {
    ne+=dos(i)/(exp(((Real)i-mu)/w)+1.);
  }
  return(ne);
}

Real calc_Selec(const Array<Real> &dos, Real mu, Real w) {
  Real S=0.;
  for (int i=0; i<dos.get_size(); i++) {
    Real f=1./(exp(((Real)i-mu)/w)+1.);
    if (f>1e-5 && (1-f)>1e-5) {
      S+=-dos(i)*(f*log(f)+(1.-f)*log(1.-f));
    }
  }
  return(S);
}

Real calc_Eelec(const Array<Real> &dos, Real mu, Real w, Real dE) {
  Real E=0.;
  for (int i=0; i<dos.get_size(); i++) {
    Real e=dE*(Real)i;
    Real f;
    if ((e-mu) < -30.*w) {
      f=1.;
    }
    else if ((e-mu) > 30.*w) {
      f=0.;
    }
    else {
      f=1./(exp((e-mu)/w)+1.);
    }
    E+=dos(i)*e*f;
  }
  return(E);
}

int file_exists(char *filename) {
  ifstream file(filename);
  if (file) {return 1;} else {return 0;}
}

extern char *helpstring;

int main(int argc, char *argv[]) {
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

  char *dosfilename="dos.out";
  Real kboltzman=1.380658e-23/1.60217733e-19; // k_B in eV/K;
  int peratom=0;
  Real perfewer=1.;
  Real smooth=0.;
  int sigdig=5;
  int nofreee=0;
  int dohelp=0;
  int dummy=0;
  AskStruct options[]={
    {"","Electronic free energy calculator " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-dos","DOS input file name (default: dos.out)",STRINGVAL,&dosfilename},
    {"-T0","Minimum temperature (default: 0)",REALVAL,&T0},
    {"-T1","Maximum temperature (default: 2000)",REALVAL,&T1},
    {"-dT","Temperature step (default: 100)",REALVAL,&dT},
    {"-kb","Boltzman's constant (default: in eV/K)",REALVAL,&kboltzman},
    {"-pa","Output free energy per atom instead of per unit cell",BOOLVAL,&peratom},
    {"-sc","Correction factor if spectator ion are present",REALVAL,&perfewer},
    {"-sd","Smooth DOS with a Gaussian this width",REALVAL,&smooth},
    {"-nf","Do not calculate free energy",BOOLVAL,&nofreee},
    {"-sig","Number of significant digits to print in output files",INTVAL,&sigdig},
    {"-h","Display more help",BOOLVAL,&dohelp},
    {"-d","Use all default values",BOOLVAL,&dummy}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    return 1;
  }
  if (dohelp) {
    cout << helpstring;
    return 1;
  }

  Real denom_atom=1.;
  if (peratom) {
    ifstream file("str.out");
    if (!file) ERRORQUIT("Unable to open str.out, needed for -pa option");
    Structure str;
    Array<Arrayint> labellookup;
    Array<AutoString> label;
    parse_lattice_file(&str.cell, &str.atom_pos, &str.atom_type, &labellookup, &label, file);
    denom_atom=(Real)(str.atom_pos.get_size());
  }

  Real nelec=0;
  Real dE=0.;
  Real scale=0.;
  Array<Real> dos;
  {
    ifstream dosfile(dosfilename);
    if (!dosfile) {ERRORQUIT("Unable to open DOS input file.");}
    dosfile >> nelec >> dE >> scale;
    LinkedList<Real> list;
    while (skip_delim(dosfile)) {
      Real dens;
      dosfile >> dens;
      list << new Real(dens*dE);
    }
    LinkedList_to_Array(&dos,list);
  }

  Real ne=0.;
  int i=0;
  while (ne<nelec && i<dos.get_size()) {
    ne+=dos(i);
    i++;
  }
  Real Ef=dE*(i-1);
  {
    Array<Real> sdos(dos.get_size());
    if (smooth>0.) {
      zero_array(&sdos);
      Real se=smooth/dE;
      Real c=1./sqrt(2*M_PI)/se;
      for (int i=0; i<dos.get_size(); i++) {
	for (int j=0; j<dos.get_size(); j++) {
	  sdos(i)+=dos(j)*c*exp(-sqr((Real)(i-j)/se));
	}
      }
    } else {
      sdos=dos;
    }
    ofstream dosout("plotdos.out");
    for (int i=0; i<sdos.get_size(); i++) {
      dosout << dE*(Real)i-Ef << '\t' << sdos(i)/dE << endl;
    }
  }

  if (nofreee) {return 0;}

  Real ne_tol=1e-7*nelec;

  Real E0=MAXFLOAT;
  {
    ofstream logfile("felec.log");
    ofstream outfile("felec");
    for (Real T=T0; T<=T1+zero_tolerance; T+=dT) {
      if (T<zero_tolerance) {
	logfile << 0 << '\t' << 0 << endl;
	outfile << 0 << endl;
      }
      else {
	Real El=Ef;
	Real Eh=Ef;
	while (calc_nelec(dos,El/dE,kboltzman*T/dE)>nelec) El-=kboltzman*T;
	while (calc_nelec(dos,Eh/dE,kboltzman*T/dE)<nelec) Eh+=kboltzman*T;
	Real ne;
	Real mu;
	do {
	  mu=(Eh+El)/2.;
	  ne=calc_nelec(dos,mu/dE,kboltzman*T/dE);
	  if (ne<nelec) El=mu; else Eh=mu;
	} while (fabs(ne-nelec)>ne_tol);
	Real f=-kboltzman*T*calc_Selec(dos,mu/dE,kboltzman*T/dE);
	if (T>0.) {
	  if (E0==MAXFLOAT) {
	    E0=calc_Eelec(dos,mu,kboltzman*dT,dE);
	  }
	  f+=calc_Eelec(dos,mu,kboltzman*T,dE)-E0;
	}
	f*=scale*perfewer/denom_atom;
	logfile << T << '\t' << f << endl;
	outfile << f << endl;
      }
    }
  }
}
