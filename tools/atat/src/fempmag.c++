#include "xtalutil.h"
#include "getvalue.h"
#include "parse.h"
#include "version.h"
#include <fstream.h>

Real calc_empirical_mag_free_energy(Real spin, Real T, Real Tc, Real kB, Real p) {
  Real tau=T/Tc;
  Real A=518./1125. + (11692./15975.)*((1./p)-1.);
  Real tauf;
  if (tau<1.) {
    tauf=tau-(1./A)*(79./(140.*p) + tau*(474./497.)*((1./p)-1)*(pow(tau,3.)/6. + pow(tau,9.)/135. + pow(tau,15.)/600.) );
  }
  else {
    tauf=(-tau/A)*(pow(tau,-5.)/10. + pow(tau,-15.)/315. + pow(tau,-25.)/1500.);
  }
  return (kB*Tc*tauf*log(spin+1));
}

extern char *helpstring;

int main(int argc, char *argv[]) {
  char *strfilename="str.out";
  char *Tcurriefilename="tcurrie";
  char *magmomfilename="magmom";
  int sigdig=5;
  int dohelp=0;
  Real kboltzman=1.380658e-23/1.60217733e-19; // k_B in eV/K;
  int peratom=0;
  Real perfewer=1.;
  int refgs=0;
  Real p=0.28;
  int dummy=0;
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
    {"","Empirical magnetic free energy calculator, version " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-s","Input file defining the structure (default: str.out)",STRINGVAL,&strfilename},
    {"-tc","Input file containing the Currie Temperature (default: tcurrie)",STRINGVAL,&Tcurriefilename},
    {"-mm","Input file containing the TOTAL magnetic moment (default: magmom)",STRINGVAL,&magmomfilename},
    {"-p","Ratio of magnetic enthalpy due to SRO to total enthalpy (for fcc or hcp, use 0.28 (default), for bcc, use 0.4)",REALVAL,&p},
    {"-T0","Minimum temperature (default: 0)",REALVAL,&T0},
    {"-T1","Maximum temperature (default: 2000)",REALVAL,&T1},
    {"-dT","Temperature step (default: 100)",REALVAL,&dT},
    {"-kb","Boltzman's constant (default: in eV/K)",REALVAL,&kboltzman},
    {"-pa","Give free energy per atom",BOOLVAL,&peratom},
    {"-sc","Correction factor if spectator ion are present",REALVAL,&perfewer},
    {"-rh","Use high-temperature enthalpy as reference",BOOLVAL,&refgs},
    {"-sig","Number of significant digits to print in output files",INTVAL,&sigdig},
    {"-d","Use all default values",BOOLVAL,&dummy},
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
  Structure str;
  Array<Arrayint> labellookup;
  Array<AutoString> label;
  {
    ifstream strfile(strfilename);
    if (!strfile) ERRORQUIT("Unable to open structure file");
    parse_lattice_file(&str.cell, &str.atom_pos, &str.atom_type, &labellookup, &label, strfile);
  }
  Real Tcurrie=0.;
  {
    ifstream file(Tcurriefilename);
    if (!file) {ERRORQUIT("Unable to open currie temperature file");}
    file >> Tcurrie;
  }
  Real magmom=0.;
  Real denom=(Real)(str.atom_pos.get_size());
  {
    ifstream file(magmomfilename);
    if (!file) {ERRORQUIT("Unable to open magnetic moment file");}
    file >> magmom;
    magmom=fabs(magmom)/denom;
  }
  {
    ofstream file("fmag");
    file.setf(ios::fixed);
    file.precision(sigdig);
    for (Real T=T0; T<=T1+zero_tolerance; T+=dT) {
      Real F=calc_empirical_mag_free_energy(magmom,T,Tcurrie,kboltzman,p)*(peratom ? 1. : denom )*perfewer;
      file << F+(refgs ? -kboltzman*T*log(magmom+1) : 0) << endl;
    }
  }
}
