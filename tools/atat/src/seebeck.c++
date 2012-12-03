#include <fstream.h>
#include "parse.h"
#include "getvalue.h"
#include "version.h"
#include "lstsqr.h"
#include "xtalutil.h"

class KBand {
public:
  KBand(void): k(), band(), v() {}
  KBand(const rVector3d &_k, const Array<Real> &_band): k(_k), band(_band), v() {}
  rVector3d k;
  Array<Real> band;
  Array<rVector3d> v;
};

int main(int argc, char *argv[]) {
  char *delim="\t";
  int dohelp=0;
  char *strfilename="str.out";
  char *bandfilename="bands.out";
  int sigdig=5;
  Real kr=0.5;
  int nosym=0;
  Real fermi=0;
  Real T=0.;
  Real Emin=0.;
  Real Emax=0.;
  Real dE=0.;
  AskStruct options[]={
    {"","Seebeck coefficient calculator" MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    //    {"-h","Display more help",BOOLVAL,&dohelp},
    {"-s","Input file defining the structure (default: str.out)",STRINGVAL,&strfilename},
    {"-bf","Input file defining the band structure (default: bands.out)",STRINGVAL,&bandfilename},
    //    {"-ns","Turn off symmetry",BOOLVAL,&nosym},
    {"-z","Tolerance for finding symmetry operations (default: 1e-3)",REALVAL,&zero_tolerance},
    {"-kr","K-space radius",REALVAL,&kr},
    //    {"-ef","Fermi energy",REALVAL,&fermi},
    //    {"-T","Temperature",REALVAL,&T},
    {"-E1","Energy min",REALVAL,&Emin},
    {"-E2","Energy max",REALVAL,&Emax},
    {"-dE","Energy step",REALVAL,&dE},
    {"-sig","Number of significant digits printed (default: 5)",INTVAL,&sigdig}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    return 1;
  }
  ofstream logfile("seebeck.out");
  logfile.setf(ios::fixed);
  logfile.precision(sigdig);
  cout.setf(ios::fixed);
  cout.precision(sigdig);

  Structure str;
  Array<Arrayint> labellookup;
  Array<AutoString> label;
  rMatrix3d axes;
  {
    ifstream strfile(strfilename);
    if (!strfile) ERRORQUIT("Unable to open structure file");
    parse_lattice_file(&str.cell, &str.atom_pos, &str.atom_type, &labellookup, &label, strfile, &axes);
  }
  SpaceGroup spacegroup;
  spacegroup.cell=str.cell;
  find_spacegroup(&spacegroup.point_op,&spacegroup.trans,str.cell,str.atom_pos,str.atom_type);
  rMatrix3d rec_lat=2.*M_PI*!(~str.cell);
  rMatrix3d inv_rec_lat=!rec_lat;

  int nb_bands=0;
  ifstream bandfile(bandfilename);
  if (!bandfile) ERRORQUIT("Unable to open band file");
  bandfile >> nb_bands;
  LinkedList<KBand> E_list;
  while (1) {
    KBand Eonek;
    bandfile >> Eonek.k;
    if (bandfile.eof()) break;
    Eonek.k=flip_into_brillouin_1(rec_lat*Eonek.k,rec_lat);
    Eonek.band.resize(nb_bands);
    for (int b=0; b<nb_bands; b++) {
      bandfile >> Eonek.band(b);
    }
    for (int op=0; op<spacegroup.point_op.get_size(); op++) {
      for (Real s=-1.; s<=1.+zero_tolerance; s+=2.) {
	KBand Eonek2;
	Eonek2.k=s*spacegroup.point_op(op)*Eonek.k;
	Eonek2.band=Eonek.band;
	LinkedListIterator<KBand> i(E_list);
	for (; i; i++) {
	  if (cylinder_norm(inv_rec_lat*(Eonek2.k-(i->k)))<zero_tolerance) break;
	}
	if (!i) {
	  E_list << new KBand(Eonek2);
	}
      }
    }
  }

  rVector3d vmin(-1.,-1.,-1.);
  rVector3d vmax(1.,1.,1.);
  LinkedListIterator<KBand> i(E_list);
  for (; i; i++) {
    LinkedList<KBand> local_E_list;
    LinkedListIterator<KBand> j(E_list);
    for (; j; j++) {
      MultiDimIterator<rVector3d> img(vmin,vmax);
      for (; img; img++) {
	rVector3d k=j->k+rec_lat*(rVector3d &)img;
	if (norm(k-i->k)<kr) {
	  local_E_list << new KBand(k,j->band);
	}
      }
    }

    cerr << i->k << endl;
    i->v.resize(nb_bands);
    for (int b=0; b<nb_bands; b++) {
      int nb_nei=local_E_list.get_size();
      Array2d<Real> x(nb_nei,3);
      Array<Real> y(nb_nei);
      LinkedListIterator<KBand> j(local_E_list);
      for (int n=0; j; j++, n++) {
	y(n)=(j->band(b)) - (i->band(b));
	for (int d=0; d<3; d++) {
	  x(n,d)=(j->k(d)) - (i->k(d));
	}
      }
      cerr << x << endl;
      Array<Real> der;
      calc_ols(&der,x,y);
      {
	  Array<Real> resid,yhat;
	  predict_ols(&yhat,x,y);
	  diff(&resid,y,yhat);
	  Real nsigma2=inner_product(resid,resid);
	  Real R2=inner_product(yhat,yhat)/inner_product(y,y);
	  Real cv=calc_cv(x,y);
	  cerr << "band= " << b << " n= " << nb_nei << " R2= " << R2 << " cv= " << cv << " mse= " << sqrt(nsigma2/(Real)(y.get_size())) << endl;
      }
      for (int d=0; d<3; d++) { i->v(b)(d)=der(d); }
    }
  }

  int nE=(int)((Emax-Emin)/dE);
  Array<rMatrix3d> sigma(nE);
  Array<Real> trsigma(nE);
  Array<rMatrix3d> vel(nE);
  Array<Real> dos(nE);
  for (int n=0; n<nE; n++) {
    sigma(n).zero();
    vel(n).zero();
    trsigma(n)=0.;
    dos(n)=0.;
  }
  i.init(E_list);
  for (; i; i++) {
    for (int b=0; b<nb_bands; b++) {
      int n=(int)((i->band(b)-Emin)/dE);
      if (n>=0 && n<nE) {
	trsigma(n)+=norm2(i->v(b));
	for (int d1=0; d1<3; d1++) {
	  for (int d2=0; d2<3; d2++) {
	    sigma(n)(d1,d2)+=i->v(b)(d1)*i->v(b)(d2);
	    dos(n)+=1.;
	  }
	}
      }
    }
  }
  for (int n=0; n<nE; n++) {
    logfile << Emin+dE*(Real)n << " " << dos(n) << " " << trsigma(n) << endl;
  }
}
