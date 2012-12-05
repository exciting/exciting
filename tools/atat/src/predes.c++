#include "mrefine.h"
#include "getvalue.h"
#include "parse.h"
#include "phonlib.h"
#include <fstream.h>

Real calc_ewald(const rMatrix3d &cell, const Array<rVector3d> &atom_pos, const Array<Real> &charges, Real eta=-1.0, Real prec=1e-25) {
    rMatrix3d rcell=2.*M_PI*!(~cell);
    LatticePointIterator lat(cell,1);
    Real rmin=norm((rVector3d)lat);
    LatticePointIterator rlat(rcell,1);
    Real kmin=norm((rVector3d)rlat);
    Real omega=fabs(det(cell));
    if (eta<0) {
      eta=sqrt(0.5*(kmin/rmin)*sqrt(log(prec)/log((prec*omega)/(4*M_PI))));
      //      cerr << eta << endl;
    }
    Real E=0.;
    for (int i=0; i<atom_pos.get_size(); i++) {
	for (int j=0; j<=i; j++) {
	    Real qq=charges(i)*charges(j);
	    Real doublecount=(j==i ? 0.5 : 1.);
	    rVector3d dr=atom_pos(i)-atom_pos(j);
	    lat.init(cell,(j==i));
	    Real maxr=sqrt(-log(prec))/eta+norm(dr);
	    int nt=0;
	    while (1) {
		Real r=norm((rVector3d)lat+dr);
		if (r>maxr) break;
		E+=qq*doublecount*erfc(eta*r)/r;
		lat++;
		nt++;
	    }
	    //cerr << "nReal=" << nt << endl;
	    rlat.init(rcell,1);
	    Real maxk2=4.*sqr(eta)*(-log(omega*prec/(4.*M_PI)));
	    nt=0;
	    while (1) {
		rVector3d k=(rVector3d)rlat;
		Real k2=norm2(k);
		if (k2>maxk2) break;
		E+=qq*cos(k*dr)*doublecount*exp(-k2/(4.*sqr(eta)))*4.*M_PI/omega/k2;
		rlat++;
		nt++;
	    }
	    //cerr << "nReciprocal=" << nt << endl;
	    E-=doublecount*qq*M_PI/sqr(eta)/omega;
	    if (i==j) {
		E-=qq*eta/sqrt(M_PI);
	    }
	}
    }
    return E*14.39964886; // Coulomb constant in eV*A/e^2;
}

class ElectroStaticPredictor: public EnergyPredictor {
  static Structure lattice;                            // the parent lattice;
  static Array<Arrayint> site_type_list;
  static Array<Real> species_charges;
  static Real scrn;

  int done;
  Real es_energy;
 public:
  // Constructor: don't touch;
  ElectroStaticPredictor(void);

  // init static members of class;
  void static_init(const ClusterExpansion &ce);

  // calculate electrostatic energy for structure str;
  Real get_energy(const StructureInfo &str);
};

Structure ElectroStaticPredictor::lattice;
Array<Arrayint> ElectroStaticPredictor::site_type_list;
Array<Real> ElectroStaticPredictor::species_charges;
Real ElectroStaticPredictor::scrn;

ElectroStaticPredictor::ElectroStaticPredictor(void) {done=0; es_energy=0.;}

void ElectroStaticPredictor::static_init(const ClusterExpansion &ce) {
  ElectroStaticPredictor::lattice=ce.get_lattice();
  ElectroStaticPredictor::site_type_list=ce.get_site_type_list();
  ElectroStaticPredictor::scrn=1.;
  ElectroStaticPredictor::species_charges.resize(ce.get_atom_label().get_size());
  {
    ifstream file("es.in");
    if (!file) ERRORQUIT("Unable to open es.in file.");
    while (!file.eof()) {
      AutoString specie;
      skip_delim(file," \t=\n");
      get_string(&specie,file," \t=");
      skip_delim(file," \t=");
      if (file.eof()) break;
      if (specie==AutoString("scrn")) {
	  file >> ElectroStaticPredictor::scrn;
	} else {
	  int wspecie=index_in_array(ce.get_atom_label(),specie);
        if (wspecie==-1) ERRORQUIT("Unknown specie in electrostatic predictor.");
        file >> ElectroStaticPredictor::species_charges(wspecie);
      }
    }
  }
}

Real ElectroStaticPredictor::get_energy(const StructureInfo &str) {
  if (!done) {
    rMatrix3d icell=!(lattice.cell);
    Array<Real> charges(str.atom_pos.get_size());
    Real total_charge=0.;
    int num_positive=0;
    int num_negative=0;
    for (int i=0; i<str.atom_pos.get_size(); i++) {
      int inlat=which_atom(lattice.atom_pos,str.atom_pos(i),icell);
      charges(i)=ElectroStaticPredictor::species_charges(ElectroStaticPredictor::site_type_list(ElectroStaticPredictor::lattice.atom_type(inlat))(str.atom_type(i)));
      total_charge+=charges(i);
      if (charges(i)>0) {num_positive++;}
      if (charges(i)<0) {num_negative++;}
    }
    Real d_charge=(total_charge==0. ? 0. : total_charge/(Real)(total_charge>0. ? num_positive : num_negative));
    for (int i=0; i<str.atom_pos.get_size(); i++) {
      if (sgn(charges(i))==sgn(total_charge)) {charges(i)-=d_charge;}
    }
    es_energy=calc_ewald(str.cell,str.atom_pos,charges)*(Real)(lattice.atom_pos.get_size())/(Real)(str.atom_pos.get_size());
    done=1;
  }
  return ElectroStaticPredictor::scrn*es_energy;
}

SpecificPlugIn<EnergyPredictor, ElectroStaticPredictor> ElectroStaticPredictorRegister("es"); // register plug-in under name "es";
