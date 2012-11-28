#include <fstream.h>
#include "xtalutil.h"
#include "kmeci.h"
#include "linalg.h"
#include "plugin.h"
#include "calccorr.h"

class ES_KSpaceECI: public KSpaceECI {
  Structure lattice;
  iVector3d supercell;
  Real scrn;
  Array<Array<int> > site_type_list;
  Array<Real> species_charges;
  Array<Array<Array<Real> > > corrfunc;
 public:
  // Constructor: don't touch;
  ES_KSpaceECI(void): KSpaceECI(), lattice(), supercell(), site_type_list(), species_charges() {}
  // read in parameters from the file es.in;
  void init(const Structure &_lattice, const Array<Array<int> > &_site_type_list, const Array<AutoString> &atom_label, const iVector3d &_supercell, const Array<Array<Array<Real> > > &_corrfunc) {
    lattice=_lattice;
    site_type_list=_site_type_list;
    supercell=_supercell;
    corrfunc=_corrfunc;

    species_charges.resize(atom_label.get_size());
    zero_array(&species_charges);
    // open file and error checking;
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
	    file >> scrn;
	  } else {
	    int wspecie=index_in_array(atom_label,specie);
          if (wspecie==-1) ERRORQUIT("Unknown specie in electrostatic predictor.");
          file >> species_charges(wspecie);
        }
      }
    }
  }
  Real FT_trunc_coulomb(Real k, Real c) {
    if (k==0.) {return 0.;}
    Real kc=k*c;
    Real k2=k*k;
    Real coskc=cos(kc);
    Real sinkc=sin(kc);
    return (-15.*coskc*pow(kc,-2.) + 75.*sinkc*pow(kc,-3.) +180.*coskc*pow(kc,-4.) -180.*sinkc*pow(kc,-5.))*4.*M_PI/k2;
  }
  void get_k_space_eci(Array<Array<Array<Array<Array<Complex> > > > > *p_ft_eci, const Array<Real> &x) {
    Real kelec=14.39964886; // Coulomb constant in eV*A/e^2;
    Array<Array<Real> > corr_to_chg(lattice.atom_pos.get_size());
    for (int s=0; s<lattice.atom_pos.get_size(); s++) {
      Array2d<Real> corrmat,icorrmat;
      int t=lattice.atom_type(s);
      int comp=site_type_list(t).get_size();
      make_square_corr_matrix(&corrmat,corrfunc(comp-2));
      invert_matrix(&icorrmat,corrmat);
      Array<Real> chg(comp);
      for (int ss=0; ss<comp; ss++) {
        chg(ss)=species_charges(site_type_list(t)(ss));
      }
      product(&corr_to_chg(s), chg,icorrmat);
    }

    rMatrix3d mat_supercell;
    mat_supercell.diag(to_real(supercell));
    rMatrix3d recip_lat=!(~(lattice.cell));
    rMatrix3d k_mesh=!(~(lattice.cell*mat_supercell));
    rMatrix3d k_big_mesh=!(~(lattice.cell));
    Real omega=det(lattice.cell);

    AtomPairIterator nn(lattice.cell,lattice.atom_pos);
    Real c=nn.length();
    iVector3d maxh(6,6,6);

    int readdump=0;
    ofstream outdumpfile;
    ifstream indumpfile("esdump.in");
    if (indumpfile) {
      readdump=1;
      cerr << "Reading esdump.in" << endl;
    }
    else {
      outdumpfile.open("esdump.out");
    }

    for (int s=0; s<lattice.atom_pos.get_size(); s++) {
      for (int t=0; t<=s; t++) {
        for (int ss=0; ss<corr_to_chg(s).get_size(); ss++) {
          for (int tt=0; tt<corr_to_chg(t).get_size(); tt++) {
            MultiDimIterator<iVector3d> cell(supercell);
            iVector3d &vcell=(iVector3d &)cell;
            for ( ; cell; cell++) {
              //int offset=((vcell(0)*supercell(1) + vcell(1))*supercell(2) + vcell(2)) //permbug;
              int offset=((vcell(2)*supercell(1) + vcell(1))*supercell(0) + vcell(0));
	      //rVector3d k0=flip_into_brillouin_1(k_mesh*to_real(cell),recip_lat);
	      rVector3d k0=k_mesh*to_real(cell);
	      Complex sum=0;
	      if (readdump) {
		indumpfile >> sum;
		if (!indumpfile) {
		  readdump=0;
		  cerr << "Stopped reading esdump.in and started writing esdump.out" << endl;
		  outdumpfile.open("esdump.out", ios::app);
		}
	      }
	      if (!readdump) {
		MultiDimIterator<iVector3d> h(-maxh,maxh);
		for ( ; h; h++) {
		  rVector3d k=2.*M_PI*(k0+k_big_mesh*to_real(h));
		  Real phase=(k*(lattice.atom_pos(t)-lattice.atom_pos(s)));
		  Complex exp_phase(cos(phase),sin(phase));
		  Complex deci=kelec*scrn*exp_phase*FT_trunc_coulomb(norm(k),c)*corr_to_chg(s)(ss)*corr_to_chg(t)(tt)/omega;
//cerr << ((iVector3d &)h) << " " << deci << endl;
		  sum+=deci/2.;
		}
		outdumpfile << sum << endl;
	      }
	      (*p_ft_eci)(s)(t)(ss)(tt)(offset)=sum;
//cerr << "end" << endl;
            }
          }
        }
      }
    }
  }
};

SpecificPlugIn<KSpaceECI,ES_KSpaceECI> ES_KSpaceECIRegister("es"); // register plug-in under name "es";
