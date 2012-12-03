#include "xtalutil.h"
#include "clus_str.h"
#include "linalg.h"
#include "lstsqr.h"
#include "multipoly.h"

class PairSpring {
 public:
  int whichatom[2];
  rVector3d dr;
  rVector3d cell_shift;
  rMatrix3d forcek;
  PairSpring(void): dr(0.,0.,0.), forcek() {whichatom[0]=0; whichatom[1]=0; forcek.zero();}
  PairSpring(const PairSpring &ps): dr(ps.dr), cell_shift(ps.cell_shift), forcek(ps.forcek) {whichatom[0]=ps.whichatom[0]; whichatom[1]=ps.whichatom[1];}
};

class SpringSvsL {
 public:
  int atom_type[2];
  Array<Array<Real> > forcek;
};

class ExtendableSpringSvsL {
 public:
  ExtendableSpringSvsL(void) {pfunc=NULL;}
  ~ExtendableSpringSvsL() {delete pfunc;}
  int atom_type[2];
  ArrayFunctionArray<Real> *pfunc;
};

void remove_vacancies(Structure *p_str, const Array<AutoString> &label, const char *specie="Vac");

void add_rVector3d_to_ArrayReal(Array<Real> *pa, int at, const rVector3d &v);

void make_nn_forcek(LinkedList<PairSpring> *p_pair, const Structure &str_ideal, const Structure &str_relax, const LinkedList<SpringSvsL> &spring_data);

void make_nn_forcek(LinkedList<PairSpring> *p_pair, const Structure &str_ideal, const Structure &str_relax, const LinkedList<SpringSvsL> &spring_data, const LinkedList<rMatrix3d> &dirdep_mat, const Array<Real> &conc, const MultiDimPoly &multipoly);
void make_nn_forcek(LinkedList<PairSpring> *p_pair, const Structure &str_ideal, const Structure &str_relax, const LinkedList<SpringSvsL> &spring_data, const LinkedList<rMatrix3d> &dirdep_mat, const Array<Real> &conc, const MultiDimPoly &multipoly, LinkedList<Array<rVector3d> > *pold_pair);

void make_nn_forcek(LinkedList<PairSpring> *p_pair, const Structure &str_ideal, const Structure &str_relax, const LinkedList<ExtendableSpringSvsL> &spring_data, const LinkedList<rMatrix3d> &dirdep_mat, const Array<Real> &conc);

int get_basis_size(const LinkedList<PairSpring> &pair);

void calc_dynmat(Array2d<Complex> *p_dynmat, const rVector3d &k, const LinkedList<PairSpring> &pair, Array<Real> mass);

int list_phonon_freq(LinkedList<Real> *p_freq, LinkedList<Real> *p_weight, const iVector3d &kmesh, const rVector3d &kshift, const LinkedList<PairSpring> &pair, const rMatrix3d &cell, const Array<Real> &mass, Real convfact);


void calc_normal_modes(Array<Real> *p_freq, Array2d<Complex> *p_eigvect, const rVector3d &k, const LinkedList<PairSpring> &pair, const Array<Real> &mass, Real convfact);

void calc_dispersion_curve(LinkedList<Array<Real> > *p_freq, LinkedList<Array2d<Complex> > *p_eigvect, const rVector3d &k1, const rVector3d &k2, int nstep, const LinkedList<PairSpring> &pair, const Array<Real> &mass, Real convfact);

Real real_erf(Real x);
Real ln_real_erf(Real x);

Real calc_vib_free_energy_robust(LinkedList<Real> freq, LinkedList<Real> weight, Real kBT, Real hplanck, Real length);
Real calc_vib_free_energy(LinkedList<Real> freq, LinkedList<Real> weight, Real kBT, Real hplanck);

Real calc_vib_entropy(LinkedList<Real> freq, LinkedList<Real> weight, Real kBT=MAXFLOAT, Real hplanck=1.);

void calc_k_mesh(iVector3d *p_mesh, const rMatrix3d &cell, Real nbkpts);

void reorder_atoms(Structure *rel_str, const Structure &str, Array<int> *pcopy_from=NULL, Array<iVector3d> *pshift=NULL);

Real calc_bulk_modulus(const LinkedList<PairSpring> &pair, Real volume);

void find_special_direction(Array<rVector3d> *pspecial_dir, Array<int> *pminus_equiv,
			    const Array<rMatrix3d> &point_op);

void find_pertubations(Array<int> *pwhich_atom, Array<rVector3d> *pdirection, Array<int> *pdo_minus,
		       const Structure &str,
		       const SpaceGroup &spacegroup);

void stretch_str(Structure *pstr, const Structure &str, Real strain);
void calc_pert_str(Structure *ppert_str, const rMatrix3d &supercell, const Structure &str, const rVector3d &pos, const rVector3d &displ, Real strain);

void read_force_vector(Array<Real> *pforce, istream &file, const Array<int> &copyfrom);

Real smooth_kernel(Real x);
void smooth_density(Real *xmin, Real *xmax, Array<Real> *pf, const Array<Real> &x);
void smooth_density(Real *xmin, Real *xmax, Array<Real> *pf, const Array<Real> &x, const Array<Real> &w);

void find_sym_const_tensor(LinkedList<rMatrix3d> *pfc_list, const Array<rMatrix3d> &pointgroup, int forcesym);
