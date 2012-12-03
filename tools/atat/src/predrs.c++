#include "predrs.h"

RSPredictorBase::RSPredictorBase(void): EnergyPredictor() {
  done=0;
}

// init static members of class (those shared among all objects of type RSPredictor;
void RSPredictorBase::static_init(const ClusterExpansion &ce) {
  // copy parent lattice from the cluster expansion object (excluding spectator atoms);
  lattice=ce.get_lattice_only();
  // set atom_type(i) to the number of possible species that can sit at site i (1 or 2);
  for (int i=0; i<lattice.atom_type.get_size(); i++) {
    lattice.atom_type(i)=ce.get_site_type_list()(lattice.atom_type(i)).get_size();
  }
  get_kspace_eci_object()->static_init(lattice);
}

// calculate constituent strain energy for structure str;
Real RSPredictorBase::get_energy(const StructureInfo &str) {
  if (!done) {
    LinkedList<rVector3d> kpts_list; // will contain a growing list of kpoints;
    // calc all lattices needed: ~=transpose, !=inverse;
    rMatrix3d rec_lat=~(!lattice.cell); // reciprocal lattice of parent lattice;
    rMatrix3d i_rec_lat=!rec_lat;       // inverse matrix of reciprocal of parent lattice;
    rMatrix3d i_lat=!lattice.cell;      // inverse of parent lattice;
    rMatrix3d rec_lat_str=~(!str.cell); // reciprocal lattice of unit cell of the structure;
    int ik=0; // include the gamma point;
    // int ik=1; // exclude the gamma point;
    LatticePointIterator k(rec_lat_str,ik); // k will scan all superlattice points of rec_lat_str;
    vol=(int)rint(det(rec_lat)/det(rec_lat_str)); // volume in multiple of unit cell of parent lattice;
    // find all (vol-1) k-points closest to the origin that are not equivalent by a translation (i.e. the 1st Brillouin zone);
    while (ik<vol) {
      LinkedListIterator<rVector3d> old_k(kpts_list); // scan all previously chosen k-points;
      for (; old_k; old_k++) {                        // in this loop;
	if (equivalent_mod_cell(*old_k,rVector3d(k),i_rec_lat)) break; // if we already have that k-point exit loop;
      }
      if (!old_k) { // if we don;t already have that k-point, store it;
	kpts_list << new rVector3d(k);
	ik++;
      }
      k++; // go to next k-point;
    }
    LinkedList_to_Array(&kpts,kpts_list); // convert the list of k-point into an array of k-points;
    // cerr << "k=" << kpts << endl;
    // compute structure factors;
    str_fact.resize(kpts.get_size());
    for (int i=0; i<kpts.get_size(); i++) { // loop over k-points;
      str_fact(i).resize(lattice.atom_pos.get_size());
      zero_array(&(str_fact(i)));
      for (int j=0; j<str.atom_pos.get_size(); j++) { // loop over atoms in unit cell of structure;
	int w_atom=which_atom(lattice.atom_pos,str.atom_pos(j),i_lat); // find index of atom within unit cell of lattice;
	if (w_atom!=-1) { // skip spectator atoms;
	  Real sigma=(str.atom_type(j)==1 ? 1. : -1.);
	  Real phase=2.*M_PI*kpts(i)*str.atom_pos(j); // convention: atom position (and not cell corner) determines phase shift;
	  str_fact(i)(w_atom)+=sigma*Complex(cos(phase),sin(phase));
	}
      }
      for (int w_atom=0; w_atom<lattice.atom_pos.get_size(); w_atom++) {
	str_fact(i)(w_atom)/=((Real)vol); // normalize so that structure factor is independent of cell size;
      }
    }
    // calc concentration;
    conc=0.;
    for (int j=0; j<str.atom_pos.get_size(); j++) {
      if (which_atom(lattice.atom_pos,str.atom_pos(j),i_lat)!=-1) {  // skip spectator atoms;
	conc+=(str.atom_type(j)==1 ? 1. : 0.);
      }
    }
    conc/=(Real)vol;
    done=1; // we won't need to do this again;
  }
  // compute energy;
  Real e=0.;
  for (int i=0; i<kpts.get_size(); i++) { // loop over k-points;
    Array2d<Complex> keci;
    get_kspace_eci_object()->get_k_space_eci(&keci,kpts(i),conc);
    e+=quadratic_form(keci,str_fact(i));
  }
  return e;
}

Structure RSPredictorBase::lattice;           // the parent lattice;

GenericPlugIn<KSpaceECI> *GenericPlugIn<KSpaceECI>::list=NULL;   // plugin list for kspace eci objects;

