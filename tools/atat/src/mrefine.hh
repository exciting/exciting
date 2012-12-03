#ifndef _RAFFINE_H_
#define _RAFFINE_H_

#include "clus_str.h"
#include "linalg.h"
#include "chull.h"
#include "lstsqr.h"
#include "plugin.h"

typedef LinkedList<Real> LinkedListReal;

class StructureInfo;
class ClusterExpansion;

class EnergyPredictor {
 public:
  EnergyPredictor(void) {}
  virtual void static_init(const ClusterExpansion &ce) {}
  virtual Real get_energy(const StructureInfo &str) {return 0.;}
};

class StructureInfo: public Structure {
 public:
  enum Status {unknown=1,calculated=2,busy=4,error=8,anystatus=15};
  Real energy;
  Status status;
  Real cost;
  AutoString label;
  LinkedList<EnergyPredictor> predictor;
  LinkedList<LinkedListReal> correlations;
  StructureInfo(void): Structure(),
    label(), correlations(), predictor()
    {
      energy=0;
      status=unknown;
      cost=0.;
    }
  StructureInfo(const Structure &str): Structure(str),
    label(), correlations(), predictor()
    {
      energy=0;
      status=unknown;
      cost=0.;
    }
};

Real calc_structure_cost(const Structure &str, const rMatrix3d &cell, const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans, Real complexity_exp);

// Describes the status of a cluster expansion fit;
class CEFitInfo {
 public:
  enum Status {fit_ok=0, fit_impossible=1, gs_problem=2, new_gs=4};
  // Respectively means:
  //  no problem;
  //  no enough data to fit;
  //  predicted ground states do not agree with currently known ground states;
  //  new ground states (whose exact energy is not known yet) have been predicted;
  Status status;
  Array<StructureInfo *> pstr; // Pointers to structures included in the fit;
  Array<Real> pure_energy; // Energies of the end members (e.g. pure elements);
  Array<Array<Real> > concentration; // Concentrations     \;
  Array<Real> fitted_energy; // Predicted energies | of the str used in the fit;
  Array<Real> energy;        // True energies      /;

  Array<Real> weight; // weights needed to get the right GS line;
  Real cv;            // crossvalidation score (when weights turned off);
  Array<int> gs_list; // indices (into any of the 4 previous arrays) of the true ground states
  LinkedList<PolytopeFace> hull;
  Array2d<Real> corr_to_fullconc; // converts correlation to redundant concentration;
  Array2d<Real> corr_to_conc; // converts correlation to non redundant concentration;
  Array2d<Real> conc_to_fullconc; // converts nonredundant concentration to full redundant concentration;
  Array<Real> conc_to_fullconc_c; // constant terms of previous matrix;
  int max_volume; // largest volume considered when looking for new ground states;
  Array<MultiCluster> cluster; // clusters used in the optimal fit;
  Array<Real> multiplicity; // corresponding multiplicities;
  Array<Real> eci; // corresponding ECI;
};

class ClusterExpansion {
 protected:
  SpaceGroup spacegroup; // see findsym.h;
  StructureBank<StructureInfo> structures; // structure generation object     \ see clus_str.h;
  Array<MultiClusterBank *> pclusters; // cluster generation object /;
                  // *(pclusters(0)): points, *(pclusters(1)): pairs, etc.;
  Array<LinkedList<Array<MultiCluster> > * > equiv_cluster; // image of each cluster by all elements of the factor group;
  Array<Real> weight; // weights used in the optimal fit;
  Array<Real> eci_mult; // fitted ECI (times multiplicity);
  Structure parent_lattice; // the lattice;
  Structure lattice_only;   // the lattice, with fixed atoms removed;
  Array<Arrayint> site_type_list; // as read by parse_lattice_file;
  Array<AutoString> atom_label; //  as read by parse_lattice_file;
  int minimal_nc; // number of cluster in a minimal nearest neighbor cluster expansion;
  int nb_calculated_str; // number of calculated structures;
  int nb_clusters; // number of clusters used in the fit;
  Real radius_1nn; // radius of the nearest neighbor shell;

  CorrFuncTable *pcorrfunc;
  Array2d<Real> corr_to_fullconc; // converts correlation to redundant concentration;
  Array2d<Real> corr_to_conc; // converts correlation to non redundant concentration;
  Array2d<Real> conc_to_fullconc; // converts nonredundant concentration to full redundant concentration;
  Array<Real> conc_to_fullconc_c; // constant terms of previous matrix;
  Array<Real> down_e;
  LinkedList<LinearInequality> conc_range; // range of concentration sampled;
  int ignore_gs;

  LinkedList<int> label_list;
  LinkedListIterator<int> i_label_list;
  int cur_label;
  const char *predictor_labels;
  Real high_energy;
 public:
  int max_multiplet; // maximum number of points in cluster (can be increased later on);
  Real complexity_exp;
  ClusterExpansion(const Structure &_parent_lattice,
                   const Array<Arrayint> &_site_type_list,
                   const Array<AutoString> &_atom_label,
                   const SpaceGroup &_spacegroup,
		   CorrFuncTable *_pcorrfunc);
  ~ClusterExpansion(void);

  // Low-level functions;
  // update corr. of one structure;
  void update_correlations(StructureInfo *str);
  // update corr. of all structures of a given status (multiple statuses can be or-ed together);
  void update_correlations(StructureInfo::Status select);
  // update list of equivalent clusters;
  void update_equivalent_clusters(void);
  // initialize the static members of the user-specified energy predictors;
  void init_predictors(void);
  // calculate the matrix that maps point correlation into concentrations;
  void calc_point_to_concentration(Array<Real> *multiplicity, int nb_cols);
  // calculate concentration from correlation vector;
  void calc_concentration(Array<Real> *pconcentration, const Array<Real> &corr_matrix);
  // calculate concentration from correlation matrix;
  void calc_concentration(Array<Array<Real> > *pconcentration, const Array2d<Real> &corr_matrix);
  // create a matrix of corr. and a vector of energies for structures maching a given status;
  void calc_regression_matrices( Array2d<Real> *pcorr_matrix, Array<Real> *penergy, StructureInfo::Status select=StructureInfo::calculated);
  // create a vector of predictor energies (using the user-specified predictor objects) for structures maching a given status;
  void calc_predictor_energy( Array<Real> *penergy, StructureInfo::Status select=StructureInfo::calculated);
  Real calc_predictor_energy(StructureInfo *str);
  // restart from minimal cluster expansion (only 1nn pairs); 
  void reset_cluster_choice(void);
  // iterate to the next choice of clusters;
  int find_next_cluster_choice(void);
  // select a given choice of clusters: choice(1) gives the nunber of pairs, choice(2) gives the number of triplets, etc.;
  void set_cluster_choice(const Array<int> &choice);
  // opposite of the above;
  void get_cluster_choice(Array<int> *pchoice);
  // make an array of correlations from the list of corr. stored in a StructureInfo object;
  void make_correlation_array(Array<Real> *array_corr, const LinkedList<LinkedListReal> &list_corr);
  // restart from first structure;
  void reset_structure(void);
  // iterate to the next structure;
  void find_next_structure(void);

  // High-level functions;
  // Set the concentration range where the ground states must be right;
  void set_concentration_range(const LinkedList<LinearInequality> &_conc_range, int _ignore_gs);
  // select predictors;
  void set_predictor_labels(const char *labels);
  // set max energy above GS hull allowed when generating new structures;
  void set_highest_energy_allowed(Real _high_energy) {high_energy=_high_energy;};

  // the following 4 routines can be overridden in derived class to enable new fitting algorithms;

  // Main routine calling the cluster enumeration and crossvalidation routines;
  // Initialize *pfitinfo with results (see class CEFitInfo above);
  virtual void find_best_cluster_choice(CEFitInfo *pfitinfo);
  // Main routine calling the structure enumeration routines and the variance calculation codes.;
  // Return pointer to structure found (or NULL upon failure);
  virtual StructureInfo * find_best_structure(void);
  // Finds the first few structures needed to initiate the refinement process.;
  // Return pointer to structure found (or NULL upon failure);
  virtual StructureInfo * find_initial_structures(void);
  // Returns the least computationally intensive structure still marked as unknown.
  // Return pointer to structure (never NULL);
  virtual StructureInfo * find_first_unknown_structure(void);

  // various functions to access private members;
  const Structure & get_lattice(void) const {
    return parent_lattice;
  }
  const Structure & get_lattice_only(void) const {
    return lattice_only;
  }
  const Array<Arrayint> & get_site_type_list(void) const {
    return site_type_list;
  }
  const Array<AutoString> & get_atom_label(void) const {
    return atom_label;
  }
  const LinkedList<StructureInfo> & get_structure_list(void) {
    return structures.get_structure_list();
  }
  LinkedList<StructureInfo> & access_structure_list(void) {
    return structures.get_structure_list();
  }
  StructureBank<StructureInfo> & access_structure_bank(void) {
    return structures;
  }
  const LinkedList<MultiCluster> & get_cluster_list(int ntuple) const {
    return pclusters(ntuple)->get_cluster_list();
  }
};

// Setup plugins to allow for alternate cluster expansion algorithms;

class ClusterExpansionCreator { // for the default algorithm;
 public:
  virtual ClusterExpansion* create(const Structure &_parent_lattice, const Array<Arrayint> &_site_type_list, const Array<AutoString> &_atom_label, const SpaceGroup &_spacegroup, CorrFuncTable *_pcorrfunc) {
    return new ClusterExpansion(_parent_lattice,_site_type_list,_atom_label,_spacegroup,_pcorrfunc);
  }
};

template<class T>
class CustomClusterExpansionCreator: public ClusterExpansionCreator { // for alternate algorithm;
  virtual ClusterExpansion* create(const Structure &_parent_lattice, const Array<Arrayint> &_site_type_list, const Array<AutoString> &_atom_label, const SpaceGroup &_spacegroup, CorrFuncTable *_pcorrfunc) {
    return new T(_parent_lattice,_site_type_list,_atom_label,_spacegroup,_pcorrfunc);
  }
};

// end setting up plugins;

#endif
