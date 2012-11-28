#ifndef _CALCCORR_H_
#define _CALCCORR_H_

#include "findsym.h"
#include "gceutil.h"

void find_equivalent_clusters(Array<ArrayrVector3d> *pclusters, const Array<rVector3d> &cluster,
                              const rMatrix3d &cell, const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans);

int calc_multiplicity(const Array<rVector3d> &cluster,
                              const rMatrix3d &cell, const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans);

Real calc_correlation(const Structure &str, const Array<ArrayrVector3d> &clusters, const rMatrix3d &cell);

Real calc_correlation(const Structure &str, const Array<rVector3d> &cluster,
                      const rMatrix3d &cell, const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans);

void find_clusters_overlapping_site(Array<ArrayrVector3d> *o_cluster, Array<Real> *o_eci,
				    const rVector3d &site,
				    const Structure &lattice, const SpaceGroup &space_group,
				    const LinkedList<ArrayrVector3d> &cluster_list, const LinkedList<Real> &eci_list, int skipzeros=0);


Real find_empty_eci(const LinkedList<ArrayrVector3d> &c, const LinkedList<Real> &e);

class CorrFuncTable: public Array<Array<Array<Real> > > {
 public:
  CorrFuncTable(): Array<Array<Array<Real> > >() {}
  virtual void init(int comp) {}
  void init_from_site_type_list(const Array<Array<int> > &site_type_list);
};

void make_square_corr_matrix(Array2d<Real> *psm, const Array<Array<Real> > &rm);

class TrigoCorrFuncTable: public CorrFuncTable {
 public:
  TrigoCorrFuncTable(): CorrFuncTable() {}
  void init(int comp);
};

void find_equivalent_clusters(Array<MultiCluster> *pclusters, const MultiCluster &cluster,
                      const rMatrix3d &cell, const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans);
int calc_multiplicity(const MultiCluster &cluster,
		      const rMatrix3d &cell, const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans);
Real calc_correlation(const Structure &str, const Array<MultiCluster> &clusters, const rMatrix3d &cell, const Array<Array<Array<Real> > > &corrfunc);
Real calc_correlation(const Structure &str, const MultiCluster &cluster,
                      const rMatrix3d &cell, const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans, const Array<Array<Array<Real> > > &corrfunc);

void find_clusters_overlapping_site(Array<MultiCluster> *o_cluster, Array<int> *o_which,
				    const rVector3d &site,
				    const SpaceGroup &space_group,
				    const LinkedList<MultiCluster> &cluster_list);

void calc_corr_to_conc(Array2d<Real> *pmat, const Structure &lat, const Array<Array<int> > &site_type_list, const SpaceGroup &spacegroup, const Array<Array<Array<Real> > > &corrfunc, int usemultcorr=0);
void extract_nonredundant(Array2d<Real> *pnored, Array2d<Real> *pnoredtored, Array<Real> *pconst, const Array2d<Real> &mat);

void calc_concentration(Array<Real> *pconc, const Array<int> &atom_type, int nbatomtype);
int calc_concentration(Array<Real> *pconc, const Structure &lattice, const Array<Array<int> > &site_type_list, const Structure &str);
void calc_noredtored(Array2d<Real> *pmat, Array<Real> *pvec, const Structure &lat, const Array<Array<int> > &site_type_list);
void calc_redtonored(Array2d<Real> *pmat, const Array2d<Real> &noredtored, const Array<Real> &noredtored_c, int keepconst=0);

void find_clusters_symmetry(SpaceGroup *psubgroup, const MultiCluster &cluster, const SpaceGroup &space_group);

#endif
