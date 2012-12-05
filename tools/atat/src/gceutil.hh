#ifndef __GCEUTIL_H__
#define __GCEUTIL_H__

#include "plugin.h"
#include "findsym.h"
#include "tensor.h"

class SymmetryObeyingObject {
public:
  SymmetryObeyingObject(void) {}
  virtual void find_basis(LinkedList<SymmetryObeyingObject> *plistobj, const SpaceGroup &space_group) {}
  virtual void zero(void) {}
  virtual SymmetryObeyingObject* clone(void) {return NULL;}
  virtual void apply_symmetry(const rMatrix3d &point_op, const rVector3d &trans) {}
  virtual void add(SymmetryObeyingObject *pobj, Real mult) {}
  virtual void get_data(Array<Real> *pdata) {}
  virtual void write(ostream &file) {}
  virtual void read(istream &file) {}
};

class GeneralizedCluster {
public:
  MultiCluster clus;
  SymmetryObeyingObject *psymobj;
public:
  GeneralizedCluster(void): clus() {psymobj=NULL;}
  GeneralizedCluster(const GeneralizedCluster &gc): clus(gc.clus) {if (!gc.psymobj) {psymobj=NULL;} else {psymobj=gc.psymobj->clone();}}
  void operator=(const GeneralizedCluster &gc) {clus=gc.clus; if (!gc.psymobj) {psymobj=NULL;} else {psymobj=gc.psymobj->clone();}}
  ~GeneralizedCluster() {delete psymobj;}
  GeneralizedCluster(const MultiCluster &_clus, SymmetryObeyingObject *_psymobj): clus(_clus) {psymobj=_psymobj;}
};

class TensorSymmetryObeyingObject: public SymmetryObeyingObject {
public:
  static int rank;
  static Array<Array<int> > flip;
  rTensor tensor;
public:
  TensorSymmetryObeyingObject(void);
  TensorSymmetryObeyingObject(const rTensor &t);
  void find_basis(LinkedList<SymmetryObeyingObject> *plistobj, const SpaceGroup &space_group);
  void zero(void);
  SymmetryObeyingObject* clone(void);
  void apply_symmetry(const rMatrix3d &point_op, const rVector3d &trans);
  void add(SymmetryObeyingObject *pobj, Real mult=1);
  void get_data(Array<Real> *pdata);
  void write(ostream &file);
  void read(istream &file);
};

void find_equivalent_clusters(Array<GeneralizedCluster> *pclusters, const GeneralizedCluster &cluster,
			      const rMatrix3d &cell, const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans);
int calc_multiplicity(const GeneralizedCluster &cluster,
		      const rMatrix3d &cell, const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans);
void calc_correlation(SymmetryObeyingObject *gencorr, const Structure &str, const Array<GeneralizedCluster> &clusters, const rMatrix3d &cell, const Array<Array<Array<Real> > > &corrfunc);
void calc_correlation(SymmetryObeyingObject *gencorr, const Structure &str, const GeneralizedCluster &cluster,
                      const rMatrix3d &cell, const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans, const Array<Array<Array<Real> > > &corrfunc);

Real read_clusters_and_eci(LinkedList<GeneralizedCluster> *clusterlist, LinkedList<Real> *ecilist,
                           istream &clusterfile, istream &ecifile, const rMatrix3d &axes);

#endif
