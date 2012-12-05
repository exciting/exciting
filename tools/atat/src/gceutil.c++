#include "gceutil.h"
#include <fstream.h>
#include "tensorsym.h"
#include "parse.h"

GenericPlugIn<SymmetryObeyingObject> *GenericPlugIn<SymmetryObeyingObject>::list=NULL;
SpecificPlugIn<SymmetryObeyingObject,TensorSymmetryObeyingObject> TensorSymmetryObeyingObjectPlugIn("tensor");

int TensorSymmetryObeyingObject::rank=-1;
Array<Array<int> > TensorSymmetryObeyingObject::flip(0);

TensorSymmetryObeyingObject::TensorSymmetryObeyingObject(void): SymmetryObeyingObject(), tensor() {
  if (rank==-1) {
    ifstream file("gcetensor.in");
    if (!file) {ERRORQUIT("Cannot open gcetensor.in to read rank of tensor");}
    file >> rank;
    skip_delim(file);
    LinkedList<Array<int> > list;
    while (!file.eof()) {
      skip_delim(file);
      LinkedList<int> rowcollist;
      while (!is_eol(file)) {
	int row=-1,col=-1;
	file >> row >> col;
	if (row==-1 || col==-1) {ERRORQUIT("Error reading gcetensor.in");}
        rowcollist << new int(row) << new int(col);
      }
      Array<int> rowcol;
      LinkedList_to_Array(&rowcol,rowcollist);
      if (rowcol.get_size()==0) break;
      list << new Array<int>(rowcol);
    }
    LinkedList_to_Array(&flip,list);
  }
  Array<int> size(rank);
  fill_array(&size,3);
  tensor.resize(size);
  tensor.zero();
}

TensorSymmetryObeyingObject::TensorSymmetryObeyingObject(const rTensor &t): SymmetryObeyingObject(), tensor() {
  tensor=t;
}

void TensorSymmetryObeyingObject::find_basis(LinkedList<SymmetryObeyingObject> *plistobj, const SpaceGroup &space_group) {
  Array<rTensor> basis;
  find_symmetric_basis(&basis,rank,space_group,flip);
  for (int i=0; i<basis.get_size(); i++) {
    (*plistobj) << new TensorSymmetryObeyingObject(basis(i));
  }
}

void TensorSymmetryObeyingObject::zero(void) {
  tensor.zero();
}

SymmetryObeyingObject* TensorSymmetryObeyingObject::clone(void) {
  return new TensorSymmetryObeyingObject(tensor);
}

void TensorSymmetryObeyingObject::apply_symmetry(const rMatrix3d &point_op, const rVector3d &trans) {
  rTensor tmp=tensor;
  ::apply_symmetry(&tensor, point_op,tmp);
}

void TensorSymmetryObeyingObject::add(SymmetryObeyingObject *pobj, Real mult) {
  rTensor t=((TensorSymmetryObeyingObject*)pobj)->tensor;
  t*=mult;
  tensor+=t;
}

void TensorSymmetryObeyingObject::get_data(Array<Real> *pdata) {
  *pdata=tensor.vectorize();
}

void TensorSymmetryObeyingObject::write(ostream &file) {
  file << tensor.get_size();
  file << tensor.vectorize();
}

void TensorSymmetryObeyingObject::read(istream &file) {
  Array<int> size;
  file >> size;
  tensor.resize(size);
  file >> tensor.vectorize();
}

// gce begin

void find_equivalent_clusters(Array<GeneralizedCluster> *pclusters, const GeneralizedCluster &cluster,
                      const rMatrix3d &cell, const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans) {
  rMatrix3d inv_cell=!cell;
  LinkedList<GeneralizedCluster> list;
  for (int op=0; op<point_op.get_size(); op++) {
    GeneralizedCluster *ptransfo_cluster=new GeneralizedCluster();
    apply_symmetry(&(ptransfo_cluster->clus),point_op(op),trans(op),cluster.clus);
    ptransfo_cluster->psymobj=cluster.psymobj->clone();
    ptransfo_cluster->psymobj->apply_symmetry(point_op(op),trans(op));
    LinkedListIterator<GeneralizedCluster> i(list);
    for ( ; i; i++) {
      if (equivalent_mod_cell(i->clus,ptransfo_cluster->clus,inv_cell)) break;
    }
    if (!i) {
      list << ptransfo_cluster;
    }
    else {
      delete ptransfo_cluster;
    }
  }
  LinkedList_to_Array(pclusters,list);
}

int calc_multiplicity(const GeneralizedCluster &cluster,
                              const rMatrix3d &cell, const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans) {
  Array<GeneralizedCluster> equiv;
  find_equivalent_clusters(&equiv, cluster,cell,point_op,trans);
  return equiv.get_size();
}

void calc_correlation(SymmetryObeyingObject *gencorr, const Structure &str, const Array<GeneralizedCluster> &clusters, const rMatrix3d &cell, const Array<Array<Array<Real> > > &corrfunc) {
  SymmetryObeyingObject *accum=gencorr->clone();
  accum->zero();
  int count=0;
  rMatrix3d inv_strcell=!str.cell;
  LatticePointInCellIterator t(cell,str.cell);
  for ( ;t; t++) {
    for (int c=0; c<clusters.get_size(); c++) {
      Real sigma=1.;
      for (int at=0; at<clusters(c).clus.clus.get_size(); at++) {
	Real spin=corrfunc(clusters(c).clus.site_type(at))(clusters(c).clus.func(at))(str.atom_type(which_atom(str.atom_pos,t+clusters(c).clus.clus(at),inv_strcell)));
        sigma=sigma*spin;
      }
      accum->add(clusters(c).psymobj,sigma);
      count++;
    }
  }
  gencorr->zero();
  gencorr->add(accum,1./(Real)count);
}

void calc_correlation(SymmetryObeyingObject *gencorr, const Structure &str, const GeneralizedCluster &cluster,
                      const rMatrix3d &cell, const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans, const Array<Array<Array<Real> > > &corrfunc) {
  Array<GeneralizedCluster> clusters;
  find_equivalent_clusters(&clusters, cluster, cell, point_op, trans);
  calc_correlation(gencorr,str,clusters,cell,corrfunc);
}

// gce end

Real read_clusters_and_eci(LinkedList<GeneralizedCluster> *clusterlist, LinkedList<Real> *ecilist,
                           istream &clusterfile, istream &ecifile, const rMatrix3d &axes) {
   Real maxlen=0.;
   while (1) {
      Real mult;
      clusterfile >> mult;
      if (clusterfile.eof()) break;
      Real len;
      clusterfile >> len;
      int nbpt;
      clusterfile >> nbpt;
      clusterfile.get();
      MultiCluster cluster(nbpt);
      for (int j=0; j<nbpt; j++) {
        char buf[MAX_LINE_LEN];
	buf[0]=0;
	while (strlen(buf)==0) {
	  clusterfile.get(buf,MAX_LINE_LEN-1);
	  clusterfile.get();
	}
	istrstream line(buf);
        line >> cluster.clus(j);
        cluster.clus(j)=axes*(cluster.clus(j));
        line >> cluster.site_type(j) >> cluster.func(j);
	  if (!line) {
          cluster.site_type(j)=0;
          cluster.func(j)=0;
	  }
      }
      AutoString gce_label;
      get_string(&gce_label,clusterfile,"\n");
      SymmetryObeyingObject *psymobj=GenericPlugIn<SymmetryObeyingObject>::create(gce_label);
      psymobj->read(clusterfile);
      GeneralizedCluster *pcluster=new GeneralizedCluster(cluster,psymobj);
      (*clusterlist) << pcluster;
      maxlen=MAX(maxlen,get_cluster_length(cluster.clus));
      if (ecilist) {
	Real cur_eci;
	ecifile >> cur_eci;
	(*ecilist) << new Real(cur_eci);
      }
    }
  return maxlen;
}
