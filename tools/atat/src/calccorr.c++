#include "calccorr.h"
#include "clus_str.h"
#include "lstsqr.h"
#include "plugin.h"

//#include "parse.h"

void find_equivalent_clusters(Array<ArrayrVector3d> *pclusters, const Array<rVector3d> &cluster,
                      const rMatrix3d &cell, const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans) {
  rMatrix3d inv_cell=!cell;
  LinkedList<ArrayrVector3d> list;
  for (int op=0; op<point_op.get_size(); op++) {
    Array<rVector3d> *ptransfo_cluster=new Array<rVector3d>(cluster.get_size());
    apply_symmetry(ptransfo_cluster,point_op(op),trans(op),cluster);
    LinkedListIterator<ArrayrVector3d> i(list);
    for ( ; i; i++) {
      if (equivalent_mod_cell(*i,*ptransfo_cluster,inv_cell)) break;
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

int calc_multiplicity(const Array<rVector3d> &cluster,
                              const rMatrix3d &cell, const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans) {
  Array<ArrayrVector3d> equiv;
  find_equivalent_clusters(&equiv, cluster,cell,point_op,trans);
  return equiv.get_size();
}

Real calc_correlation(const Structure &str, const Array<ArrayrVector3d> &clusters, const rMatrix3d &cell) {
  Real accum=0.;
  int count=0;
  rMatrix3d inv_strcell=!str.cell;
  LatticePointInCellIterator t(cell,str.cell);
  for ( ;t; t++) {
    for (int c=0; c<clusters.get_size(); c++) {
      Real sigma=1.;
      for (int at=0; at<clusters(c).get_size(); at++) {
	int spin=str.atom_type(which_atom(str.atom_pos,t+clusters(c)(at),inv_strcell));
        sigma=sigma*(2*spin-1);
      }
      accum+=sigma;
      count++;
    }
  }
  return accum/(Real)count;
}

Real calc_correlation(const Structure &str, const Array<rVector3d> &cluster,
                      const rMatrix3d &cell, const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans) {
  Array<ArrayrVector3d> clusters;
  find_equivalent_clusters(&clusters, cluster, cell, point_op, trans);
  return calc_correlation(str,clusters,cell);
}

void find_clusters_overlapping_site(Array<ArrayrVector3d> *o_cluster, Array<Real> *o_eci,
                  const rVector3d &site,
                  const Structure &lattice, const SpaceGroup &space_group,
                  const LinkedList<ArrayrVector3d> &cluster_list, const LinkedList<Real> &eci_list, int skipzeros) {
  rMatrix3d inv_cell=!lattice.cell;
  LinkedList<ArrayrVector3d> long_cluster_list;
  LinkedList<Real> long_eci_list;
  LinkedListIterator<ArrayrVector3d> c(cluster_list);
  LinkedListIterator<Real> e(eci_list);
  for ( ; c; c++,e++) {
    if (c->get_size()>0) {
      if (!(skipzeros && (*e)==0)) {
	Array<ArrayrVector3d> clusters;
	find_equivalent_clusters(&clusters, *c,space_group.cell,space_group.point_op,space_group.trans);
	for (int ec=0; ec<clusters.get_size(); ec++) {
	  for (int center_atom=0; center_atom<clusters(ec).get_size(); center_atom++) {
	    rVector3d lat_shift=site-clusters(ec)(center_atom);
	    if (is_int(inv_cell*lat_shift)) {
	      ArrayrVector3d *to_add=new ArrayrVector3d(clusters(ec).get_size());
	      int j=0;
	      (*to_add)(j)=site;
	      j++;
	      for (int i=0; i<clusters(ec).get_size(); i++) {
		if (i!=center_atom) {
		  (*to_add)(j)=clusters(ec)(i)+lat_shift;
		  j=j+1;
		}
	      }
	      long_cluster_list << to_add;
	      long_eci_list << new Real(*e);
	    }
	  }
	}
      }
    }  
  }
  LinkedList_to_Array(o_cluster,long_cluster_list);
  LinkedList_to_Array(o_eci,long_eci_list);
}

/*
void calc_correlation(Array<Real> *correlation, int components, const Structure &str, const Array<ArrayrVector3d> &clusters, const rMatrix3d &cell) {
  Real accum=0.;
  int count=0;
  LatticePointInCellIterator t(cell,str.cell);
  for ( ;t; t++) {
    for (int c=0; c<clusters.get_size(); c++) {
      Real sigma=1.;
      for (int at=0; at<clusters(c).get_size(); at++) {
        sigma=sigma*str.atom_type(which_atom(t+clusters(c)(at),str.atom_pos));
      }
      accum+=sigma;
      count++;
    }
  }
  return accum/(Real)count;
}
*/

Real find_empty_eci(const LinkedList<ArrayrVector3d> &cluster_list, const LinkedList<Real> &eci_list) {
  Real E_ref=0.;
  LinkedListIterator<ArrayrVector3d> ic(cluster_list);
  LinkedListIterator<Real> ie(eci_list);
  for (;ic; ic++, ie++) {
    if (ic->get_size()==0) E_ref=*ie;
  }
  return E_ref;
}

void CorrFuncTable::init_from_site_type_list(const Array<Array<int> > &site_type_list) {
  int c=0;
  for (int i=0; i<site_type_list.get_size(); i++) {
    c=MAX(c,site_type_list(i).get_size());
  }
  init(c);
}

void make_square_corr_matrix(Array2d<Real> *psm, const Array<Array<Real> > &rm) {
  int comp=rm(0).get_size();
  psm->resize(iVector2d(comp,comp));
  for (int i=0; i<comp; i++) {
    (*psm)(0,i)=1.;
    for (int j=1; j<comp; j++) {
      (*psm)(j,i)=rm(j-1)(i);
    }
  }
}

void TrigoCorrFuncTable::init(int comp) {
  Array<Array<Array<Real> > > &table=*this;
  table.resize(comp-1);
  for (int m=2; m<=comp; m++) {
    table(m-2).resize(m-1);
    for (int t=0; t<m-1; t++) {
      table(m-2)(t).resize(m);
    }
    for (int s=0; s<m; s++) {
      for (int t=1; t<=(m/2); t++) {
	table(m-2)(2*t-2)(s)=-cos(2*M_PI*s*t/m);
      }
      for (int t=1; t<=((m+1)/2-1); t++) {
	table(m-2)(2*t-1)(s)=-sin(2*M_PI*s*t/m);
      }
    }
  }
}

GenericPlugIn<CorrFuncTable> *GenericPlugIn<CorrFuncTable>::list=NULL;
SpecificPlugIn<CorrFuncTable,TrigoCorrFuncTable> TrigoCorrFuncPlugIn("trigo");

void find_equivalent_clusters(Array<MultiCluster> *pclusters, const MultiCluster &cluster,
                      const rMatrix3d &cell, const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans) {
  rMatrix3d inv_cell=!cell;
  LinkedList<MultiCluster> list;
  for (int op=0; op<point_op.get_size(); op++) {
    MultiCluster *ptransfo_cluster=new MultiCluster(cluster.clus.get_size());
    apply_symmetry(ptransfo_cluster,point_op(op),trans(op),cluster);
    LinkedListIterator<MultiCluster> i(list);
    for ( ; i; i++) {
      if (equivalent_mod_cell(*i,*ptransfo_cluster,inv_cell)) break;
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

int calc_multiplicity(const MultiCluster &cluster,
                              const rMatrix3d &cell, const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans) {
  Array<MultiCluster> equiv;
  find_equivalent_clusters(&equiv, cluster,cell,point_op,trans);
  return equiv.get_size();
}

Real calc_correlation(const Structure &str, const Array<MultiCluster> &clusters, const rMatrix3d &cell, const Array<Array<Array<Real> > > &corrfunc) {
  Real accum=0.;
  int count=0;
  rMatrix3d inv_strcell=!str.cell;
  LatticePointInCellIterator t(cell,str.cell);
  for ( ;t; t++) {
    for (int c=0; c<clusters.get_size(); c++) {
      Real sigma=1.;
      for (int at=0; at<clusters(c).clus.get_size(); at++) {
	  Real spin=corrfunc(clusters(c).site_type(at))(clusters(c).func(at))(str.atom_type(which_atom(str.atom_pos,t+clusters(c).clus(at),inv_strcell)));
        sigma=sigma*spin;
      }
      accum+=sigma;
      count++;
    }
  }
  return accum/(Real)count;
}

Real calc_correlation(const Structure &str, const MultiCluster &cluster,
                      const rMatrix3d &cell, const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans, const Array<Array<Array<Real> > > &corrfunc) {
  Array<MultiCluster> clusters;
  find_equivalent_clusters(&clusters, cluster, cell, point_op, trans);
  return calc_correlation(str,clusters,cell,corrfunc);
}

void find_clusters_overlapping_site(Array<MultiCluster> *o_cluster, Array<int> *o_which,
				    const rVector3d &site,
				    const SpaceGroup &space_group,
				    const LinkedList<MultiCluster> &cluster_list) {
  rMatrix3d inv_cell=!space_group.cell;
  LinkedList<MultiCluster> long_cluster_list;
  LinkedList<int> long_which_list;
  LinkedListIterator<MultiCluster> c(cluster_list);
  int e=0;
  for ( ; c; c++,e++) {
    if (c->clus.get_size()>0) {
      Array<MultiCluster> clusters;
      find_equivalent_clusters(&clusters, *c,space_group.cell,space_group.point_op,space_group.trans);
      for (int ec=0; ec<clusters.get_size(); ec++) {
	for (int center_atom=0; center_atom<clusters(ec).clus.get_size(); center_atom++) {
	  rVector3d lat_shift=site-clusters(ec).clus(center_atom);
	  if (is_int(inv_cell*lat_shift)) {
	    MultiCluster *to_add=new MultiCluster(clusters(ec).clus.get_size());
	    int j=0;
	    to_add->clus(j)=site;
	    to_add->site_type(j)=clusters(ec).site_type(center_atom);
	    to_add->func(j)=clusters(ec).func(center_atom);
	    j++;
	    for (int i=0; i<clusters(ec).clus.get_size(); i++) {
	      if (i!=center_atom) {
		to_add->clus(j)=clusters(ec).clus(i)+lat_shift;
		to_add->site_type(j)=clusters(ec).site_type(i);
		to_add->func(j)=clusters(ec).func(i);
		j=j+1;
	      }
	    }
	    long_cluster_list << to_add;
	    long_which_list << new int(e);
	  }
	}
      }
    }  
  }
  LinkedList_to_Array(o_cluster,long_cluster_list);
  LinkedList_to_Array(o_which,long_which_list);
}

void calc_corr_to_conc(Array2d<Real> *pmat, const Structure &lat, const Array<Array<int> > &site_type_list, const SpaceGroup &spacegroup, const Array<Array<Array<Real> > > &corrfunc, int usemultcorr) {
  int total_species=0;
  for (int i=0; i<site_type_list.get_size(); i++) {
    total_species=MAX(total_species,1+max(site_type_list(i)));
  }
  /*
  Array<Real> max_conc_vector(total_species);
  zero_array(&max_conc_vector);
  for (int i=0; i<lat.atom_type.get_size(); i++) {
    for (int j=0; j<site_type_list(lat.atom_type(i)).get_size(); j++) {
      max_conc_vector(site_type_list(lat.atom_type(i))(j))+=1.;
    }
  }
  */
  Structure strnb;
  strnb=lat;
  for (int i=0; i<strnb.atom_type.get_size(); i++) {
    strnb.atom_type(i)=site_type_list(strnb.atom_type(i)).get_size();
  }
  MultiClusterBank point(strnb,1,spacegroup);
  rMatrix3d invcell=!(lat.cell);

  LinkedListIterator<MultiCluster> c(point.get_cluster_list());
  Array<int> nb_species(point.get_cluster_list().get_size());
  Array<Array<int> > equiv_site(nb_species.get_size());
  for (int i=0; i<nb_species.get_size(); i++, c++) {
    nb_species(i)=site_type_list(lat.atom_type(which_atom(lat.atom_pos,c->clus(0),invcell))).get_size();
    Array<MultiCluster> eq_clusters;
    find_equivalent_clusters(&eq_clusters, *c, lat.cell, spacegroup.point_op,  spacegroup.trans);
    equiv_site(i).resize(eq_clusters.get_size());
    for (int j=0; j<eq_clusters.get_size(); j++) {
      equiv_site(i)(j)=which_atom(lat.atom_pos,eq_clusters(j).clus(0),invcell);
    }
  }

  MultiDimIterator<Array<int> > endpt(nb_species);
  Array<Array<Real> > corr_basis;
  LinkedList<Array<Real> > conc_basis_list;
  Structure str;
  str=lat;
  zero_array(&(str.atom_type));
  
  int nb_in_basis=0;
  for ( ; endpt; endpt++) {
    for (int i=0; i<nb_species.get_size(); i++) {
      for (int j=0; j<equiv_site(i).get_size(); j++) {
	str.atom_type(equiv_site(i)(j))=((Array<int> &)endpt)(i);
      }
    }
    /*    cout << (Array<int> &)endpt << endl;
    {
      Array<AutoString> atom_label(5);
      atom_label(0)=AutoString("A");
      atom_label(1)=AutoString("B");
      atom_label(2)=AutoString("C");
      atom_label(3)=AutoString("D");
      atom_label(4)=AutoString("E");
      write_structure(str, lat, site_type_list, atom_label, lat.cell, cout);
      cout << "end" << endl << endl;
    }
    */
    LinkedList<Real> corr_list;
    LinkedListIterator<MultiCluster> c(point.get_cluster_list());
    corr_list << new Real(1.);
    for ( ; c; c++) {
      Real corr=calc_correlation(str,*c,spacegroup.cell,spacegroup.point_op,spacegroup.trans,corrfunc);
      if (usemultcorr) {corr*=calc_multiplicity(*c,spacegroup.cell,spacegroup.point_op,spacegroup.trans)/(Real)(lat.atom_type.get_size());}
      corr_list << new Real(corr);
    }
    Array<Real> cur_corr;
    LinkedList_to_Array(&cur_corr,corr_list);
    if (build_basis(&corr_basis,&nb_in_basis,cur_corr)) {
      Array<Real> conc_vector(total_species);
      zero_array(&conc_vector);
      for (int i=0; i<str.atom_type.get_size(); i++) {
	conc_vector(site_type_list(lat.atom_type(i))(str.atom_type(i)))+=1.;
      }
      for (int i=0; i<conc_vector.get_size(); i++) {
	conc_vector(i)/=(Real)(lat.atom_type.get_size());
      }
      conc_basis_list << new Array<Real>(conc_vector);
    }
  }

  /*
  MultiDimIterator<Array<int> > endpt(str.atom_type);
  Array<Array<Real> > corr_basis;
  LinkedList<Array<Real> > conc_basis_list;
  
  int nb_in_basis=0;
  for ( ; endpt; endpt++) {

    for (int i=0; i<str.atom_type.get_size(); i++) {
      str.atom_type(i)=((Array<int> &)endpt)(lat.atom_type(i));
    }

    str.atom_type=(Array<int> &)endpt;
    LinkedList<Real> corr_list;
    LinkedListIterator<MultiCluster> c(point.get_cluster_list());
    corr_list << new Real(1.);
    for ( ; c; c++) {
      Real corr=calc_correlation(str,*c,spacegroup.cell,spacegroup.point_op,spacegroup.trans,corrfunc);
      if (usemultcorr) {corr*=calc_multiplicity(*c,spacegroup.cell,spacegroup.point_op,spacegroup.trans)/(Real)(lat.atom_type.get_size());}
      corr_list << new Real(corr);
    }
    Array<Real> cur_corr;
    LinkedList_to_Array(&cur_corr,corr_list);
    if (build_basis(&corr_basis,&nb_in_basis,cur_corr)) {
      Array<Real> conc_vector(total_species);
      zero_array(&conc_vector);
      for (int i=0; i<str.atom_type.get_size(); i++) {
	conc_vector(site_type_list(lat.atom_type(i))(str.atom_type(i)))+=1.;
      }
      for (int i=0; i<conc_vector.get_size(); i++) {
	conc_vector(i)/=(Real)(lat.atom_type.get_size());
      }
      conc_basis_list << new Array<Real>(conc_vector);
    }
  }
  */

  Array2d<Real> corrmat(nb_in_basis,nb_in_basis);
  for (int i=0; i<nb_in_basis; i++) {
    for (int j=0; j<nb_in_basis; j++) {
      corrmat(i,j)=corr_basis(j)(i);
    }
  }
  Array2d<Real> icorrmat;
  invert_matrix(&icorrmat,corrmat);
  LinkedListIterator<Array<Real> > ic(conc_basis_list);
  Array2d<Real> concmat(ic->get_size(),nb_in_basis);
  for (int i=0; i<nb_in_basis; i++, ic++) {
    for (int j=0; j<ic->get_size(); j++) {
      concmat(j,i)=(*ic)(j);
    }
  }
  //  cerr << corrmat << endl;
  //  cerr << concmat << endl;
  product(pmat,concmat,icorrmat);
}

void extract_nonredundant(Array2d<Real> *pnored, Array2d<Real> *pnoredtored, Array<Real> *pconst, const Array2d<Real> &mat) {
  int dim=mat.get_size()(1);
  Array<Array<Real> > basis(dim);
  int nb_in_basis=1;
  basis(0).resize(dim);
  zero_array(&basis(0));
  basis(0)(0)=1.;
  for (int i=0; i<mat.get_size()(0); i++) {
    Array<Real> v;
    extract_row(&v,mat,i);
    v(0)=0.;
    build_basis(&basis,&nb_in_basis,v);
  }
  pnored->resize(nb_in_basis-1,dim-1);
  for (int i=0; i<nb_in_basis-1; i++) {
    for (int j=0; j<dim-1; j++) {
      (*pnored)(i,j)=basis(i+1)(j+1);
    }
  }
  Array2d<Real> bignored(nb_in_basis,dim);
  for (int i=0; i<nb_in_basis; i++) {
    for (int j=0; j<dim; j++) {
      bignored(i,j)=basis(i)(j);
    }
  }
  Array2d<Real> bignoredt,mat_bignoredt,bignored_bignoredt,inv,bignoredtored;
  transpose_matrix(&bignoredt,bignored);
  product(&mat_bignoredt,mat,bignoredt);
  product(&bignored_bignoredt,bignored,bignoredt);
  invert_matrix(&inv,bignored_bignoredt);
  product(&bignoredtored,mat_bignoredt,inv);
  extract_column(pconst,bignoredtored,0);
  Array<int> cols(bignoredtored.get_size()(1)-1);
  for (int i=0; i<cols.get_size(); i++) {
    cols(i)=i+1;
  }
  extract_columns(pnoredtored,bignoredtored,cols);
/*
  pnored->resize(dim-1,dim-1);
  for (int i=0; i<dim-1; i++) {
    for (int j=0; j<dim-1; j++) {
      (*pnored)(i,j)=basis(i+1)(j+1);
    }
  }
  Array2d<Real> inv,tmp;
  ArrayArray_to_Array2d(&tmp,basis,0);
  invert_matrix(&inv,tmp);
  product(&tmp,mat,inv);
  extract_column(pconst,tmp,0);
  Array<int> cols(dim-1);
  for (int i=0; i<cols.get_size(); i++) {
    cols(i)=i+1;
  }
  extract_columns(pnoredtored,tmp,cols);
*/
}

void calc_concentration(Array<Real> *pconc, const Array<int> &atom_type, int nbatomtype) {
  Real dc=1./(Real)atom_type.get_size();
  pconc->resize(nbatomtype);
  zero_array(pconc);
  for (int i=0; i<atom_type.get_size(); i++) {
    (*pconc)(atom_type(i))+=dc;
  }
}

int calc_concentration(Array<Real> *pconc, const Structure &lattice, const Array<Array<int> > &site_type_list, const Structure &str) {
  int comp=0;
  for (int i=0; i<site_type_list.get_size(); i++) {
    comp=MAX(comp,max(site_type_list(i)));
  }
  comp++;
  pconc->resize(comp);
  zero_array(pconc);
  rMatrix3d invcell=!lattice.cell;
  Real tot=0.;
  for (int i=0; i<str.atom_type.get_size(); i++) {
    int s=which_atom(lattice.atom_pos,str.atom_pos(i),invcell);
    if (s==-1) {return 0;}
    (*pconc)(site_type_list(lattice.atom_type(s))(str.atom_type(i)))+=1.;
    tot+=1.;
  }
  product(pconc,*pconc,1./tot);
  return 1;
}

void calc_noredtored(Array2d<Real> *pmat, Array<Real> *pvec, const Structure &lat, const Array<Array<int> > &site_type_list) {
  int total_species=0;
  for (int i=0; i<site_type_list.get_size(); i++) {
    total_species=MAX(total_species,1+max(site_type_list(i)));
  }
  Structure str;
  str=lat;
  Array<int> nb_species(site_type_list.get_size());
  for (int i=0; i<site_type_list.get_size(); i++) {
    nb_species(i)=site_type_list(i).get_size();
  }
  MultiDimIterator<Array<int> > endpt(nb_species);
  Array<Array<Real> > conc_basis;
  
  int nb_in_basis=0;
  Array<Real> ref(0);
  for ( ; endpt; endpt++) {
    for (int i=0; i<str.atom_type.get_size(); i++) {
      str.atom_type(i)=((Array<int> &)endpt)(lat.atom_type(i));
    }
    Array<Real> cur_conc;
    calc_concentration(&cur_conc,lat,site_type_list,str);
    if (ref.get_size()==0) {
      ref=cur_conc;
    }
    else {
      diff(&cur_conc,cur_conc,ref);
      build_basis(&conc_basis,&nb_in_basis,cur_conc);
    }
  }
  *pvec=ref;
  pmat->resize(ref.get_size(),nb_in_basis);
  for (int i=0; i<ref.get_size(); i++) {
    for (int j=0; j<nb_in_basis; j++) {
      (*pmat)(i,j)=conc_basis(j)(i);
    }
  }
}

void calc_redtonored(Array2d<Real> *pmat, const Array2d<Real> &noredtored, const Array<Real> &noredtored_c, int keepconst) {
  int dropconst=1-keepconst;
  Array2d<Real> x(noredtored.get_size()(0),1+noredtored.get_size()(1));
  for (int i=0; i<noredtored.get_size()(0); i++) {
    x(i,0)=noredtored_c(i);
    for (int j=0; j<noredtored.get_size()(1); j++) {
      x(i,1+j)=noredtored(i,j);
    }
  }
  Array2d<Real> xx,ixx;
  inner_product(&xx,x,x);
  invert_matrix(&ixx,xx);
  Array2d<Real> mat;
  product(&mat,x,ixx);
  pmat->resize(mat.get_size()(1)-dropconst,mat.get_size()(0));
  for (int i=0; i<pmat->get_size()(0); i++) {
    for (int j=0; j<pmat->get_size()(1); j++) {
      (*pmat)(i,j)=mat(j,dropconst+i);
    }
  }
}

void find_clusters_symmetry(SpaceGroup *psubgroup, const MultiCluster &cluster, const SpaceGroup &space_group) {
  rMatrix3d inv_cell=!(space_group.cell);
  LinkedList<rMatrix3d> point_list;
  LinkedList<rVector3d> trans_list;
  for (int op=0; op<space_group.point_op.get_size(); op++) {
    MultiCluster transfo_cluster;
    apply_symmetry(&transfo_cluster,space_group.point_op(op),space_group.trans(op),cluster);
    if (equivalent_mod_cell(transfo_cluster,cluster,inv_cell)) {
      rVector3d shift(0.,0.,0.);
      for (int i=0; i<cluster.clus.get_size(); i++) {
	shift+=(cluster.clus(i)-transfo_cluster.clus(i));
      }
      shift=shift/(Real)(cluster.clus.get_size());
      point_list << new rMatrix3d(space_group.point_op(op));
      trans_list << new rVector3d(space_group.trans(op)+shift);
    }
  }
  psubgroup->cell=space_group.cell;
  LinkedList_to_Array(&(psubgroup->point_op),point_list);
  LinkedList_to_Array(&(psubgroup->trans),trans_list);
}
