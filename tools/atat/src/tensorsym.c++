#include "tensorsym.h"
#include "lstsqr.h"
#include "phonlib.h"

void apply_symmetry(rTensor *ptt, const rMatrix3d &op, const rTensor &t) {
  ptt->resize(t.get_size());
  ptt->zero();
  int rank=t.get_size().get_size();
  MultiDimIterator<Array<int> > i(t.get_size());
  for (; i; i++) {
    MultiDimIterator<Array<int> > j(t.get_size());
    for (; j; j++) {
      Real acc=1.;
      for (int k=0; k<rank; k++) {
	acc*=op(((Array<int> &)i)(k),((Array<int> &)j)(k));
      }
      (*ptt)(i)+=t(j)*acc;
    }
  }
}

int symmetrize(rTensor *psymt, const rTensor &t, const SpaceGroup &space_group) {
  psymt->resize(t.get_size());
  psymt->zero();
  for (int i=0; i<space_group.point_op.get_size(); i++) {
    rTensor tt;
    apply_symmetry(&tt,space_group.point_op(i),t);
    (*psymt)+=tt;
  }
  Real n=norm((*psymt));
  if (fabs(n)>zero_tolerance) {
    (*psymt)/=n;
    return 1;
  }
  else {
    return 0;
  }
}

void transpose(rTensor *psymt, const rTensor &t, int row, int col) {
  psymt->resize(t.get_size());
  MultiDimIterator<Array<int> > i(t.get_size());
  for (; i; i++) {
    Array<int> ti(i);
    swap(&(ti(row)),&(ti(col)));
    (*psymt)(ti)=t(i);
  }
}

void transpose(rTensor *pt, int row, int col) {
  rTensor tmp(*pt);
  transpose(pt,tmp,row,col);
}

void transpose(rTensor *psymt, const rTensor &t, const Array<int> &map) {
  psymt->resize(t.get_size());
  Array<int> ti(map.get_size());
  MultiDimIterator<Array<int> > i(t.get_size());
  for (; i; i++) {
    Array<int> &ai=(Array<int> &)i;
    for (int k=0; k<map.get_size(); k++) {ti(k)=ai(map(k));}
    (*psymt)(ti)=t(ai);
  }
}

void calc_flip_group(Array<Array<int> > *pgroup, int rank, const Array<Array<int> > &flip) {
  LinkedList<Array<int> > group;
  Array<int> map(rank);
  for (int j=0; j<rank; j++) {map(j)=j;}
  group << new Array<int>(map);
  for (int f=0; f<flip.get_size(); f++) {
    LinkedListIterator<Array<int> > i(group);
    int n=group.get_size();
    for (; n>0; i++,n--) {
      map=*i;
      for (int p=0; p<flip(f).get_size(); p+=2) {
	swap(&(map(flip(f)(p))),&(map(flip(f)(p+1))));
      }
      LinkedListIterator<Array<int> > j(group);
      for ( ; j; j++) {
	if ((*j)==map) break;
      }
      if (!j) {
	group << new Array<int>(map);
      }
    }
  }
  LinkedList_to_Array(pgroup,group);
}

void symmetrize_flip(rTensor *psymt, const Array<Array<int> > &flipgroup) { 
  rTensor t(*psymt);
  for (int f=1; f<flipgroup.get_size(); f++) {
    rTensor ft;
    transpose(&ft,t,flipgroup(f));
    (*psymt)+=ft;
  }
}

/*
void symmetrize(rTensor *psymt, const Array<Array<int> > &flip) { 
  for (int f=0; f<flip.get_size(); f++) {
    rTensor ft(*psymt);
    for (int i=0; i<flip(f).get_size(); i+=2) {
      transpose(&ft,flip(f)(i),flip(f)(i+1));
    }
    (*psymt)+=ft;
  }
}
*/

/*
int symmetrize(rTensor *psymt, const rTensor &t, const SpaceGroup &space_group, const Array<Array<int> > &flip) { 
  if (!symmetrize(psymt, t,space_group)) {
    return 0;
  }
  symmetrize(psymt, flip);
  Real n=norm((*psymt));
  if (fabs(n)>zero_tolerance) {
    (*psymt)/=n;
    return 1;
  }
  else {
    return 0;
  }
}
*/

void find_symmetric_basis(Array<rTensor> *pbasis, int rank, const SpaceGroup &space_group, const Array<Array<int> > &flip) {
  Array<Array<int> > flipgroup;
  calc_flip_group(&flipgroup, rank,flip);  
  int dim=3;
  Array<Array<Real> > basis;
  int nb_in_basis=0;
  Array<int> size(rank);
  for (int i=0; i<size.get_size(); i++) {size(i)=dim;}
  rTensor t(size);
  rTensor done(size);
  done.zero();
  MultiDimIterator<Array<int> > i(size);
  for (; i; i++) {
    if (done(i)==0.) {
      t.zero();
      t(i)=1.;
      symmetrize_flip(&t,flipgroup);
      rTensor st;
      if (symmetrize(&st,t,space_group)) {
	build_basis(&basis,&nb_in_basis,st.vectorize());
      }
      done+=t;
    }
  }
  gram_schmidt(&basis,basis);
  pbasis->resize(nb_in_basis);
  for (int j=0; j<nb_in_basis; j++) {
    (*pbasis)(j).resize(size);
    (*pbasis)(j).vectorize()=basis(j);
  }
}

void print_harm(const rTensor &t) {
  MultiDimIterator<Array<int> > j(t.get_size());
  for (; j; j++) {
    Array<int> &aj=(Array<int> &)j;
    if (fabs(t(aj))>zero_tolerance) {
      Array<int> p(3);
      zero_array(&p);
      int s=1;
      for (int k=0; k<aj.get_size(); k++) {p(aj(k))++;}
      for (int k=0; k<aj.get_size()-1; k++) {if (aj(k)>aj(k+1)) {s=0;}}
      char c[]="xyz";
      if (s) {
	cerr << " + " << t(aj)*factorial(t.get_size().get_size())/(factorial(p(0))*factorial(p(1))*factorial(p(2)));
	for (int k=0; k<3; k++) {if (p(k)>0) {cerr << "*" << c[k] << "^" << p(k);}}
      }
    }
  }
  cerr << endl;
}

void find_symmetric_basis_projection(Array<rTensor> *pbasis, int rank, const SpaceGroup &space_group, const Array<Array<int> > &flip, const Array<rTensor> &proj, int paral ) {
  if (proj.get_size()==0) {
    if (paral==0) {find_symmetric_basis(pbasis, rank,space_group,flip);}
    return;
  }
  Array<Array<int> > flipgroup;
  calc_flip_group(&flipgroup, rank,flip);
  int dim=3;
  Array<Array<Real> > basis;
  int nb_in_basis=0;
  Array<int> size(rank);
  for (int i=0; i<size.get_size(); i++) {size(i)=dim;}
  rTensor t(size);
  rTensor done(size);
  done.zero();

  Array2d<Real> mat(t.vectorize().get_size(),proj.get_size());
  for (int j=0; j<proj.get_size(); j++) {
    for (int k=0; k<proj(j).vectorize().get_size(); k++) {
      mat(k,j)=proj(j).vectorize()(k);
    }
  }
  // cerr << mat;
  Array2d<Real> projmat,imat;
  invert_matrix_tall(&imat,mat);
  product(&projmat, mat,imat);
  if (!paral) {
    product(&projmat,projmat,-1.);
    for (int j=0; j<projmat.get_size()(0); j++) {
      projmat(j,j)+=1.;
    }
  }
  //  cerr << projmat;
  MultiDimIterator<Array<int> > i(size);
  for (; i; i++) {
    if (done(i)==0.) {
      t.zero();
      t(i)=1.;
      symmetrize_flip(&t,flipgroup);
      rTensor st;
      if (symmetrize(&st,t,space_group)) {
	rTensor pt(size);
	product(&(pt.vectorize()),projmat,st.vectorize());
	if (norm(pt)>zero_tolerance) {
	  build_basis(&basis,&nb_in_basis,pt.vectorize());
	}
      }
      done+=t;
    }
  }
  gram_schmidt(&basis,basis);
  pbasis->resize(nb_in_basis);
  for (int j=0; j<nb_in_basis; j++) {
    (*pbasis)(j).resize(size);
    (*pbasis)(j).vectorize()=basis(j);
  }
}

void all_flip(Array<Array<int> > *pflip, int rank) {
  pflip->resize(rank*(rank-1)/2);
  int f=0;
  for (int i=0; i<rank; i++) {
    for (int j=0; j<i; j++) {
      (*pflip)(f).resize(2);
      (*pflip)(f)(0)=i;
      (*pflip)(f)(1)=j;
      f++;
    }
  }
}

void calc_sym_harmonics(Array<rTensor> *pbasis, int rank, const SpaceGroup &space_group) {
  int dim=3;
  Array<int> size(rank);
  for (int i=0; i<size.get_size(); i++) {size(i)=dim;}

  Array<Array<int> > flip;
  all_flip(&flip,rank);
  Array<Array<int> > flipgroup;
  calc_flip_group(&flipgroup, rank,flip);
  Array<rTensor> proj_2;
  if (rank>2) {
    Array<Array<int> > flip_2;
    all_flip(&flip_2,rank-2);
    Array<Array<int> > flipgroup_2;
    calc_flip_group(&flipgroup_2, rank,flip_2);
    Array<rTensor> basis_2;
    find_symmetric_basis(&basis_2, rank-2,space_group,flip_2);
    proj_2.resize(basis_2.get_size());
    for (int i=0; i<proj_2.get_size(); i++) {
      rTensor tmp(size);
      MultiDimIterator<Array<int> > j(size);
      for (; j; j++) {
	Array<int> subj;
	extract_elements(&subj,(Array<int> &)j,2,rank);
	Array<int> &aij=(Array<int> &)j;
	tmp(j)=basis_2(i)(subj)*(aij(0)==aij(1));
      }
      //cerr << flip << endl;
      //print_harm(tmp);
      //      apply_symmetry(&(proj_2(i)),space_group.point_op(1),tmp);
      //      cerr << space_group.point_op(1) << endl;
      //      print_harm(proj_2(i));
      //      exit(1);
      //      cerr << space_group.point_op.get_size() << endl;
      symmetrize_flip(&tmp,flipgroup);
      symmetrize(&(proj_2(i)),tmp,space_group);
      //print_harm(proj_2(i));
      //exit(1);
      //cerr << proj_2(i).get_size() << proj_2(i).vectorize();
    }
  }
  else if (rank==2) {
    proj_2.resize(1);
    proj_2(0).resize(size);
    MultiDimIterator<Array<int> > j(size);
    for (; j; j++) {
      Array<int> &aij=(Array<int> &)j;
      proj_2(0)(j)=(aij(0)==aij(1));
    }
  }

  //    find_symmetric_basis(pbasis, rank,space_group,flip);
  find_symmetric_basis_projection(pbasis, rank,space_group,flip,proj_2,0);
}

void find_symmetry_breaking_basis(Array<rTensor> *pbasis, int rank, const SpaceGroup &space_group, const Array<Array<int> > &flip) {
  Array<Array<int> > flipgroup;
  calc_flip_group(&flipgroup, rank,flip);
  Array<rVector3d> special_dir;
  Array<rVector3d> all_dir(3);
  find_special_direction(&special_dir,NULL,space_group.point_op);
  for (int i=0; i<special_dir.get_size(); i++) {
    all_dir(i)=special_dir(i);
  }
  if (special_dir.get_size()==1) {
    all_dir(1)=find_perpendicular(all_dir(0));
  }
  if (special_dir.get_size()<3) {
    all_dir(2)=all_dir(0)^all_dir(1);
  }

  int dim=3;
  LinkedList<rTensor> tensor_list;
  Array<Array<Real> > basis;
  int nb_in_basis=0;
  Array<int> size(rank);
  for (int i=0; i<size.get_size(); i++) {size(i)=dim;}
  rTensor t(size);
  MultiDimIterator<Array<int> > i(size);
  for (; i; i++) {
    MultiDimIterator<Array<int> > j(size);
    for (; j; j++) {
      Real a=1.;
      for (int k=0; k<rank; k++) {
	a*=all_dir(((Array<int>&)i)(k))(((Array<int>)j)(k));
      }
      t(j)=a;
    }
    symmetrize_flip(&t,flipgroup);
    int added=0;
    for (int op=0; op<space_group.point_op.get_size(); op++) {
      rTensor st;
      apply_symmetry(&st,space_group.point_op(op),t);
      if (build_basis(&basis,&nb_in_basis,st.vectorize())) {
	added=1;
      }
    }
    if (added) {
      normalize(&(t.vectorize()));
      tensor_list << new rTensor(t);
    }
  }
  LinkedList_to_Array(pbasis,tensor_list);
}

