#include "findsym.h"

#define HARTFORCADE

void apply_symmetry(Array<rVector3d> *b, const rMatrix3d &point_op, const rVector3d &trans, const Array<rVector3d> &a) {
  b->resize(a.get_size());
  for (int i=0; i<a.get_size(); i++) {
    (*b)(i)=point_op*a(i)+trans;
  }
}

void apply_symmetry(MultiCluster *b, const rMatrix3d &point_op, const rVector3d &trans, const MultiCluster &a) {
  b->clus.resize(a.clus.get_size());
  b->site_type=a.site_type;
  b->func=a.func;
  for (int i=0; i<a.clus.get_size(); i++) {
    (b->clus)(i)=point_op*a.clus(i)+trans;
  }
}

int equivalent_by_symmetry(const rVector3d &a, const rVector3d &b,
                           const rMatrix3d &cell,
                           const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans) {
  rMatrix3d inv_cell=!cell;
  for (int op=0; op<point_op.get_size(); op++) {
    rVector3d t_b=point_op(op)*b+trans(op);
    if (equivalent_mod_cell(a,t_b,inv_cell)) return 1;
  }
  return 0;
}

int equivalent_by_symmetry(const Array<rVector3d> &a, const Array<rVector3d> &b,
                           const rMatrix3d &cell,
                           const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans) {
  rMatrix3d inv_cell=!cell;
  for (int op=0; op<point_op.get_size(); op++) {
    Array<rVector3d> t_b;
    apply_symmetry(&t_b,point_op(op),trans(op),b);
    if (equivalent_mod_cell(a,t_b,inv_cell)) return 1;
  }
  return 0;
}

int equivalent_by_symmetry(const MultiCluster &a, const MultiCluster &b,
                           const rMatrix3d &cell,
                           const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans) {
  rMatrix3d inv_cell=!cell;
  for (int op=0; op<point_op.get_size(); op++) {
    MultiCluster t_b;
    apply_symmetry(&t_b.clus,point_op(op),trans(op),b.clus);
    t_b.func=b.func;
    if (equivalent_mod_cell(a,t_b,inv_cell)) return 1;
  }
  return 0;
}

void calc_concentration(int *pmin_type, int *pmax_type, Array<Real> *pconc, const Array<int> &atom_type) {
  Real dc=1./(Real)atom_type.get_size();
  *pmin_type=min(atom_type);
  *pmax_type=max(atom_type);
  pconc->resize((*pmax_type)-(*pmin_type)+1);
  zero_array(pconc);
  for (int i=0; i<atom_type.get_size(); i++) {
    (*pconc)(atom_type(i)-(*pmin_type))+=dc;
  }
}

Real check_sum_structure(const Array<int> &atom_type) {
  int sum=0;
  for (int i=0; i<atom_type.get_size(); i++) {sum+=atom_type(i);}
  return (Real)sum/(Real)atom_type.get_size();
}

int equivalent_by_symmetry(const Structure &a, const Structure &b,
                           const rMatrix3d &cell,
                           const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans) {
  if (a.atom_pos.get_size() > b.atom_pos.get_size()) return equivalent_by_symmetry(b,a,cell,point_op,trans);
  if ((b.atom_pos.get_size() % a.atom_pos.get_size())!=0) return 0;
  //  if (!near_zero(check_sum_structure(a.atom_type)-check_sum_structure(b.atom_type))) return 0;

  Array<Real> conca,concb;
  int min_typea,max_typea,min_typeb,max_typeb;
  calc_concentration(&min_typea,&max_typea,&conca,a.atom_type);
  calc_concentration(&min_typeb,&max_typeb,&concb,b.atom_type);
  if (min_typea!=min_typeb || max_typea!=max_typeb) return 0;
  for (int i=0; i<conca.get_size(); i++) {
    if (fabs(conca(i)-concb(i))>zero_tolerance) return 0;
  }

  rMatrix3d inv_acell=!a.cell;
  for (int op=0; op<point_op.get_size(); op++) {
    Structure tb;
    tb.cell=point_op(op)*b.cell;
    apply_symmetry(&tb.atom_pos,point_op(op),trans(op),b.atom_pos);
    if (is_int(inv_acell*tb.cell)) {
      LatticePointInCellIterator shift(cell,a.cell);
      for (;shift; shift++) {
        int at;
        for (at=0; at<b.atom_pos.get_size(); at++) {
          int at2=which_atom(a.atom_pos,tb.atom_pos(at)+shift,inv_acell);
          if (at2==-1) break;
          if (a.atom_type(at2)!=b.atom_type(at)) break;
        }
        if (at==b.atom_pos.get_size()) return 1;
      }
    }
  }
  return 0;
}

void find_pointgroup(Array<rMatrix3d> *point_op, const rMatrix3d &cell) {
  LinkedList<rMatrix3d> point_op_list;
  rVector3d u1=cell.get_column(0);
  rVector3d u2=cell.get_column(1);
  rVector3d u3=cell.get_column(2);
  rMatrix3d inv_cell=!cell;
  rMatrix3d id;
  id.identity();

  LatticePointIterator p1(cell);
  while (norm(p1)<norm(u1)-zero_tolerance) p1++;
  for (;near_zero(norm(p1)-norm(u1)); p1++) {
    LatticePointIterator p2(cell);
    while (norm(p2)<norm(u2)-zero_tolerance) p2++;
    for (;near_zero(norm(p2)-norm(u2)); p2++) {
      if (near_zero(u1*u2-p1*p2)) {
        LatticePointIterator p3(cell);
        while (norm(p3)<norm(u3)-zero_tolerance) p3++;
        for (;near_zero(norm(p3)-norm(u3)); p3++) {
	  if (near_zero(u2*u3-p2*p3) && near_zero(u1*u3-p1*p3)) {
	    rMatrix3d new_cell;
	    new_cell.set_column(0,p1);
	    new_cell.set_column(1,p2);
	    new_cell.set_column(2,p3);
	    rMatrix3d t=new_cell*(inv_cell);
            point_op_list << new rMatrix3d(t);
          }
        }
      }
    }
  }
  LinkedList_to_Array(point_op,point_op_list);
}

void pointgroup_from_spacegroup(Array<rMatrix3d> *pointgroup, const Array<rMatrix3d> &point_op) {
  LinkedList<rMatrix3d> list;
  for (int op=0; op<point_op.get_size(); op++) {
    LinkedListIterator<rMatrix3d> i(list);
    while (i) {
      if (max_norm((*i)-point_op(op))<zero_tolerance) break;
      i++;
    }
    if (!i) {
      list << new rMatrix3d(point_op(op));
    }
  }
  LinkedList_to_Array(pointgroup,list);
}

void pointgroup_from_spacegroup(const Structure &str, Array<rMatrix3d> *pointgroup, const rMatrix3d &cell, const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans) {
  LinkedList<rMatrix3d> list;
  rMatrix3d inv_cell=!str.cell;
  for (int op=0; op<point_op.get_size(); op++) {
    if (is_int((inv_cell)*point_op(op)*str.cell)) {
      LatticePointInCellIterator shift(cell,str.cell);
      for ( ; shift; shift++) {
	Array<rVector3d> transfo_pos;
	apply_symmetry(&transfo_pos,point_op(op),trans(op)+shift,str.atom_pos);
	int at;
	for (at=0; at<str.atom_pos.get_size(); at++) {
	  if (str.atom_type(which_atom(str.atom_pos,transfo_pos(at),inv_cell))!=str.atom_type(at)) break;
	}
	if (at==str.atom_pos.get_size()) break;
      }
      if (shift) {
	LinkedListIterator<rMatrix3d> i(list);
	while (i) {
	  if (max_norm((*i)-point_op(op))<zero_tolerance) break;
	  i++;
	}
	if (!i) {
	  list << new rMatrix3d(point_op(op));
	}
      }
    }
  }
  LinkedList_to_Array(pointgroup,list);
}

void find_primitive_unit_cell(rMatrix3d *pcell, const Structure &lat) {
  rVector3d a[3];
  Structure t_lat;
  t_lat.cell=lat.cell;
  t_lat.atom_type=lat.atom_type;
  rMatrix3d point_op;
  point_op.identity();
  Array<rMatrix3d> point_ops(1);
  Array<rVector3d> trans(1);
  point_ops(0)=point_op;
  trans(0)=rVector3d(0.,0.,0.);
  AtomPairIterator p(lat.cell,lat.atom_pos(0),lat.atom_pos);
  for (; ; p++) {
    rVector3d t=p(1)-p(0);
    apply_symmetry(&t_lat.atom_pos,point_op,t,lat.atom_pos);
    if (equivalent_by_symmetry(t_lat,lat,lat.cell,point_ops,trans)) {
      a[0]=t;
      break;
    }
  }
  for (; ; p++) {
    rVector3d t=p(1)-p(0);
    apply_symmetry(&t_lat.atom_pos,point_op,t,lat.atom_pos);
    if (equivalent_by_symmetry(t_lat,lat,lat.cell,point_ops,trans)) {
      if (!near_zero(norm(a[0]^t))) {
	a[1]=t;
	break;
      }
    }
  }
  for (; ; p++) {
    rVector3d t=p(1)-p(0);
    apply_symmetry(&t_lat.atom_pos,point_op,t,lat.atom_pos);
    if (equivalent_by_symmetry(t_lat,lat,lat.cell,point_ops,trans)) {
      if ((a[0]^a[1])*t>zero_tolerance) {
	a[2]=t;
	break;
      }
    }
  }
  for (int i=0; i<3; i++) {pcell->set_column(i,a[i]);}
}

void find_spacegroup(Array<rMatrix3d> *point_op, Array<rVector3d> *trans,
                      const rMatrix3d &cell,
                      const Array<rVector3d> &atom_pos, const Array<int>
&atom_type) {
  rMatrix3d inv_cell=!cell;
  int max_type,min_type;
  Array<Real> atom_type_conc;
  calc_concentration(&min_type,&max_type,&atom_type_conc,atom_type);
  int rare_type=-1;
  int min_conc=MAXINT;
  for (int i=0; i<atom_type_conc.get_size(); i++) {
    if (atom_type_conc(i)!=0) {
      if (atom_type_conc(i)<min_conc) {
	rare_type=min_type+i;
      }
    }
  }
  int ref_atom=0;
  while (atom_type(ref_atom)!=rare_type) {ref_atom++;}
  Array<rMatrix3d> pointgroup;
  {
    rMatrix3d prim_cell;
    Structure tmplat;
    tmplat.cell=cell;
    tmplat.atom_pos=atom_pos;
    tmplat.atom_type=atom_type;
    find_primitive_unit_cell(&prim_cell,tmplat); 
    find_pointgroup(&pointgroup,prim_cell);
  }
  //  find_pointgroup(&pointgroup,cell);
  LinkedList<rMatrix3d> point_op_list;
  LinkedList<rVector3d> trans_list;
  for (int i=ref_atom; i<atom_pos.get_size(); i++) {
    if (atom_type(i)==rare_type) {
      for (int op=0; op<pointgroup.get_size(); op++) {
        rVector3d trans=-pointgroup(op)*atom_pos(ref_atom)+atom_pos(i);
	  int at;
        for (at=0; at<atom_pos.get_size(); at++) {
          int at2=which_atom(atom_pos,pointgroup(op)*atom_pos(at)+trans,inv_cell);
          if (at2<0) break;
          if (atom_type(at2)!=atom_type(at)) break;
        }
        if (at==atom_pos.get_size()) {
          point_op_list << new rMatrix3d(pointgroup(op));
          trans_list << new rVector3d(trans);
        }
      }
    }
  }
  LinkedList_to_Array(point_op,point_op_list);
  if (trans) {
    LinkedList_to_Array(trans,trans_list);
  }
  else {
    pointgroup_from_spacegroup(point_op,*point_op);
  }
}

rMatrix3d find_almost_reduced_cell(const rMatrix3d &cell) {
  LatticePointIterator lp(cell);
  rVector3d a,b,c;
  a=lp;
  do {
    lp++;
  } while (norm(a^lp)<zero_tolerance);
  b=lp;
  do {
    lp++;
  } while (fabs((a^b)*lp)<zero_tolerance);
  c=lp;
  rMatrix3d rcell;
  rcell.set_column(0,a);
  rcell.set_column(1,b);
  rcell.set_column(2,c);
  iVector3d signs((b*c>0),(c*a>0),(a*b>0));
  int sign_to_flip=(l1_norm(signs) % 2);
  for (int i=0; i<3; i++) {
    if (signs(i)==sign_to_flip) {
      rcell.set_column(i,-rcell.get_column(i));
    }
  }
  if (det(rcell)<0) rcell=-rcell;
  return rcell;
}

void find_all_equivalent_cell(Array<rMatrix3d> *psupercell, const rMatrix3d &cell, Real radius) {
  LinkedList<rMatrix3d> cell_list;
  Real d=fabs(det(cell));
  LatticePointIterator lp[3];
  rMatrix3d curcell;
  Real radius2=radius*radius+zero_tolerance;
  lp[2].init(cell,1);
  while (norm2(lp[2])<=radius2) {
    lp[1].init(cell,1);
    while (norm2(lp[1])<=radius2) {
      lp[0].init(cell,1);
      while (norm2(lp[0])<=radius2) {
	for  (int i=0; i<3; i++) {
	  curcell.set_column(i,lp[i]);
	}
	if (near_zero(fabs(det(curcell))-d)) {
	  cell_list << new rMatrix3d(curcell);
	}
	lp[0]++;
      }
      lp[1]++;
    }
    lp[2]++;
  }
  LinkedList_to_Array(psupercell,cell_list);
}



void find_supercells(Array<rMatrix3d> *supercell, int min_volume, int max_volume, const rMatrix3d &unitcell, const Array<rMatrix3d> &pointgroup) {
  Array<iMatrix3d> A(pointgroup.get_size());
  for (int op=0; op<A.get_size(); op++) {
    A(op)=to_int ( (!unitcell)*pointgroup(op)*unitcell );
  }

  LinkedList<iMatrix3d> list_abc;

  for (int volume=min_volume; volume<=max_volume; volume++) {
    Real maxa=pow((Real)volume*sqrt(2.),1./3.);
    int imaxa=(int)ceil(maxa);
    LinkedList<iVector3d> list_a;
    MultiDimIterator<iVector3d> a(iVector3d(0,-imaxa,-imaxa),iVector3d(imaxa,imaxa,imaxa));
    for ( ; a; a++) {
      if (norm((iVector3d &)a)<=maxa+zero_tolerance) {
	LinkedListIterator<iVector3d> cur_a(list_a);
	while (cur_a) {
	  int sym;
	  for (sym=0; sym<A.get_size(); sym++) {
	    if ((A(sym)*(iVector3d &)a) == (*cur_a) || (A(sym)*(iVector3d &)a) == -(*cur_a)) break;
	  }
	  if (sym<A.get_size()) break;
	  cur_a++;
	}
	if (!cur_a && norm2((iVector3d &)a)!=0) {
	  list_a << new iVector3d((iVector3d &)a);
	}
      }
    }

    LinkedListIterator<iVector3d> cur_a(list_a);
    for ( ; cur_a; cur_a++) {
      iVector3d a=*cur_a;
      Real len_a=norm(a);
      Real maxb=sqrt(sqr(len_a)/6.+sqrt(sqr(sqr(len_a)/6.) + 4./3.*sqr((Real)volume/len_a)));
      int imaxb=(int)ceil(maxb);
      MultiDimIterator<iVector3d> b(iVector3d(0,-imaxb,-imaxb),iVector3d(imaxb,imaxb,imaxb));
      for ( ; b; b++) {
        if (norm((iVector3d &)b)<=maxb+zero_tolerance && norm((iVector3d &)a)>=len_a-zero_tolerance) {
	  Real proj=(Real)((iVector3d &)b*(iVector3d &)a)/(Real)((iVector3d &)a*(iVector3d &)a);
	  iVector3d axb=(iVector3d &)a^(iVector3d &)b;
	  int axb2=norm2(axb);
	  if (proj>=(-0.5+   zero_tolerance) &&
	      proj<=( 0.5+2.*zero_tolerance) &&
	      axb2!=0) {
	    int minw=-axb2/2;
	    int maxw=-axb2/2+axb2;
	    for (int wa=minw; wa<maxw; wa++) {
	      for (int wb=minw; wb<maxw; wb++) {
		iVector3d num=volume*axb+wa*(iVector3d &)a+wb*(iVector3d &)b;
		int i;
		for (i=0; i<3; i++) {
		  if ( (num(i) % axb2) != 0 ) break;
		}
		if (i==3) {
		  iVector3d c=num/axb2;
		  iMatrix3d abc;
		  abc.set_column(0,(iVector3d &)a);
		  abc.set_column(1,(iVector3d &)b);
		  abc.set_column(2,(iVector3d &)c);
		  LinkedListIterator<iMatrix3d> cur_abc(list_abc);
		  while (cur_abc) {
		    int sym;
		    for (sym=0; sym<A.get_size(); sym++) {
		      rMatrix3d new_abc=to_real(A(sym)*abc);
		      rMatrix3d old_abc=to_real(*cur_abc);
		      if (is_int((!new_abc)*old_abc) && is_int((!old_abc)*new_abc)) break;
		    }
		    if (sym<A.get_size()) break;
		    cur_abc++;
		  }
		  if (!cur_abc) {
		    list_abc << new iMatrix3d(abc);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  supercell->resize(list_abc.get_size());
  LinkedListIterator<iMatrix3d> cur_abc(list_abc);
  for (int i=0; i<supercell->get_size(); i++,cur_abc++) {
    (*supercell)(i)=unitcell*to_real(*cur_abc);
  }
}

void find_supercells_2D(Array<rMatrix3d> *supercell, int min_volume, int max_volume, const rMatrix3d &unitcell, const Array<rMatrix3d> &pointgroup) {
  rVector3d rec=unitcell.get_column(0)^unitcell.get_column(1);
  rVector3d rc=unitcell.get_column(2);
  if (!near_zero(rec*rc-norm(rec)*norm(rc))) {
    ERRORQUIT("c axis is not perpendicular to a,b plane.");
  }
  Array<iMatrix3d> A(pointgroup.get_size());
  for (int op=0; op<A.get_size(); op++) {
    A(op)=to_int ( (!unitcell)*pointgroup(op)*unitcell );
  }
  iVector3d c(0,0,1);

  LinkedList<iMatrix3d> list_abc;

  for (int volume=min_volume; volume<=max_volume; volume++) {
    Real maxa=pow((Real)volume*2./sqrt(3.),1./2.);
    int imaxa=(int)ceil(maxa);
    LinkedList<iVector3d> list_a;
    MultiDimIterator<iVector3d> a(iVector3d(0,-imaxa,0),iVector3d(imaxa,imaxa,0));
    for ( ; a; a++) {
      if (norm((iVector3d &)a)<=maxa+zero_tolerance) {
	LinkedListIterator<iVector3d> cur_a(list_a);
	while (cur_a) {
	  int sym;
	  for (sym=0; sym<A.get_size(); sym++) {
	    if ((A(sym)*(iVector3d &)a) == (*cur_a) || (A(sym)*(iVector3d &)a) == -(*cur_a)) break;
	  }
	  if (sym<A.get_size()) break;
	  cur_a++;
	}
	if (!cur_a && norm2((iVector3d &)a)!=0) {
	  list_a << new iVector3d((iVector3d &)a);
	}
      }
    }

    LinkedListIterator<iVector3d> cur_a(list_a);
    for ( ; cur_a; cur_a++) {
      iVector3d a=*cur_a;
      iVector3d cxa=c^a;
      int cxa2=norm2(cxa);
      int minw=-cxa2/2;
      int maxw=-cxa2/2+cxa2;
      for (int wa=minw; wa<maxw; wa++) {
	iVector3d num=volume*cxa+wa*a;
	int i;
	for (i=0; i<3; i++) {
	  if ( (num(i) % cxa2) != 0 ) break;
	}
	if (i==3) {
	  iVector3d b=num/cxa2;
	  iMatrix3d abc;
	  abc.set_column(0,(iVector3d &)a);
	  abc.set_column(1,(iVector3d &)b);
	  abc.set_column(2,(iVector3d &)c);
	  LinkedListIterator<iMatrix3d> cur_abc(list_abc);
	  while (cur_abc) {
	    int sym;
	    for (sym=0; sym<A.get_size(); sym++) {
	      rMatrix3d new_abc=to_real(A(sym)*abc);
	      rMatrix3d old_abc=to_real(*cur_abc);
	      if (is_int((!new_abc)*old_abc) && is_int((!old_abc)*new_abc)) break;
	    }
	    if (sym<A.get_size()) break;
	    cur_abc++;
	  }
	  if (!cur_abc) {
	    list_abc << new iMatrix3d(abc);
	  }
	}
      }
    }
  }
  supercell->resize(list_abc.get_size());
  LinkedListIterator<iMatrix3d> cur_abc(list_abc);
  for (int i=0; i<supercell->get_size(); i++,cur_abc++) {
    (*supercell)(i)=unitcell*to_real(*cur_abc);
  }
}

void sort_supercells(Array<rMatrix3d> *p_supercell, const Array<rMatrix3d> &pointgroup) {
  Array<int> nb_sym_op(p_supercell->get_size());
  zero_array(&nb_sym_op);
  for (int i=0; i<nb_sym_op.get_size(); i++) {
    rMatrix3d &cell=(*p_supercell)(i);
    rMatrix3d icell=!cell;
    for (int j=0; j<pointgroup.get_size(); j++) {
      if (is_int(icell*pointgroup(j)*cell)) nb_sym_op(i)++;
    }
  }
  for (int i=1; i<nb_sym_op.get_size(); i++) {
    for (int j=nb_sym_op.get_size()-1; j>=i; j--) {
      if (nb_sym_op(j)>nb_sym_op(j-1)) {
	swap(&nb_sym_op(j),&nb_sym_op(j-1));
	swap(&(*p_supercell)(j),&(*p_supercell)(j-1));
      }
    }
  }
}

int find_point_op_type(rVector3d *pspecial_dir, const rMatrix3d &point_op) {
  rVector3d &special_dir=*pspecial_dir;
  Real rdet=det(point_op);
  int det=iround(rdet);
  Real rtrace=trace(point_op);
  int trace=iround(rtrace);
  if (fabs(rdet-(Real)det)>zero_tolerance || fabs(rtrace-(Real)trace)>zero_tolerance ||
      abs(trace)>3 || abs(det)>1) {
    ERRORQUIT("get_point_op_type: Not a crystallographic symmetry");
  }
  int sym_type;
  switch (abs(trace)) {
  case 3:
    if (det==1) {sym_type=1;} else {sym_type=-1;}
    break;
  case 2:
    if (det==1) {sym_type=6;} else {sym_type=-6;}
    break;
  case 0:
    if (det==1) {sym_type=3;} else {sym_type=-3;}
    break;
  case 1:
    if (trace==1) {
      if (det==1) {sym_type=4;} else {sym_type=-2;}
    }
    else {
      if (det==1) {sym_type=2;} else {sym_type=-4;}
    }
  }
  if (sym_type==-2) {
    for (int axis=0; axis<3; axis++) {
      rVector3d v(0.,0.,0.);
      v(axis)=1.;
      special_dir=v-point_op*v;
      if (norm(special_dir)>zero_tolerance) {
	special_dir.normalize();
	break;
      }
    }
  }
  else {
    int step=1,period=0;
    switch (sym_type) {
    case 1:
    case -1:
      break;
    case 2:
      period=2;
      break;
    case 4:
    case -4:
      step=2;
      period=4;
      break;
    case 3:
      period=3;
      break;
    case 6:
      step=3;
      period=6;
      break;
    case -3:
    case -6:
      step=2;
      period=6;
      break;
    }
    for (int axis=0; axis<3; axis++) {
      special_dir=rVector3d(0.,0.,0.);
      rVector3d v(0.,0.,0.);
      v(axis)=1.;
      for (int i=0; i<period; i++) {
	if ((i % step) == 0) special_dir+=v;
	v=point_op*v;
      }
      if (norm(special_dir)>zero_tolerance) {
	special_dir.normalize();
	break;
      }
    }
  }
  return sym_type;
}

int number_of_point_op(int point_op_type, const Array<rMatrix3d> &point_op) {
  int nb=0;
  rVector3d dir;
  for (int op=0; op<point_op.get_size(); op++) {
    if (find_point_op_type(&dir,point_op(op))==point_op_type) nb++;
  }
  return nb;
}  

PointGroupType find_point_group_type(const Array<rMatrix3d> &point_op) {
  PointGroupType point_group;
  rVector3d dir;
  switch (point_op.get_size()) {
  case 1:
    point_group=pt_grp_1;
    break;
  case 2:
    if (number_of_point_op(-1,point_op)>0) {
      point_group=pt_grp_1_;
    }
    else if (number_of_point_op(2,point_op)>0) {
      point_group=pt_grp_2;
    }
    else {
      point_group=pt_grp_m;
    }
    break;
  case 3:
    point_group=pt_grp_3;
    break;
  case 4:
    if (number_of_point_op(4,point_op)>0) {
      point_group=pt_grp_4;
    }
    else if (number_of_point_op(-4,point_op)>0) {
      point_group=pt_grp_4_;
    }
    else if (number_of_point_op(2,point_op)==3) {
      point_group=pt_grp_222;
    }
    else if (number_of_point_op(-2,point_op)==1) {
      point_group=pt_grp_2om;
    }
    else {
      point_group=pt_grp_2mm;
    }
    break;
  case 6:
    if (number_of_point_op(6,point_op)>0) {
      point_group=pt_grp_6;
    }
    else if (number_of_point_op(-6,point_op)>0) {
      point_group=pt_grp_6_;
    }
    else if (number_of_point_op(-1,point_op)==1) {
      point_group=pt_grp_3_;
    }
    else if (number_of_point_op(2,point_op)>0) {
      point_group=pt_grp_32;
    }
    else {
      point_group=pt_grp_3m;
    }
    break;
  case 8:
    if (number_of_point_op(-2,point_op)==0) {
      point_group=pt_grp_422;
    }
    else if (number_of_point_op(-2,point_op)==1) {
      point_group=pt_grp_4om;
    }
    else if (number_of_point_op(-2,point_op)==3) {
      point_group=pt_grp_mmm;
    }
    else if (number_of_point_op(4,point_op)>0) {
      point_group=pt_grp_4mm;
    }
    else {
      point_group=pt_grp_4_2m;
    }
    break;
  case 12:
    if (number_of_point_op(-2,point_op)==0) {
      if (number_of_point_op(6,point_op)>0) {
	point_group=pt_grp_622;
      }
      else {
	point_group=pt_grp_23;
      }
    }
    else if (number_of_point_op(6,point_op)>0) {
      if (number_of_point_op(-2,point_op)==1) {
	point_group=pt_grp_6om;
      }
      else {
	point_group=pt_grp_6mm;
      }
    }
    else {
      if (number_of_point_op(-6,point_op)>0) {
	point_group=pt_grp_6_m2;
      }
      else {
	point_group=pt_grp_3_m;
      }
    }
    break;
  case 16:
    point_group=pt_grp_4ommm;
    break;
  case 24:
    if (number_of_point_op(6,point_op)>0) {
      point_group=pt_grp_6ommm;
    }
    else if (number_of_point_op(-2,point_op)==0) {
      point_group=pt_grp_432;
    }
    else if (number_of_point_op(-4,point_op)>0) {
      point_group=pt_grp_4_3m;
    }
    else {
      point_group=pt_grp_m3_;
    }
    break;
  case 48:
    point_group=pt_grp_m3_m;
    break;
  default:
    ERRORQUIT("Error determining point group. Adjust tolerance.");
  }
  return point_group;
}

void find_atom_point_group(Array<rMatrix3d> *patom_point_group,
			   const rVector3d &atom_pos,
			   const rMatrix3d &cell, const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans) {
  rMatrix3d inv_cell=!cell;
  int nb_op=0;
  Array<int> op_index(point_op.get_size());
  for (int op=0; op<point_op.get_size(); op++) {
    if (equivalent_mod_cell(point_op(op)*atom_pos+trans(op),atom_pos,inv_cell)) {
      int op2;
      for (op2=0; op2<nb_op; op2++) {
	  if (max_norm(point_op(op)-point_op(op_index(op2)))<zero_tolerance) break;
      }
      if (op2==nb_op) {
       op_index(nb_op)=op;
	 nb_op++;
      }
    }
  }
  Array<rMatrix3d> &atom_point_group=*patom_point_group;
  atom_point_group.resize(nb_op);
  for (int op=0; op<nb_op; op++) {
    atom_point_group(op)=point_op(op_index(op));
  }
}

int contains_pure_translations(const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans) {
  rMatrix3d id;
  id.identity();
  int nb_id=0;
  for (int op=0; op<point_op.get_size(); op++) {
    if (max_norm(point_op(op)-id)<zero_tolerance) {nb_id++;}
  }
  return (nb_id>1);
}

/* BAD
int contains_pure_translations(const Structure &str) {
  Array<rMatrix3d> point_ops(1);
  Array<rVector3d> trans(1);
  point_ops(0).identity();
  trans(0)=rVector3d(0.,0.,0.);
  Structure t_str;
  t_str.cell=str.cell;
  t_str.atom_type=str.atom_type;
  for (int p=1; p<str.atom_pos.get_size() ; p++) {
    if (str.atom_type(p)==str.atom_type(0)) {
      rVector3d t=str.atom_pos(p)-str.atom_pos(0);
      apply_symmetry(&t_str.atom_pos,point_ops(0),t,str.atom_pos);
      if (equivalent_by_symmetry(t_str,str,str.cell,point_ops,trans)) {
	return 1;
      }
    }
  }
  return 0;
}
*/

int contains_pure_translations(const Structure &str, rMatrix3d unitcell) {
  rMatrix3d isupcel=!(str.cell);
  LatticePointInCellIterator l(unitcell,str.cell);
  for (; l; l++) {
    if (!near_zero(norm(l))) {
      int i;
      for (i=str.atom_pos.get_size()-1; i>=0; i--) {
	int j=which_atom(str.atom_pos,str.atom_pos(i)+(rVector3d)l,isupcel);
	if (str.atom_type(i)!=str.atom_type(j)) break;
      }
      if (i<0) return 1;
    }
  }
  return 0;
}

int contains_pure_translations_or_lexico_successor(const Structure &str, rMatrix3d unitcell) {
  rMatrix3d isupcel=!(str.cell);
  LatticePointInCellIterator l(unitcell,str.cell);
  for (; l; l++) {
    if (!near_zero(norm(l))) {
      int i;
      for (i=str.atom_pos.get_size()-1; i>=0; i--) {
	int j=which_atom(str.atom_pos,str.atom_pos(i)+(rVector3d)l,isupcel);
#ifdef HARTFORCADE
	if (str.atom_type(i)>str.atom_type(j)) return 1;
#endif
	if (str.atom_type(i)!=str.atom_type(j)) break;
      }
      if (i<0) return 1;
    }
  }
  return 0;
}

class SymmetryOp {
public:
  SymmetryOp(const rMatrix3d &_point_op, const rVector3d &_trans):point_op(_point_op),trans(_trans) {}
  rMatrix3d point_op;
  rVector3d trans;
};

void generate_space_group(SpaceGroup *p_fullgroup, const SpaceGroup &generator) {
  rMatrix3d invcell=!generator.cell;
  LinkedList<SymmetryOp> sym_list;
  for (int i=0; i<generator.point_op.get_size(); i++) {
    sym_list << new SymmetryOp(generator.point_op(i),generator.trans(i));
  }
  int done=0;
  while (!done) {
    done=1;
    LinkedListIterator<SymmetryOp> i(sym_list);
    for (; i; i++) {
      LinkedListIterator<SymmetryOp> j(sym_list);
      for (; j; j++) {
	MultiDimIterator<iVector3d> t(iVector3d(-2,-2,-2),iVector3d(2,2,2));
	for (; t; t++) {
	  // cerr << to_real(t) << endl;
	  SymmetryOp op(j->point_op * i->point_op,wrap_inside_cell(j->point_op * (i->trans + generator.cell*to_real(t)) + j->trans,generator.cell));
	  LinkedListIterator<SymmetryOp> k(sym_list);
	  for (; k; k++) {
	    //	  cerr << norm(op.point_op - k->point_op)/9 << " " << equivalent_mod_cell(op.trans,k->trans,invcell) << endl;
	    if (norm(op.point_op - k->point_op)/9<zero_tolerance && equivalent_mod_cell(op.trans,k->trans,invcell)) break;
	  }
	  if (!k) {
	    sym_list << new SymmetryOp(op.point_op,op.trans);
	    done=0;
	  }
	}
      }
    }
  }
  p_fullgroup->cell=generator.cell;
  p_fullgroup->point_op.resize(sym_list.get_size());
  p_fullgroup->trans.resize(sym_list.get_size());
  LinkedListIterator<SymmetryOp> i(sym_list);
  for (int j=0; i; i++,j++) {
    p_fullgroup->point_op(j)=i->point_op;
    p_fullgroup->trans(j)=i->trans;
  }
}
