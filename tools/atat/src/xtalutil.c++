#include "xtalutil.h"

int in01(rVector3d v) {
  for (int i=0; i<3; i++) {
    if (v(i)<0. || v(i)>=1.) return 0;
  }
  return 1;
}

int in01included(rVector3d v) {
  for (int i=0; i<3; i++) {
    if (v(i)<=-zero_tolerance || v(i)>=1+zero_tolerance) return 0;
  }
  return 1;
}

int which_atom(const Array<rVector3d> &atom_pos, const rVector3d &pos, const rMatrix3d &inv_cell) {
  for (int i=0; i<atom_pos.get_size(); i++) {
    if (equivalent_mod_cell(atom_pos(i),pos,inv_cell)) return i;
  }
  return -1;
}

int which_atom(const Array<rVector3d> &atom_pos, const rVector3d &pos) {
  for (int i=0; i<atom_pos.get_size(); i++) {
    if (norm(atom_pos(i)-pos)<zero_tolerance) return i;
  }
  return -1;
}

rVector3d wrap_inside_cell(const rVector3d &v, const rMatrix3d &cell) {
  rVector3d shift=rVector3d(1.,1.,1.)*M_PI*zero_tolerance*0.1;
  return (cell*(mod1((!cell)*v-shift)+shift));
}

void wrap_inside_cell(Array<rVector3d> *dest, const Array<rVector3d> &src, const rMatrix3d &cell) {
  if (dest->get_size()!=src.get_size()) {dest->resize(src.get_size());}
  rVector3d shift=rVector3d(1,1,1)*M_PI*zero_tolerance*0.1;
  rMatrix3d inv_cell=!cell;
  for (int i=0; i<src.get_size(); i++) {
    (*dest)(i)=cell*(mod1(inv_cell*src(i)-shift)+shift);
    //(*dest)(i)=cell*(mod1(inv_cell*src(i)));
  }
}

rVector3d flip_into_brillouin_1(const rVector3d &k, const rMatrix3d &rec_lat) {
  iVector3d max_cell=find_sphere_bounding_box(rec_lat,2.*norm(k));
  MultiDimIterator<iVector3d> c(-max_cell,max_cell);
  rVector3d closest_k;
  Real min_norm=MAXFLOAT;
  for (;c; c++) {
    rVector3d cur_k=k+rec_lat*to_real(c);
    Real norm_k=norm(cur_k);
    if (norm_k<min_norm) {
      closest_k=cur_k;
      min_norm=norm_k;
    }
  }
  return closest_k;
}

rVector3d find_perpendicular(const rVector3d &dir) {
  rVector3d result(0.,0.,0.);
  for (int axis=0; axis<3; axis++) {
    rVector3d v(0.,0.,0.);
    v(axis)=1.;
    result=v-dir*(dir*v)/(dir*dir);
    Real l=norm(result);
    if (l>zero_tolerance) {
      result=result/l;
      break;
    }
  }
  return result;
}

iVector3d make_miller(rVector3d normal) {
  int p,q,common;
  common=1;
  for (int i=0; i<3; i++) {
    integer_ratio(&p,&q,normal(i),zero_tolerance);
    least_common_multiple(q,1);
  }
  iVector3d miller;
  for (int i=0; i<3; i++) {
    miller(i)=iround(normal(i)*(Real)common);
  }
  return miller;
}

void calc_lattice_vectors(rMatrix3d *lattice_vector, const rVector3d &cell_length, rVector3d cell_angle) {
  cell_angle=cell_angle*M_PI/180.;
  rVector3d a=rVector3d(1.                , 0.                , 0.)*cell_length(0);
  rVector3d b=rVector3d(cos(cell_angle(2)), sin(cell_angle(2)), 0.)*cell_length(1);
  rVector3d axb=a^b;
  rMatrix3d linear_system;
  linear_system.set_row(0,a);
  linear_system.set_row(1,b);
  linear_system.set_row(2,axb);
  rVector3d right_hand_side;
  right_hand_side(0)=norm(a)*cell_length(2)*cos(cell_angle(1));
  right_hand_side(1)=norm(b)*cell_length(2)*cos(cell_angle(0));
  right_hand_side(2)=0.;
  rVector3d solution=(!linear_system)*right_hand_side;
  Real x=sqrt(sqr(cell_length(2))-norm2(solution))/norm(axb);
  lattice_vector->set_column(0,a);
  lattice_vector->set_column(1,b);
  lattice_vector->set_column(2,solution+x*axb);
}

void calc_lattice_vectors(rVector3d *cell_length, rVector3d *cell_angle, const rMatrix3d &lattice_vector) {
  for (int i=0; i<3; i++) {
    (*cell_length)(i)=norm(lattice_vector.get_column(i));
  }
  for (int i=0; i<3; i++) {
    (*cell_angle)(i)=acos(lattice_vector.get_column((i+1)%3)*lattice_vector.get_column((i+2)%3)/(((*cell_length)((i+1)%3))*((*cell_length)((i+2)%3))))*180.0/M_PI;
  }
}

iVector3d find_sphere_bounding_box(const rMatrix3d &cell, Real radius) {
  iVector3d max_cell;
  rMatrix3d rec_lat=~!cell;
  for (int i=0 ; i<3; i++) {
    max_cell(i)=(int)ceil(norm(rec_lat.get_column(i))*radius+zero_tolerance);
  }
  return(max_cell);
}

int equivalent_mod_cell(const rVector3d &a, const rVector3d &b, const rMatrix3d &inv_cell) {
  rVector3d frac_a=inv_cell*a;
  rVector3d frac_b=inv_cell*b;
  rVector3d delta=frac_a-frac_b;
  for (int i=0; i<3; i++) {
    delta(i)=cylinder_norm(delta(i));
  }
  return (norm(delta)<zero_tolerance ? 1 : 0);
}

int equivalent_mod_cell(const Array<rVector3d> &a, const Array<rVector3d> &b, const rMatrix3d &inv_cell) {
  if (a.get_size()!=b.get_size()) return 0;
  if (a.get_size()==0) return 1;
  for (int i=0; i<a.get_size(); i++) {
    if (equivalent_mod_cell(a(0),b(i),inv_cell)) {
      rVector3d shift=b(i)-a(0);
      int j;
      for (j=0; j<a.get_size(); j++) {
        if (which_atom(b,a(j)+shift)==-1) break; // if did not find match
      }
      if (j==a.get_size()) return 1; // if all match
    }
  }
  return 0;
}

int equivalent_mod_cell(const MultiCluster &a, const MultiCluster &b, const rMatrix3d &inv_cell) {
  if (a.clus.get_size()!=b.clus.get_size()) return 0;
  if (a.clus.get_size()==0) return 1;
  for (int i=0; i<a.clus.get_size(); i++) {
    if (equivalent_mod_cell(a.clus(0),b.clus(i),inv_cell)) {
      rVector3d shift=b.clus(i)-a.clus(0);
      int j;
      for (j=0; j<a.clus.get_size(); j++) {
	  int at=which_atom(b.clus,a.clus(j)+shift);
        if (at==-1) break; // if did not find matching site
	  if (b.func(at)!=a.func(j)) break; // if did not find matching func
      }
      if (j==a.clus.get_size()) return 1; // if all match
    }
  }
  return 0;
}

void LatticePointIterator::find_new_neighbors(Real new_max_radius) {
  iVector3d max_cell=find_sphere_bounding_box(cell,new_max_radius);
  MultiDimIterator<iVector3d> c(-max_cell,max_cell);
  for (;c;c++) {
    rVector3d x=cell*to_real(c);
    Real r=norm(x);
    if (max_radius<r && r<=new_max_radius) {
      LinkedListIterator<rVector3d> i(list);
      while (i && norm(*i)<r) i++;
      list.add(new rVector3d(x),i);
    }
  }
  max_radius=new_max_radius;
}

void LatticePointIterator::init(const rMatrix3d &_cell, int skipzero) {
  list.delete_all();
  cell=_cell;
  inc_radius=0;
  for (int i=0; i<3; i++) {
    inc_radius=MAX(inc_radius,norm(cell.get_column(i)));
  }
  inc_radius+=zero_tolerance;
  max_radius=(skipzero ? 0. : -1.);
  find_new_neighbors(inc_radius);
}

rVector3d LatticePointIterator::operator++(int) {
  LinkedListIterator<rVector3d> i(list);
  delete list.detach(i);
  i.init(list);
  while (!i) {
    find_new_neighbors(max_radius+inc_radius);
    i.init(list);
  }
  return *i;
}

void LatticePointInCellIterator::init(const rMatrix3d &_cell, const rMatrix3d &_supercell) {
  cell=_cell;
  cell_to_super=(!_supercell)*_cell;
  shift=rVector3d(1.,1.,1.)*zero_tolerance*0.1*M_PI; //some irrational number;

  rMatrix3d super_to_cell=!cell_to_super;
  BoundingBox<Real,3> box;
  box.before_update();
  MultiDimIterator<iVector3d> take_tip(iVector3d(0,0,0),iVector3d(1,1,1));
  for (; take_tip; take_tip++) {
    box.update(super_to_cell*(to_real(take_tip)+shift));
  }
  current.init(to_int(floor(box.min)),to_int(ceil(box.max)));
  advance_to_valid();
}

void LatticePointInCellIterator::advance_to_valid(void) {
  while (current) {
    if (in01(cell_to_super*to_real(current)+shift)) break;
    current++;
  }
}

void find_all_atom_in_supercell(Array<rVector3d> *psuper_atom_pos,
				Array<int> *psuper_atom_type,
				const Array<rVector3d> &atom_pos,
				const Array<int> &atom_type,
				const rMatrix3d &cell, const rMatrix3d &supercell) {
  int nb_atom_supercell=iround(atom_pos.get_size()*fabs(det(supercell)/det(cell)));
  psuper_atom_pos->resize(nb_atom_supercell);
  psuper_atom_type->resize(nb_atom_supercell);
  LatticePointInCellIterator l(cell,supercell);
  int index=0;
  for ( ; l; l++) {
    for (int s=0; s<atom_pos.get_size(); s++) {
      (*psuper_atom_pos)(index)=l+atom_pos(s);
      (*psuper_atom_type)(index)=atom_type(s);
      index++;
    }
  }
  wrap_inside_cell(psuper_atom_pos,*psuper_atom_pos,supercell);
}

/*
void find_all_atom_in_supercell(Array<rVector3d> *psuper_atom_pos,
				Array<int> *psuper_atom_type,
				const Array<rVector3d> &atom_pos,
				const Array<int> &atom_type,
				const rMatrix3d &cell, const rMatrix3d &supercell) {
  rMatrix3d super_to_cell=(!cell)*supercell;
  rVector3d take_tip;
  BoundingBox<Real,3> box;
  box.before_update();
  for (take_tip(0)=0.; take_tip(0)<=1.; take_tip(0)+=1.) {
    for (take_tip(1)=0.; take_tip(1)<=1.; take_tip(1)+=1.) {
      for (take_tip(2)=0.; take_tip(2)<=1.; take_tip(2)+=1.) {
	box.update(super_to_cell*take_tip);
      }
    }
  }

  int nb_atom=atom_pos.get_size();
  int nb_atom_supercell=nb_atom*iround(fabs(det(supercell)/det(cell)));
  rMatrix3d axes_to_super=!supercell;

  Array<rVector3d> &super_atom_pos=*psuper_atom_pos;
  Array<int> &super_atom_type=*psuper_atom_type;
  super_atom_pos.resize(nb_atom_supercell);
  if (psuper_atom_type) super_atom_type.resize(nb_atom_supercell);
  {
    int index=0;
    rVector3d which_cell;
    for (which_cell(0)=box.min(0); which_cell(0)<=box.max(0)+zero_tolerance; which_cell(0)+=1.) {
      for (which_cell(1)=box.min(1); which_cell(1)<=box.max(1)+zero_tolerance; which_cell(1)+=1.) {
	for (which_cell(2)=box.min(2); which_cell(2)<=box.max(2)+zero_tolerance; which_cell(2)+=1.) {
	  for (int at=0; at<nb_atom; at++) {
	    rVector3d pos_axes=cell*which_cell+atom_pos(at);
	    rVector3d pos_super=axes_to_super*pos_axes;
	    int i;
	    for (i=0; i<3; i++) {
	      if (pos_super(i) < -zero_tolerance || pos_super(i) > 1.+zero_tolerance) break;
	    }
	    if (i==3) {
	      int index2;
	      for (index2=0; index2<index; index2++) {
		if (equivalent_mod_cell(super_atom_pos(index2),pos_axes,axes_to_super)) break;
	      }
	      if (index2==index) {
		for (i=0; i<3; i++) {
		  if (pos_super(i)>1.-zero_tolerance) pos_super(i)-=1.;
		}
		pos_axes=supercell*pos_super;
		super_atom_pos(index)=pos_axes;
		if (psuper_atom_type) super_atom_type(index)=atom_type(at);
		index++;
	      }
	    }
	  }
	}
      }
    }
    if (index!=super_atom_pos.get_size()) ERRORQUIT("Computed supercell does not have the right number of atoms.");
  }
}
*/

Real find_1nn_radius(const Structure &str, int atom) {
    Real dist=0.;
    LinkedList<ArrayrVector3d> pair_list;
    AtomPairIterator pair(str.cell,str.atom_pos(atom),str.atom_pos);
    while (1) {
      LinkedListIterator<ArrayrVector3d> i(pair_list);
      for ( ; i; i++) {
        if (cos_angle((*i)(1)-(*i)(0),pair(1)-(*i)(1))>-zero_tolerance) break;
      }
      if (i) {
        break;
      }
      else {
        pair_list << new Array<rVector3d>(pair);
        dist=norm(pair(1)-pair(0));
      }
      pair++;
    }
  return dist;
}

Real find_1nn_radius(const Structure &str) {
  Real max_dist=0.;
  for (int at=0; at<str.atom_pos.get_size(); at++) {
    Real dist=0.;
    LinkedList<ArrayrVector3d> pair_list;
    AtomPairIterator pair(str.cell,str.atom_pos(at),str.atom_pos);
    while (1) {
      LinkedListIterator<ArrayrVector3d> i(pair_list);
      for ( ; i; i++) {
        if (cos_angle((*i)(1)-(*i)(0),pair(1)-(*i)(1))>-zero_tolerance) break;
      }
      if (i) {
        break;
      }
      else {
        pair_list << new Array<rVector3d>(pair);
        dist=norm(pair(1)-pair(0));
      }
      pair++;
    }
    max_dist=MAX(max_dist,dist);
  }
  return max_dist;
}

void AtomPairIterator::find_new_neighbors(Real new_max_radius) {
  LatticePointIterator lat_vect(cell,0);
  for (; norm(lat_vect) < (new_max_radius+max_dist_in_cell+zero_tolerance); lat_vect++) {
    for (int i=0; i<atom_pos1->get_size(); i++) {
      for (int j=(same_pos && norm(lat_vect)<zero_tolerance ? i : 0); j<atom_pos2->get_size(); j++) {
        rVector3d atom_a=(*atom_pos1)(i);
        rVector3d atom_b=(*atom_pos2)(j)+lat_vect;
        Real r=norm(atom_b-atom_a);
        if (max_radius<r && r<=new_max_radius && r>zero_tolerance) {
          LinkedListIterator<ArrayrVector3d> i(list);
          while (i && norm((*i)(1)-(*i)(0))<r) i++;
	  Array<rVector3d> *parray=new Array<rVector3d>(2);
          (*parray)(0)=atom_a;
          (*parray)(1)=atom_b;
          list.add(parray,i);
        }
      }
    }
  }
  max_radius=new_max_radius;
}

void AtomPairIterator::init(void) {
  max_dist_in_cell=0.;
  for (int i=0; i<atom_pos1->get_size(); i++) {
    for (int j=0; j<atom_pos2->get_size(); j++) {
      max_dist_in_cell=MAX(max_dist_in_cell,norm((*atom_pos1)(i)-(*atom_pos2)(j)));
    }
  }
  inc_radius=0;
  for (int i=0; i<3; i++) {
    inc_radius=MAX(inc_radius,norm(cell.get_column(i)));
  }
  inc_radius+=zero_tolerance;
  max_radius=0.;
  LinkedListIterator<ArrayrVector3d> it(list);
  while (it) {delete list.detach(it);}
  find_new_neighbors(inc_radius);
}

void AtomPairIterator::init(const rMatrix3d &_cell, const Array<rVector3d> &_atom_pos) {
  cell=_cell;
  atom_pos1=&_atom_pos;
  atom_pos2=&_atom_pos;
  same_pos=1;
  init();
}

void AtomPairIterator::init(const rMatrix3d &_cell, const Array<rVector3d> &_atom_pos1, const Array<rVector3d> &_atom_pos2) {
  cell=_cell;
  atom_pos1=&_atom_pos1;
  atom_pos2=&_atom_pos2;
  same_pos=0;
  init();
}

void AtomPairIterator::init(const rMatrix3d &_cell, const rVector3d &_an_atom, const Array<rVector3d> &_atom_pos) {
  cell=_cell;
  an_atom.resize(1);
  an_atom(0)=_an_atom;
  atom_pos1=&an_atom;
  atom_pos2=&_atom_pos;
  same_pos=0;
  init();
}

void AtomPairIterator::operator++(int) {
  LinkedListIterator<ArrayrVector3d> i(list);
  delete list.detach(i);
  i.init(list);
  while (!i) {
    find_new_neighbors(max_radius+inc_radius);
    i.init(list);
  }
}

void AtomMultipletIterator::init(const rMatrix3d &_cell, const Array<rVector3d> &_atom_pos, int _ntuple) {
  cell=_cell;
  atom_pos=&_atom_pos;
  ntuple=_ntuple;
  multiplet.resize(ntuple);
  pairs.resize(ntuple-1);
  if (ntuple>1) {
    pairs(ntuple-2).init(cell,*atom_pos,*atom_pos);
    for (int i=0; i<ntuple-2; i++) {pairs(i).init(cell,pairs(ntuple-2)(0),*atom_pos);}
    (*this)++;
  }
}

void AtomMultipletIterator::operator++(int) {
  while (1) {
    int i=0;
    while (i<pairs.get_size()-1) {
      pairs(i)++;
      if (pairs(i).length()<=pairs(i+1).length()+zero_tolerance) {
        break;
      }
      pairs(i).init(cell,pairs(ntuple-2)(0),*atom_pos);
      i++;
    }
    if (i==pairs.get_size()-1) {
      pairs(i)++;
      for (int j=0; j<i; j++) {
	pairs(j).init(cell,pairs(ntuple-2)(0),*atom_pos);
      }
    }
    Real rmin=pairs(0).length();
    Real rmax=pairs(ntuple-2).length();
    for (i=0; i<pairs.get_size(); i++) {
      int j;
      for (j=i+1; j<pairs.get_size(); j++) {
        Real r=norm(pairs(i)(1)-pairs(j)(1));
        if (near_zero(r) || r<rmin-zero_tolerance || r>rmax+zero_tolerance) break;
      }
      if (j<pairs.get_size()) break;
    }
    if (i==pairs.get_size()) break;
  }
  multiplet(0)=pairs(ntuple-2)(0);
  for (int i=0; i<multiplet.get_size()-1; i++) {
    multiplet(i+1)=pairs(i)(1);
  }
}

Real get_cluster_length(const Array<rVector3d> &c) {
  Real max_length=0.;
  for (int i=0; i<c.get_size(); i++) {
    for (int j=i+1; j<c.get_size(); j++) {
      Real l=norm(c(i)-c(j));
      if (l>max_length) max_length=l;
    }
  }
  return max_length;
}

void find_common_supercell(rMatrix3d *psupercell, const rMatrix3d &a, const rMatrix3d &b) {
  rMatrix3d s,g;
  if (det(a)<det(b)) {
    s=a;
    g=b;
  }
  else {
    s=b;
    g=a;
  }
  rMatrix3d ig=!g;
  rVector3d v1,v2,v3;
  LatticePointIterator p(s);
  while (!is_int(ig*p)) p++;
  v1=p;
  while (!is_int(ig*p) || fabs(sqr(v1*p)-norm2(v1)*norm2(p))<zero_tolerance) p++;
  v2=p;
  while (!is_int(ig*p) || fabs((v1^v2)*p)<zero_tolerance) p++;
  v3=p;
  psupercell->set_column(0,v1);
  psupercell->set_column(1,v2);
  psupercell->set_column(2,v3);
}

void find_common_simple_supercell(iVector3d *psimple_supercell, const rMatrix3d &unitcell, const rMatrix3d &supercell) {
  rMatrix3d isupercell=!supercell;
  for (int i=0; i<3; i++) {
    int m=1;
    while (!is_int(isupercell*((Real)m*unitcell.get_column(i)))) m++;
    (*psimple_supercell)(i)=m;
  }
}

void strain_str(Structure *pstr, const Structure &str, const rMatrix3d &strain) {
  rMatrix3d id;
  id.identity();
  rMatrix3d s=strain+id;
  pstr->cell=s*(str.cell);
  pstr->atom_pos.resize(str.atom_pos.get_size());
  pstr->atom_type.resize(str.atom_pos.get_size());
  for (int at=0; at<str.atom_pos.get_size(); at++) {
    pstr->atom_pos(at)=s*str.atom_pos(at);
  }
  pstr->atom_type=str.atom_type;
}

