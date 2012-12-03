#include "phonlib.h"
#include <fstream.h>

void remove_vacancies(Structure *p_str, const Array<AutoString> &label, const char *specie) {
  Structure str=*p_str;
  int newsize=str.atom_type.get_size();
  for (int at=0; at<str.atom_type.get_size(); at++) {
    if (strcmp(label(str.atom_type(at)),specie)==0) {newsize--;}
  }
  p_str->atom_pos.resize(newsize);
  p_str->atom_type.resize(newsize);
  int at=0;
  for (int i=0; i<str.atom_type.get_size(); i++) {
    if (strcmp(label(str.atom_type(i)),specie)!=0) {
      p_str->atom_pos(at)=str.atom_pos(i);
      p_str->atom_type(at)=str.atom_type(i);
      at++;
    }
  }
}

void add_rVector3d_to_ArrayReal(Array<Real> *pa, int at, const rVector3d &v) {
  for (int i=0; i<3; i++) {
    (*pa)(at*3+i)+=v(i);
  }
}

void make_nn_forcek(LinkedList<PairSpring> *p_pair, const Structure &str_ideal, const Structure &str_relax, const LinkedList<SpringSvsL> &spring_data) {
	rMatrix3d inv_cell=!str_ideal.cell;
	LinkedList<Array<rVector3d> > old_pair;
	for (int at=0; at<str_ideal.atom_pos.get_size(); at++) {
		Real r=find_1nn_radius(str_ideal,at)+zero_tolerance;
		AtomPairIterator pair(str_ideal.cell,str_ideal.atom_pos(at),str_ideal.atom_pos);
		while (get_length_quick(pair)<=r) {
			LinkedListIterator<Array<rVector3d> > old_p(old_pair);
			for (; old_p; old_p++) {
				if (equivalent_mod_cell(pair,*old_p,inv_cell)) break;
			}
			if (!old_p) {
				old_pair << new Array<rVector3d>(pair);
				PairSpring ps;
				ps.whichatom[0]=which_atom(str_ideal.atom_pos, pair(0), inv_cell);
				ps.whichatom[1]=which_atom(str_ideal.atom_pos, pair(1), inv_cell);
				rVector3d ideal_offset=(pair(1)-pair(0))-(str_ideal.atom_pos(ps.whichatom[1])-str_ideal.atom_pos(ps.whichatom[0]));
				ps.cell_shift=inv_cell*ideal_offset;
				ps.dr=str_relax.atom_pos(ps.whichatom[1])+str_relax.cell*ps.cell_shift-str_relax.atom_pos(ps.whichatom[0]);
				LinkedListIterator<SpringSvsL> si(spring_data);
				for (; si; si++) {
					if (   str_ideal.atom_type(ps.whichatom[0])==si->atom_type[0]
					    && str_ideal.atom_type(ps.whichatom[1])==si->atom_type[1]) break;
					if (   str_ideal.atom_type(ps.whichatom[0])==si->atom_type[1]
					    && str_ideal.atom_type(ps.whichatom[1])==si->atom_type[0]) break;
				}
				if (!si) {
				  cerr << "Did not find spring in make_nn_forcek. Ignoring" << endl;
				}
				else {
				  Real l=norm(ps.dr);
				  Real s=0.;
				  for (int i=0; i<si->forcek(0).get_size(); i++) {
				    s+=si->forcek(0)(i)*ipow(l,i);
				  }
				  Real b=0.;
				  for (int i=0; i<si->forcek(1).get_size(); i++) {
				    b+=si->forcek(1)(i)*ipow(l,i);
				  }
				  rMatrix3d diag;
				  diag.zero();
				  diag(0,0)=s;
				  diag(1,1)=b;
				  diag(2,2)=b;
				  rVector3d u=ps.dr;
				  u.normalize();
				  rVector3d pu=find_perpendicular(u);
				  rMatrix3d rot;
				  rot.set_column(0,u);
				  rot.set_column(1,pu);
				  rot.set_column(2,u^pu);
				  ps.forcek=rot*diag*(~rot);
				  (*p_pair) << new PairSpring(ps);
				}
			}
			pair++;
		}
	}
}

void make_nn_forcek(LinkedList<PairSpring> *p_pair, const Structure &str_ideal, const Structure &str_relax, const LinkedList<SpringSvsL> &spring_data, const LinkedList<rMatrix3d> &dirdep_mat, const Array<Real> &conc, const MultiDimPoly &multipoly) {
	rMatrix3d inv_cell=!str_ideal.cell;
	LinkedList<Array<rVector3d> > old_pair;
	for (int at=0; at<str_ideal.atom_pos.get_size(); at++) {
		Real r=find_1nn_radius(str_ideal,at)+zero_tolerance;
		AtomPairIterator pair(str_ideal.cell,str_ideal.atom_pos(at),str_ideal.atom_pos);
		while (get_length_quick(pair)<=r) {
			LinkedListIterator<Array<rVector3d> > old_p(old_pair);
			for (; old_p; old_p++) {
				if (equivalent_mod_cell(pair,*old_p,inv_cell)) break;
			}
			if (!old_p) {
				old_pair << new Array<rVector3d>(pair);
				PairSpring ps;
				ps.whichatom[0]=which_atom(str_ideal.atom_pos, pair(0), inv_cell);
				ps.whichatom[1]=which_atom(str_ideal.atom_pos, pair(1), inv_cell);
				rVector3d ideal_offset=(pair(1)-pair(0))-(str_ideal.atom_pos(ps.whichatom[1])-str_ideal.atom_pos(ps.whichatom[0]));
				ps.cell_shift=inv_cell*ideal_offset;
				ps.dr=str_relax.atom_pos(ps.whichatom[1])+str_relax.cell*ps.cell_shift-str_relax.atom_pos(ps.whichatom[0]);
				LinkedListIterator<SpringSvsL> si(spring_data);
				for (; si; si++) {
					if (   str_ideal.atom_type(ps.whichatom[0])==si->atom_type[0]
					    && str_ideal.atom_type(ps.whichatom[1])==si->atom_type[1]) break;
					if (   str_ideal.atom_type(ps.whichatom[0])==si->atom_type[1]
					    && str_ideal.atom_type(ps.whichatom[1])==si->atom_type[0]) break;
				}
				if (!si) {
				  cerr << "Did not find spring in make_nn_forcek. Ignoring" << endl;
				}
				else {
				  Real l=norm(ps.dr);
				  rVector3d ideal_dir=pair(1)-pair(0);
				  ideal_dir.normalize();
				  Array<Real> argpoly(multipoly.get_dim_arg());
				  argpoly(0)=l;
				  int ii=1;
				  LinkedListIterator<rMatrix3d> imat(dirdep_mat);
				  for (; imat ; imat++,ii++) {
				    argpoly(ii)=ideal_dir*((*imat)*ideal_dir);
				  }
				  for (int j=0; j<conc.get_size(); j++,ii++) {
				    argpoly(ii)=conc(j);
				  }
				  
				  Real s=multipoly.eval(si->forcek(0),argpoly);
				  Real b=multipoly.eval(si->forcek(1),argpoly);
				  rMatrix3d diag;
				  diag.zero();
				  diag(0,0)=s;
				  diag(1,1)=b;
				  diag(2,2)=b;
				  rVector3d u=ps.dr;
				  u.normalize();
				  rVector3d pu=find_perpendicular(u);
				  rMatrix3d rot;
				  rot.set_column(0,u);
				  rot.set_column(1,pu);
				  rot.set_column(2,u^pu);
				  ps.forcek=rot*diag*(~rot);
				  (*p_pair) << new PairSpring(ps);
				}
			}
			pair++;
		}
	}
}

void make_nn_forcek(LinkedList<PairSpring> *p_pair, const Structure &str_ideal, const Structure &str_relax, const LinkedList<SpringSvsL> &spring_data, const LinkedList<rMatrix3d> &dirdep_mat, const Array<Real> &conc, const MultiDimPoly &multipoly, LinkedList<Array<rVector3d> > *pold_pair) {
  rMatrix3d inv_cell=!str_ideal.cell;
  LinkedList<Array<rVector3d> > old_pair;
  if (!pold_pair) {
    pold_pair=&old_pair;
  }
  if (pold_pair->get_size()==0) {
    for (int at=0; at<str_ideal.atom_pos.get_size(); at++) {
      Real r=find_1nn_radius(str_ideal,at)+zero_tolerance;
      AtomPairIterator pair(str_ideal.cell,str_ideal.atom_pos(at),str_ideal.atom_pos);
      while (get_length_quick(pair)<=r) {
	LinkedListIterator<Array<rVector3d> > old_p(*pold_pair);
	for (; old_p; old_p++) {
	  if (equivalent_mod_cell(pair,*old_p,inv_cell)) break;
	}
	if (!old_p) {
	  (*pold_pair) << new Array<rVector3d>(pair);
	}
	pair++;
      }
    }
  }
  //cerr << pold_pair->get_size() << endl;
  LinkedListIterator<Array<rVector3d> > ipair(*pold_pair);
  while (ipair) {
    PairSpring ps;
    ps.whichatom[0]=which_atom(str_ideal.atom_pos, (*ipair)(0), inv_cell);
    ps.whichatom[1]=which_atom(str_ideal.atom_pos, (*ipair)(1), inv_cell);
    rVector3d ideal_offset=((*ipair)(1)-(*ipair)(0))-(str_ideal.atom_pos(ps.whichatom[1])-str_ideal.atom_pos(ps.whichatom[0]));
    ps.cell_shift=inv_cell*ideal_offset;
    ps.dr=str_relax.atom_pos(ps.whichatom[1])+str_relax.cell*ps.cell_shift-str_relax.atom_pos(ps.whichatom[0]);
    LinkedListIterator<SpringSvsL> si(spring_data);
    for (; si; si++) {
      if (   str_ideal.atom_type(ps.whichatom[0])==si->atom_type[0]
	     && str_ideal.atom_type(ps.whichatom[1])==si->atom_type[1]) break;
      if (   str_ideal.atom_type(ps.whichatom[0])==si->atom_type[1]
	     && str_ideal.atom_type(ps.whichatom[1])==si->atom_type[0]) break;
    }
    if (!si) {
      cerr << "Did not find spring in make_nn_forcek. Ignoring" << endl;
    }
    else {
      Real l=norm(ps.dr);
      rVector3d ideal_dir=(*ipair)(1)-(*ipair)(0);
      ideal_dir.normalize();
      Array<Real> argpoly(multipoly.get_dim_arg());
      argpoly(0)=l;
      int ii=1;
      LinkedListIterator<rMatrix3d> imat(dirdep_mat);
      for (; imat ; imat++,ii++) {
	argpoly(ii)=ideal_dir*((*imat)*ideal_dir);
      }
      for (int j=0; j<conc.get_size(); j++,ii++) {
	argpoly(ii)=conc(j);
      }
      
      Real s=multipoly.eval(si->forcek(0),argpoly);
      Real b=multipoly.eval(si->forcek(1),argpoly);
      rMatrix3d diag;
      diag.zero();
      diag(0,0)=s;
      diag(1,1)=b;
      diag(2,2)=b;
      rVector3d u=ps.dr;
      u.normalize();
      rVector3d pu=find_perpendicular(u);
      rMatrix3d rot;
      rot.set_column(0,u);
      rot.set_column(1,pu);
      rot.set_column(2,u^pu);
      ps.forcek=rot*diag*(~rot);
      (*p_pair) << new PairSpring(ps);
    }
    ipair++;
  }
}

// extendable version;
void make_nn_forcek(LinkedList<PairSpring> *p_pair, const Structure &str_ideal, const Structure &str_relax, const LinkedList<ExtendableSpringSvsL> &spring_data, const LinkedList<rMatrix3d> &dirdep_mat, const Array<Real> &conc) {
#ifdef DEBUG
  ofstream debugfile("fc_debug.out");
  debugfile << "at1 at2 length #concentration concentration(s) stretching bending" << endl;
#endif
	rMatrix3d inv_cell=!str_ideal.cell;
	LinkedList<Array<rVector3d> > old_pair;
	for (int at=0; at<str_ideal.atom_pos.get_size(); at++) {
		Real r=find_1nn_radius(str_ideal,at)+zero_tolerance;
		AtomPairIterator pair(str_ideal.cell,str_ideal.atom_pos(at),str_ideal.atom_pos);
		while (get_length_quick(pair)<=r) {
			LinkedListIterator<Array<rVector3d> > old_p(old_pair);
			for (; old_p; old_p++) {
				if (equivalent_mod_cell(pair,*old_p,inv_cell)) break;
			}
			if (!old_p) {
				old_pair << new Array<rVector3d>(pair);
				PairSpring ps;
				ps.whichatom[0]=which_atom(str_ideal.atom_pos, pair(0), inv_cell);
				ps.whichatom[1]=which_atom(str_ideal.atom_pos, pair(1), inv_cell);
				rVector3d ideal_offset=(pair(1)-pair(0))-(str_ideal.atom_pos(ps.whichatom[1])-str_ideal.atom_pos(ps.whichatom[0]));
				ps.cell_shift=inv_cell*ideal_offset;
				ps.dr=str_relax.atom_pos(ps.whichatom[1])+str_relax.cell*ps.cell_shift-str_relax.atom_pos(ps.whichatom[0]);
				LinkedListIterator<ExtendableSpringSvsL> si(spring_data);
				for (; si; si++) {
					if (   str_ideal.atom_type(ps.whichatom[0])==si->atom_type[0]
					    && str_ideal.atom_type(ps.whichatom[1])==si->atom_type[1]) break;
					if (   str_ideal.atom_type(ps.whichatom[0])==si->atom_type[1]
					    && str_ideal.atom_type(ps.whichatom[1])==si->atom_type[0]) break;
				}
				if (!si) {
				  cerr << "Did not find spring in make_nn_forcek. Ignoring" << endl;
				}
				else {
				  Real l=norm(ps.dr);
				  rVector3d ideal_dir=pair(1)-pair(0);
				  ideal_dir.normalize();
				  Array<Real> argfunc(1+dirdep_mat.get_size()+conc.get_size());
				  argfunc(0)=l;
				  int ii=1;
				  // first concentration then direction (opposite of standard SpringSvsL)
				  for (int j=0; j<conc.get_size(); j++,ii++) {
				    argfunc(ii)=conc(j);
				  }
				  LinkedListIterator<rMatrix3d> imat(dirdep_mat);
				  for (; imat ; imat++,ii++) {
				    argfunc(ii)=ideal_dir*((*imat)*ideal_dir);
				  }
				  Array<Real> sb(2);
				  si->pfunc->eval(&sb,argfunc);
				  rMatrix3d diag;
				  diag.zero();
				  diag(0,0)=sb(0);
				  diag(1,1)=sb(1);
				  diag(2,2)=sb(1);
#ifdef DEBUG
				  debugfile << si->atom_type[0] << " " << si->atom_type[1] << " " << l << " " << conc.get_size() << " ";
				  for (int j=0; j<conc.get_size(); j++) {debugfile << conc(j) << " ";}
				  debugfile << sb(0) << " " << sb(1) << endl;
#endif
				  rVector3d u=ps.dr;
				  u.normalize();
				  rVector3d pu=find_perpendicular(u);
				  rMatrix3d rot;
				  rot.set_column(0,u);
				  rot.set_column(1,pu);
				  rot.set_column(2,u^pu);
				  ps.forcek=rot*diag*(~rot);
				  (*p_pair) << new PairSpring(ps);
				}
			}
			pair++;
		}
	}
}

int get_basis_size(const LinkedList<PairSpring> &pair) {
	int n_basis=0;
	LinkedListIterator<PairSpring> p(pair);
	for (; p; p++) {
		n_basis=MAX(n_basis,MAX(p->whichatom[0],p->whichatom[1]));
	}
	return (n_basis+1);
}

void calc_dynmat(Array2d<Complex> *p_dynmat, const rVector3d &k, const LinkedList<PairSpring> &pair, Array<Real> mass) {
	const int dim=3;
      int n_basis=get_basis_size(pair);
	if (mass.get_size()==0) {
		mass.resize(n_basis);
		one_array(&mass);
	}
	p_dynmat->resize(iVector2d(dim*n_basis,dim*n_basis));
	zero_array(p_dynmat);
	LinkedListIterator<PairSpring> p(pair);
	for (; p; p++) {
		Real phase=2.*M_PI*k*(p->dr);
		rMatrix3d rblock=(p->forcek)*cos(phase);
		rMatrix3d iblock=(p->forcek)*sin(phase);
		int at0=p->whichatom[0];
		int at1=p->whichatom[1];
		for (int i=0; i<dim; i++) {
			for (int j=0; j<dim; j++) {
				(*p_dynmat)(at0*dim+i,at1*dim+j)-=Complex(rblock(i,j),iblock(i,j))/sqrt(mass(at0)*mass(at1));
				(*p_dynmat)(at1*dim+j,at0*dim+i)-=Complex(rblock(i,j),-iblock(i,j))/sqrt(mass(at0)*mass(at1));
				(*p_dynmat)(at0*dim+j,at0*dim+i)+=Complex(p->forcek(i,j),0)/sqrt(mass(at0)*mass(at0));
				(*p_dynmat)(at1*dim+i,at1*dim+j)+=Complex(p->forcek(i,j),0)/sqrt(mass(at1)*mass(at1));
			}
		}
	}
}

int list_phonon_freq(LinkedList<Real> *p_freq, LinkedList<Real> *p_weight, const iVector3d &kmesh, const rVector3d &kshift, const LinkedList<PairSpring> &pair, const rMatrix3d &cell, const Array<Real> &mass, Real convfact) {
	rMatrix3d kcell=!((~cell));
	for (int i=0; i<3; i++) {
	  kcell.set_column(i,kcell.get_column(i)/(Real)kmesh(i));
	}
	Real w=1./get_basis_size(pair)/(kmesh(0)*kmesh(1)*kmesh(2));
	Real minv=MAXFLOAT;
	Real maxv=0.;
	MultiDimIterator<iVector3d> ik(kmesh);
	for (; ik; ik++) {
		Array2d<Complex> dynmat;
		calc_dynmat(&dynmat, kcell*(kshift+to_real(ik)),pair,mass);
//cerr.precision(6);
//cerr.setf(ios::fixed);
//cerr << dynmat << endl;
		Array<Real> lambda;
		diagonalize_symmetric_matrix(&lambda,NULL,dynmat);
//cerr << lambda  << endl;
		for (int i=0; i<lambda.get_size(); i++) {
			Real v=sqrt(convfact*fabs(lambda(i)))/(2.*M_PI);
			if (lambda(i)<0) {v=-v;}
			(*p_freq) << new Real(v);
			if (p_weight) (*p_weight) << new Real(w);
			maxv=MAX(v,maxv);
			minv=MIN(v,minv);
		}
	}
	return (minv < -fabs(maxv)*1e-3 ? 0 : 1);
}

void calc_normal_modes(Array<Real> *p_freq, Array2d<Complex> *p_eigvect, const rVector3d &k, const LinkedList<PairSpring> &pair, const Array<Real> &mass, Real convfact) {
  Array2d<Complex> dynmat;
  calc_dynmat(&dynmat, k,pair,mass);
  Array<Real> lambda;
  Array2d<Complex> vect;
  if (p_eigvect) {
    diagonalize_symmetric_matrix(&lambda,&vect,dynmat);
  }
  else {
    diagonalize_symmetric_matrix(&lambda,NULL,dynmat);
  }
  int n=lambda.get_size();
  Array<int> lookup(n);
  for (int i=0; i<n; i++) {lookup(i)=i;}
  for (int i=1; i<n; i++) {
    for (int j=0; j<n-i; j++) {
      if (lambda(lookup(j+1))<lambda(lookup(j))) {
	swap(&(lookup(j+1)),&(lookup(j)));
      }
    }
  }
  p_freq->resize(lambda.get_size());
  for (int i=0; i<lambda.get_size(); i++) {
    (*p_freq)(i)=sqrt(convfact*fabs(lambda(lookup(i))))/(2.*M_PI);
    if (lambda(lookup(i))<0) {(*p_freq)(i)=-(*p_freq)(i);}
  }
  if (p_eigvect) {
    p_eigvect->resize(iVector2d(n,n));
    for (int i=0; i<n; i++) {
      Real norm=0.;
      for (int j=0; j<n; j++) {
	norm+=sqr(abs(vect(j,lookup(i))));
      }
      norm=sqrt(norm);
      for (int j=0; j<n; j++) {
	(*p_eigvect)(j,i)=vect(j,lookup(i))/norm;
      }
    }
  }
}

void calc_dispersion_curve(LinkedList<Array<Real> > *p_freq, LinkedList<Array2d<Complex> > *p_eigvect, const rVector3d &k1, const rVector3d &k2, int nstep, const LinkedList<PairSpring> &pair, const Array<Real> &mass, Real convfact) {
  rVector3d dk=(k2-k1)/(Real)MAX(1,nstep-1);
  rVector3d k=k1;
  for (int i=0; i<nstep; i++, k+=dk) {
    Array2d<Complex> dynmat;
    calc_dynmat(&dynmat, k,pair,mass);
    Array<Real> lambda;
    Array<Real> nu;
    Array2d<Complex> vect;
    if (p_eigvect) {
      diagonalize_symmetric_matrix(&lambda,&vect,dynmat);
    }
    else {
      diagonalize_symmetric_matrix(&lambda,NULL,dynmat);
    }
    int n=lambda.get_size();
    Array<int> lookup(n);
    for (int i=0; i<n; i++) {lookup(i)=i;}
    for (int i=1; i<n; i++) {
      for (int j=0; j<n-i; j++) {
	if (lambda(lookup(j+1))<lambda(lookup(j))) {
	  swap(&(lookup(j+1)),&(lookup(j)));
	}
      }
    }
    nu.resize(lambda.get_size());
    for (int i=0; i<lambda.get_size(); i++) {
      nu(i)=sqrt(convfact*fabs(lambda(lookup(i))))/(2.*M_PI);
      if (lambda(lookup(i))<0) {nu(i)=-nu(i);}
    }
    (*p_freq) << new Array<Real>(nu);
    if (p_eigvect) {
      Array2d<Complex> *pa=new Array2d<Complex>(iVector2d(n,n));
      for (int i=0; i<n; i++) {
	Real norm=0.;
	for (int j=0; j<n; j++) {
	  norm+=sqr(abs(vect(j,lookup(i))));
	}
	norm=sqrt(norm);
	for (int j=0; j<n; j++) {
	  (*pa)(j,i)=vect(j,lookup(i))/norm;
	}
      }
      (*p_eigvect) << pa;
    }
  }
}

Real real_erf(Real x) {
  Real du=0.01;
  Real maxx=3.5;
  if (x>maxx) return 1.;
  Real a=2*fabs(x)/sqrt(M_PI);
  Real adu=a*du;
  Real b=-sgn(x)*sqr(x);
  Real y=0.;
  for (Real u=0; u<=1+du/2; u+=du) {
    y+=exp(b*sqr(u))*adu;
  }
  return (y);
}

Real ln_real_erf(Real x) {
  Real du=0.01;
  Real maxx=3.5;
  if (x>maxx) return 0.;
  Real a=2*fabs(x)/sqrt(M_PI);
  Real adu=a*du;
  Real b=-sgn(x)*sqr(x);
  Real y=0.;
  for (Real u=0; u<=1+du/2; u+=du) {
    y+=exp(b*(sqr(u)-1.))*adu;
  }
  return (b+log(y));
}

Real calc_vib_free_energy_robust(LinkedList<Real> freq, LinkedList<Real> weight, Real kBT, Real hplanck, Real length) {
	LinkedListIterator<Real> v(freq);
	LinkedListIterator<Real> w(weight);
	Real c=(2.*M_PI*length)/sqrt(2.*kBT);
	Real f=0.;
	for (; v; v++, w++) {
	  Real a=(2.*kBT)/(hplanck*fabs(*v));
	  if (a<zero_tolerance) {
	    f+=(*w)*hplanck*fabs(*v)/2.;
	  }
	  else if (a>(1./zero_tolerance)) {
	    f+=(*w)*kBT*log(2./a);
	  }
	  else {
	    f+=(*w)*kBT*log(2.*sinh(1./a));
	  }
	  //	  cerr << (*w)*kBT*log(2.*sinh(1./a)) << " " << -(*w)*kBT*log(real_erf(c*(*v))) << " " << real_erf(c*(*v)) << endl;
	  if (kBT>zero_tolerance) {
	    f+=-(*w)*kBT*ln_real_erf(c*(*v));
	  }
	  else {
	    if (*v<0) {f+=-(*w)*2.*sqr(M_PI*length*(*v));}
	  }
	}
	return f;
}

Real calc_vib_free_energy(LinkedList<Real> freq, LinkedList<Real> weight, Real kBT, Real hplanck) {
	LinkedListIterator<Real> v(freq);
	LinkedListIterator<Real> w(weight);
	Real f=0.;
	for (; v; v++, w++) {
	  if ((*v)>0) {
	    Real a=(2.*kBT)/(hplanck*(*v));
	    if (a<zero_tolerance) {
	      f+=(*w)*hplanck*(*v)/2.;
	    }
	    else if (a>(1./zero_tolerance)) {
  	      f+=(*w)*kBT*log(2./a);
	    }
	    else {
	      f+=(*w)*kBT*log(2.*sinh(1./a));
	    }
	  }
	}
	return f;
}

Real calc_vib_entropy(LinkedList<Real> freq, LinkedList<Real> weight, Real kBT, Real hplanck) {
	LinkedListIterator<Real> v(freq);
	LinkedListIterator<Real> w(weight);
	Real s=0.;
	Real tw=0.,ttw=0.;
	if (kBT<zero_tolerance) {
		return 0.;
	}
	if (kBT==MAXFLOAT) {
		for (; v; v++, w++) {
		  ttw+=*w;
		  if ((*v)>zero_tolerance) {
			s+=-(*w)*log((*v));
			tw+=*w;
		  }
		}
	}
	else {
		for (; v; v++, w++) {
		  if ((*v)>zero_tolerance) {
		    Real r=hplanck*(*v)/(2.*kBT);
		    s+=-(*w)*(log(2.*sinh(r))-r/tanh(r));
		  }
		}
	}
	return s*ttw/tw;
}

void calc_k_mesh(iVector3d *p_mesh, const rMatrix3d &cell, Real nbkpts) {
  if (nbkpts<=1.) {
    *p_mesh=iVector3d(1,1,1);
    return;
  }
  rMatrix3d reciprocal=~(!cell);
  rVector3d proj;
  for (int i=0; i<3; i++) {
    rVector3d a=reciprocal.get_column( i     );
    rVector3d b=reciprocal.get_column((i+1)%3);
    rVector3d c=reciprocal.get_column((i+2)%3);
    proj(i)=fabs(a*((b^c)/norm(b^c)));
  }
  Real normalizer=pow(nbkpts/(proj(0)*proj(1)*proj(2)),1./3.);
  rVector3d mesh=proj*normalizer;
  iVector3d imesh=to_int(mesh);
  rVector3d fmesh=mesh-to_real(imesh);
  while (imesh(0)*imesh(1)*imesh(2)<nbkpts) {
    Real bestf=0.;
    for (int i=0; i<3; i++) {
      if (fmesh(i)>bestf) {
        bestf=fmesh(i);
      }
    }
    for (int i=0; i<3; i++) {
      if (fabs(fmesh(i)-bestf)<zero_tolerance) {
        imesh(i)++;
        fmesh(i)=0.;
      }
    }
  }
  *p_mesh=imesh;
}

void reorder_atoms(Structure *rel_str, const Structure &str, Array<int> *pcopy_from, Array<iVector3d> *pshift) {
  {
    Real scale=pow(det(rel_str->cell)/det(str.cell),1./3.);
    iMatrix3d supercell=to_int((!(rel_str->cell))*(str.cell)*scale);
    rel_str->cell=(rel_str->cell)*to_real(supercell);
  }
  /*
  {
    Real maxr=0;
    for (int i=0; i<3; i++) {maxr=MAX(maxr,norm(str.cell.get_column(i)));}
    Real scale=pow(det(str.cell)/det(rel_str->cell),1./3.);
    maxr*=MAX(1.,pow(scale,-3.));
    Array<rMatrix3d> supercell;
    find_all_equivalent_cell(&supercell,rel_str->cell,maxr);
    Real mind=MAXFLOAT;
    int best_ss;
    for (int s=0; s<supercell.get_size(); s++) {
      Real d=norm(supercell(s)*scale-str.cell);
      if (d<mind) {
	mind=d;
	best_ss=s;
      }
    }
    rel_str->cell=supercell(best_ss);
  }
  */
  int n=str.atom_pos.get_size();
  Array<rVector3d> frac_str,frac_rel_str;
  Array<int> str_done(n),rel_str_done(n);
  zero_array(&str_done);
  zero_array(&rel_str_done);
  rVector3d t(0.,0.,0.);
  apply_symmetry(&frac_str,!str.cell,t,str.atom_pos);
  apply_symmetry(&frac_rel_str,!rel_str->cell,t,rel_str->atom_pos);
  Array<int> copy_from(n);
  Array<iVector3d> shift(n);
  for (int t=0; t<n; t++) {
    Real min_d=MAXFLOAT;
    int best_ir=-1;
    int best_i=-1;
    for (int ir=0; ir<n; ir++) {
      if (rel_str_done(ir)==0) {
        for (int i=0; i<n; i++) {
          if (str_done(i)==0) {
            Real d=cylinder_norm(frac_rel_str(ir)-frac_str(i));
            if (d<min_d) {
              min_d=d;
              best_ir=ir;
              best_i=i;
            }
          }
        }
      }
    }
    if (str.atom_type(best_i)!=rel_str->atom_type(best_ir)) {
      ERRORQUIT("Too large shift between relaxed and unrelaxed positions.");
    }
    copy_from(best_i)=best_ir;
    rVector3d s=frac_str(best_i)-frac_rel_str(best_ir);
    shift(best_i)=to_int(s);
    rel_str_done(best_ir)=1;
    str_done(best_i)=1;
  }
  Structure tmpstr(*rel_str);
  for (int i=0; i<n; i++) {
    rel_str->atom_pos(i)=tmpstr.atom_pos(copy_from(i))+rel_str->cell*to_real(shift(i));
    rel_str->atom_type(i)=tmpstr.atom_type(copy_from(i));
  }
  if (pcopy_from) (*pcopy_from)=copy_from;
  if (pshift) (*pshift)=shift;
}

Real calc_bulk_modulus(const LinkedList<PairSpring> &pair, Real volume) {
  LinkedListIterator<PairSpring> ip(pair);
  Real b=0.;
  for (; ip; ip++) {
//cerr << ip->dr << endl;
//cerr << ip->forcek << endl;
//cerr << (ip->dr)*((ip->forcek)*(ip->dr)) << endl;
    b+=(ip->dr)*((ip->forcek)*(ip->dr));
  }
  return (b/volume/9);
}

void find_special_direction(Array<rVector3d> *pspecial_dir, Array<int> *pminus_equiv,
                      const Array<rMatrix3d> &point_op) {
  PointGroupType point_group_type=find_point_group_type(point_op);
  Array<rVector3d> &special_dir=*pspecial_dir;
  rVector3d dir,pert_dir,pert_pert_dir;
  switch (point_group_type) {
  case pt_grp_1:
  case pt_grp_1_:
    special_dir.resize(3);
    special_dir(0)=rVector3d(1.,0.,0.);
    special_dir(1)=rVector3d(0.,1.,0.);
    special_dir(2)=rVector3d(0.,0.,1.);
    break;
  case pt_grp_3:
  case pt_grp_4:
  case pt_grp_6:
  case pt_grp_3_:
  case pt_grp_4_:
  case pt_grp_6_:
  case pt_grp_4om:
  case pt_grp_6om:
    special_dir.resize(1);
    for (int op=0; op<point_op.get_size(); op++) {
      int type=find_point_op_type(&dir,point_op(op));
      if (abs(type)>=3 && abs(type)<=6) break;
    }
    pert_dir=find_perpendicular(dir);
    special_dir(0)=dir+pert_dir;
    break;
  case pt_grp_32:
  case pt_grp_422:
  case pt_grp_622:
  case pt_grp_4ommm:
  case pt_grp_6ommm:
  case pt_grp_6_m2:
  case pt_grp_3_m:
  case pt_grp_4_2m:
    special_dir.resize(1);
    for (int op=0; op<point_op.get_size(); op++) {
      int type=find_point_op_type(&dir,point_op(op));
      if (abs(type)>=3 && abs(type)<=6) break;
    }
    for (int op=0; op<point_op.get_size(); op++) {
      int type=find_point_op_type(&pert_dir,point_op(op));
      if (type==2) {
	if (norm( pert_dir*(pert_dir*dir)-dir ) > zero_tolerance) break;
      }
    }
    special_dir(0)=dir+pert_dir;
    break;
  case pt_grp_3m:
  case pt_grp_4mm:
  case pt_grp_6mm:
    special_dir.resize(1);
    for (int op=0; op<point_op.get_size(); op++) {
      int type=find_point_op_type(&dir,point_op(op));
      if (type>=3 && type<=6) break;
    }
    for (int op=0; op<point_op.get_size(); op++) {
      int type=find_point_op_type(&pert_dir,point_op(op));
      if (type==-2) break;
    }
    pert_pert_dir=(dir^pert_dir);
    special_dir(0)=dir+pert_pert_dir;
    break;
  case pt_grp_2mm:
    special_dir.resize(1);
    for (int op=0; op<point_op.get_size(); op++) {
      int type=find_point_op_type(&dir,point_op(op));
      if (type==2) break;
    }
    for (int op=0; op<point_op.get_size(); op++) {
      int type=find_point_op_type(&pert_dir,point_op(op));
      if (type==-2) break;
    }
    pert_pert_dir=(dir^pert_dir);
    special_dir(0)=dir+pert_dir+pert_pert_dir;
    break;
  case pt_grp_2:
  case pt_grp_m:
  case pt_grp_2om:
    special_dir.resize(2);
    for (int op=0; op<point_op.get_size(); op++) {
      int type=find_point_op_type(&dir,point_op(op));
      if (abs(type)==2) break;
    }
    pert_dir=find_perpendicular(dir);
    special_dir(0)=dir+pert_dir;
    pert_pert_dir=(dir^pert_dir);
    special_dir(1)=pert_pert_dir;
    break;
  case pt_grp_222:
  case pt_grp_mmm:
    {
      special_dir.resize(1);
      special_dir(0)=rVector3d(0.,0.,0.);
      for (int op=0; op<point_op.get_size(); op++) {
	int type=find_point_op_type(&dir,point_op(op));
	if (type==2) {
	  special_dir(0)+=dir;
	}
      }
    }
    break;
  case pt_grp_23:
    special_dir.resize(1);
    for (int op=0; op<point_op.get_size(); op++) {
      int type=find_point_op_type(&dir,point_op(op));
      if (type==3) break;
    }
    special_dir(0)=dir;
    break;
  case pt_grp_m3_:
    special_dir.resize(1);
    for (int op=0; op<point_op.get_size(); op++) {
      int type=find_point_op_type(&dir,point_op(op));
      if (type==-3) break;
    }
    special_dir(0)=dir;
    break;
  case pt_grp_432:
  case pt_grp_m3_m:
    special_dir.resize(1);
    for (int op=0; op<point_op.get_size(); op++) {
      int type=find_point_op_type(&dir,point_op(op));
      if (type==4) break;
    }
    special_dir(0)=dir;
    break;
  case pt_grp_4_3m:
    special_dir.resize(1);
    for (int op=0; op<point_op.get_size(); op++) {
      int type=find_point_op_type(&dir,point_op(op));
      if (type==-4) break;
    }
    special_dir(0)=dir;
    break;
  }
  for (int d=0; d<special_dir.get_size(); d++) {
    special_dir(d).normalize();
  }
  if (pminus_equiv) {
    pminus_equiv->resize(special_dir.get_size());
    for (int d=0; d<special_dir.get_size(); d++) {
      int sym;
      for (sym=0; sym<point_op.get_size(); sym++) {
        if (norm(point_op(sym)*special_dir(d)+special_dir(d))<zero_tolerance) break;
      }
      if (sym==point_op.get_size()) {
        (*pminus_equiv)(d)=0;
      }
      else {
        (*pminus_equiv)(d)=1;
      }
    }
  }
}

void find_pertubations(Array<int> *pwhich_atom, Array<rVector3d> *pdirection, Array<int> *pdo_minus,
		       const Structure &str,
		       const SpaceGroup &spacegroup) {
  LinkedList<int> which_atom_list;
  LinkedList<rVector3d> dir_list;
  LinkedList<int> minus_list;
  LinkedList<Array<rVector3d> > pointlist;
  for (int at=0; at<str.atom_pos.get_size(); at++) {
    Array<rVector3d> p(1);
    p(0)=str.atom_pos(at);
    if (add_unique(&pointlist,p,spacegroup)) {
      Array<rMatrix3d> atom_point_group;
      find_atom_point_group(&atom_point_group,str.atom_pos(at),str.cell,spacegroup.point_op,spacegroup.trans);
      Array<rVector3d> special_dir;
      Array<int> minus_equiv;
      find_special_direction(&special_dir, &minus_equiv, atom_point_group);
      for (int i=0; i<special_dir.get_size(); i++) {
        which_atom_list << new int(at);
        dir_list << new rVector3d(special_dir(i));
        minus_list << new int(1-minus_equiv(i));
      }
    }
  }
  LinkedList_to_Array(pwhich_atom,which_atom_list);
  LinkedList_to_Array(pdirection,dir_list);
  if (pdo_minus) {
    LinkedList_to_Array(pdo_minus,minus_list);
  }
}


void stretch_str(Structure *pstr, const Structure &str, Real strain) {
  pstr->cell=str.cell*(1.+strain);
  pstr->atom_pos.resize(str.atom_pos.get_size());
  for (int at=0; at<str.atom_pos.get_size(); at++) {
    pstr->atom_pos(at)=str.atom_pos(at)*(1.+strain);
  }
  pstr->atom_type=str.atom_type;
}

void calc_pert_str(Structure *ppert_str, const rMatrix3d &supercell, const Structure &str, const rVector3d &pos, const rVector3d &displ, Real strain) {
  Real scale=1.+strain;
  find_all_atom_in_supercell(&(ppert_str->atom_pos), &(ppert_str->atom_type), str.atom_pos, str.atom_type, str.cell, supercell);
  if (norm(displ)!=0.) {
    int displ_atom=which_atom(ppert_str->atom_pos,pos,!supercell);
    ppert_str->atom_pos(displ_atom)+=displ/scale;
  }
  ppert_str->cell=supercell*scale;
  for (int at=0; at<ppert_str->atom_pos.get_size(); at++) {
    ppert_str->atom_pos(at)=ppert_str->atom_pos(at)*scale;
  }

}

void read_force_vector(Array<Real> *pforce, istream &file, const Array<int> &copyfrom) {
  Array<Real> force(3*copyfrom.get_size());
  for (int at=0; at<3*copyfrom.get_size(); at++) {
    file >> force(at);
  }
  pforce->resize(force.get_size());
  for (int at=0; at<copyfrom.get_size(); at++) {
    for (int i=0; i<3; i++) {
      (*pforce)(3*at+i)=force(3*copyfrom(at)+i);
    }
  }
}

Real smooth_kernel(Real x) {
  return exp(-x*x/2)/sqrt(2.*M_PI);
}

void smooth_density(Real *xmin, Real *xmax, Array<Real> *pf, const Array<Real> &x) {
  *xmin=MIN(min(x),0.);
  *xmax=max(x);
  Real h=0.1*((*xmax)-(*xmin))*pow(x.get_size(),-1./5.);
  *xmax+=3.*h;
  int n=10*(int)(((*xmax)-(*xmin))/h);
  pf->resize(n);
  zero_array(pf);
  for (int i=0; i<n; i++) {
    Real y=*xmin+((*xmax)-(*xmin))*(Real)i/(Real)(n-1);
    for (int j=0; j<x.get_size(); j++) {
      (*pf)(i)+=smooth_kernel((x(j)-y)/h)/h;
    }
  }
}

void smooth_density(Real *xmin, Real *xmax, Array<Real> *pf, const Array<Real> &x, const Array<Real> &w) {
  *xmin=MIN(min(x),0.);
  *xmax=max(x);
  Real h=0.1*((*xmax)-(*xmin))*pow(x.get_size(),-1./5.);
  *xmax+=3.*h;
  int n=10*(int)(((*xmax)-(*xmin))/h);
  pf->resize(n);
  zero_array(pf);
  for (int i=0; i<n; i++) {
    Real y=*xmin+((*xmax)-(*xmin))*(Real)i/(Real)(n-1);
    for (int j=0; j<x.get_size(); j++) {
      (*pf)(i)+=w(j)*smooth_kernel((x(j)-y)/h)/h;
    }
  }
}

void find_sym_const_tensor(LinkedList<rMatrix3d> *pfc_list, const Array<rMatrix3d> &pointgroup, int forcesym) {
  Array<Array<Real> > basis;
  int nb=0;
  for (int ifc=0; ifc<3; ifc++) {
    for (int jfc=0; jfc<3; jfc++) {
      rMatrix3d fc,acc;
      fc.zero();
      fc(ifc,jfc)=1.;
      acc.zero();
      for (int op=0; op< pointgroup.get_size(); op++) {
	rMatrix3d tfc=pointgroup(op)*fc*(~pointgroup(op));
	acc=acc+tfc;
      }
      if (forcesym) {
	acc=(acc+~acc)/2;
      }
      Real l=max_norm(acc);
      if (l>zero_tolerance) {
	acc=acc/l;
	Array<Real> new_comp(9);
	int k=0;
	for (int i=0; i<3; i++) {
	  for (int j=0; j<3; j++) {
	    new_comp(k)=acc(i,j);
	    k++;
	  }
	}
	if (build_basis(&basis,&nb,new_comp)) {
	  (*pfc_list) << new rMatrix3d(acc);
	}
      }
    }
  }
}

