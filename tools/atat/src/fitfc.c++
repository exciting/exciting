#include <fstream.h>
#include <sys/stat.h>
#include "phonlib.h"
#include "getvalue.h"
#include "parse.h"
#include "version.h"
#include "lattype.h"
#include "arraylist.h"

template<class T>
class NormEpsilon {
  Real epsilon;
 public:
  NormEpsilon(Real _epsilon) {epsilon=_epsilon;}
  int operator () (const T& a, const T& b) const {
    return (norm(a-b)<epsilon);
  }
};

template<class T>
class NormArrayEpsilon {
  Real epsilon;
 public:
  NormArrayEpsilon(Real _epsilon) {epsilon=_epsilon;}
  int operator () (const Array<T>& a, const Array<T>& b) const {
    Real max=0.;
    for (int i=0; i<a.get_size(); i++) {
      max=MAX(max,norm(a(i)-b(i)));
    }
    return (max<epsilon);
  }
};

template<class T>
class NormPairEpsilon {
  Real epsilon;
 public:
  NormPairEpsilon(Real _epsilon) {epsilon=_epsilon;}
  int operator () (const Array<T>& a, const Array<T>& b) const {
    if (norm(a(0)-b(0))<epsilon && norm(a(1)-b(1))<epsilon) return 1;
    if (norm(a(0)-b(1))<epsilon && norm(a(1)-b(0))<epsilon) return 1;
    return 0;
  }
};

void find_allowed_forcek(LinkedList<rMatrix3d> *pfc_list, const Array<rVector3d> &pair, const SpaceGroup &spacegroup) {
  rMatrix3d inv_cell=!spacegroup.cell;
  Array<Array<Real> > basis;
  int nb=0;
  for (int ifc=0; ifc<3; ifc++) {
    for (int jfc=0; jfc<3; jfc++) {
      rMatrix3d fc,acc;
      fc.zero();
      fc(ifc,jfc)=1.;
      acc.zero();
      for (int op=0; op< spacegroup.point_op.get_size(); op++) {
        Array<rVector3d> tpair;
        apply_symmetry(&tpair,spacegroup.point_op(op),spacegroup.trans(op),pair);
        if (equivalent_mod_cell(tpair,pair,inv_cell)) {
          rMatrix3d tfc=spacegroup.point_op(op)*fc*(~spacegroup.point_op(op));
          if ((tpair(1)-tpair(0))*(pair(1)-pair(0))<0) {
            acc=acc+(~tfc);
          }
          else {
            acc=acc+(tfc);
          }
        }
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

class SpringGeneral {
 public:
  Array<rVector3d> pair;
  LinkedList<rMatrix3d> fc_list;
 public:
  SpringGeneral(void): pair(2), fc_list() {};
};

class EquivalentSpringGeneral {
  const SpaceGroup &spacegroup;
 public:
  EquivalentSpringGeneral(const SpaceGroup &_spacegroup): spacegroup(_spacegroup) {}
  int operator () (const SpringGeneral& a, const SpringGeneral& b) const {
    return spacegroup(a.pair,b.pair);
  }
};

void find_distinct_forcek(LinkedList<SpringGeneral> *pfc_list, Real radius, const Structure &str_ideal, const Structure &str, const SpaceGroup &spacegroup) {
  EquivalentSpringGeneral equiv(spacegroup);
  rMatrix3d inv_cell_ideal=!str_ideal.cell;
  AtomPairIterator pair(str_ideal.cell,str_ideal.atom_pos);
  while (get_length_quick(pair)<=radius+zero_tolerance) {
    Array<rVector3d> tmp=pair;
    int whichatom[2];
    whichatom[0]=which_atom(str_ideal.atom_pos, pair(0), inv_cell_ideal);
    whichatom[1]=which_atom(str_ideal.atom_pos, pair(1), inv_cell_ideal);
    rVector3d offset=str.cell*inv_cell_ideal*((pair(1)-pair(0))-(str_ideal.atom_pos(whichatom[1])-str_ideal.atom_pos(whichatom[0])));
     SpringGeneral spr;
     spr.pair(0)=str.atom_pos(whichatom[0]);
     spr.pair(1)=str.atom_pos(whichatom[1])+offset;
     add_unique(pfc_list,spr,equiv);
     pair++;
  }
  LinkedListIterator<SpringGeneral> s(*pfc_list);
  for ( ; s; s++) {
    //    cerr << (s->pair(0))/3.57635 << " " << ((s->pair)(1)-(s->pair)(0))/3.57635 << endl;
    find_allowed_forcek(&(s->fc_list),s->pair,spacegroup);
  }
}

void calc_force_vector(Array<Real> *forces, const Structure &unpert_str, const Structure &pert_str, const SpringGeneral &spring, int which_fc, SpaceGroup spacegroup) {
  forces->resize(unpert_str.atom_pos.get_size()*3);
  zero_array(forces);
  LinkedListIterator<rMatrix3d> fi(spring.fc_list);
  for (int i=0; i<which_fc; i++) {
    fi++;
  }
  NormPairEpsilon<rVector3d> identical(zero_tolerance);
  rMatrix3d fc=*fi;
  rMatrix3d inv_cell=!spacegroup.cell;
  rMatrix3d inv_bigcell=!unpert_str.cell;
  //  cerr << "begin" << endl;
  for (int at=0; at<pert_str.atom_pos.get_size(); at++) {
    LinkedList<Array<rVector3d> > pair_done;
    rVector3d d=pert_str.atom_pos(at)-unpert_str.atom_pos(at);
    if (norm(d)>zero_tolerance) {
      // cerr << pert_str.atom_pos(at) << unpert_str.atom_pos(at) << endl;
      for (int op=0; op<spacegroup.point_op.get_size(); op++) {
        rVector3d force(0.,0.,0.);
        Array<rVector3d> tpair;
        apply_symmetry(&tpair,spacegroup.point_op(op),spacegroup.trans(op),spring.pair);
	for (int e=0; e<2; e++) {
	  int ontop=equivalent_mod_cell(unpert_str.atom_pos(at),tpair(0),inv_cell);
	  if (ontop) {
	    tpair(1)+=unpert_str.atom_pos(at)-tpair(0);
	    tpair(0)=unpert_str.atom_pos(at);
	    if (add_unique(&pair_done,tpair,identical)) {
	      // cerr << tpair(0)/5.1631 << " " << (tpair(1)-tpair(0))/5.1631 << " " << which_fc << endl;
	      force=(spacegroup.point_op(op))*((~fc)*((~spacegroup.point_op(op))*d));
	      //                             here;
	      // cerr << spacegroup.point_op(op)*(fc*(~spacegroup.point_op(op))) << endl;
	      add_rVector3d_to_ArrayReal(forces,which_atom(unpert_str.atom_pos, tpair(0), inv_bigcell),-force);
	      add_rVector3d_to_ArrayReal(forces,which_atom(unpert_str.atom_pos, tpair(1), inv_bigcell),force);
	    }
	  }
	  fc=~fc;
	  swap(&(tpair(0)),&(tpair(1)));
        }
      }
    }
  }
}

Real calc_diff_energy(const Structure &unpert_str, const Structure &pert_str, const Array<Real> &force3) {
  Real sum=0.;
  for (int i=0; i<unpert_str.atom_pos.get_size(); i++) {
    rVector3d f(force3(3*i),force3(3*i+1),force3(3*i+2));
    sum+=(pert_str.atom_pos(i)-unpert_str.atom_pos(i))*f;
  }
  return -sum/2.;
}

void calc_equation_matrix(Array<Array<Real> > *peqn_mat, const Structure &unpert_str, const Structure &pert_str, const LinkedList<SpringGeneral> &spring_list, const SpaceGroup &spacegroup) {
  int nbcol=0;
  LinkedListIterator<SpringGeneral> i_s(spring_list);
  for (;i_s; i_s++) {
    nbcol+=i_s->fc_list.get_size();
  }
  int nbrow=pert_str.atom_pos.get_size()*3;
  resize(peqn_mat,nbrow,nbcol);
  int c=0;
  i_s.init(spring_list);
  for (;i_s; i_s++) {
    for (int ifc=0; ifc<i_s->fc_list.get_size(); ifc++) {
      Array<Real> force;
      calc_force_vector(&force, unpert_str,pert_str,*i_s,ifc,spacegroup);
      for (int i=0; i<nbrow; i++) {
        (*peqn_mat)(i)(c)=force(i);
      }
      c++;
    }
  }
}

class IdenticalPairSpring {
  Real epsilon;
 public:
  IdenticalPairSpring(Real _epsilon) {epsilon=_epsilon;}
  int operator () (const PairSpring& a, const PairSpring& b) const {
    if (a.whichatom[0]==b.whichatom[0] && a.whichatom[1]==b.whichatom[1] && norm(a.dr-b.dr)<epsilon) {return 1;}
    return 0;
  }
};


void print_fc(ostream &file, const LinkedList<SpringGeneral> &spring_data, const Array<Real> &beta, const Structure &str, const Array<AutoString> &label, const rMatrix3d &axes) {
  LinkedListIterator<SpringGeneral> i_sg(spring_data);
  int c=0;
  for (; i_sg; i_sg++) {
    rMatrix3d fc;
    fc.zero();
    LinkedListIterator<rMatrix3d> i_m(i_sg->fc_list);
    int bc=c;
    for ( ; i_m; i_m++, c++) {
      fc=fc+(*i_m)*beta(c);
    }
    file << "svsl ";
    file << label(str.atom_type(which_atom(str.atom_pos,i_sg->pair(0),!str.cell))) << " ";
    file << label(str.atom_type(which_atom(str.atom_pos,i_sg->pair(1),!str.cell))) << " ";
    rVector3d u=i_sg->pair(1)-i_sg->pair(0);
    file << norm(u) << " ";;
    u.normalize();
    file << u*(fc*u) << " ";
    rVector3d a=find_perpendicular(u);
    rVector3d b=a^u;
    file << (a*(fc*a)+b*(fc*b))/2. << endl;
    file << (!axes)*(i_sg->pair(0)) << " " << (!axes)*(i_sg->pair(1)-i_sg->pair(0)) << endl << endl;
    file << i_sg->fc_list.get_size() << endl;
    c=bc;
    for ( i_m.init(i_sg->fc_list); i_m; i_m++,c++) {
      file << beta(c) << endl << (*i_m) << endl;
    }
    file << fc << endl;
    file << "---------------" << endl << endl;
  }
}

void make_forcek(LinkedList<PairSpring> *p_pair, const Structure &str_relax, const LinkedList<SpringGeneral> &spring_data, const Array<Real> &beta, const SpaceGroup &spacegroup, int simplify=0) {
  rMatrix3d inv_cell=!str_relax.cell;
  //  NormPairEpsilon<rVector3d> identical(zero_tolerance);
  LinkedList<Array<rVector3d> > pair_done;
  LinkedListIterator<SpringGeneral> i_sg(spring_data);
  int c=0;
  SpaceGroup idmodcell;
  idmodcell.cell=str_relax.cell;
  idmodcell.point_op.resize(1);
  idmodcell.point_op(0).identity();
  idmodcell.trans.resize(1);
  idmodcell.trans(0)=rVector3d(0.,0.,0.);
  //  cerr << "mkforcek" << endl;
  for (; i_sg; i_sg++) {
    rMatrix3d fc;
    fc.zero();
    LinkedListIterator<rMatrix3d> i_m(i_sg->fc_list);
    for ( ; i_m; i_m++, c++) {
      fc=fc+(*i_m)*beta(c);
    }
    //    cerr << i_sg->pair(0)/3.57636 << " " << (i_sg->pair(1)-i_sg->pair(0))/3.57636 << endl;
    //    cerr << fc << endl;
    for (int op=0; op<spacegroup.point_op.get_size(); op++) {
      Array<rVector3d> tpair;
      apply_symmetry(&tpair,spacegroup.point_op(op),spacegroup.trans(op),i_sg->pair);
      if (add_unique(&pair_done,tpair,idmodcell)) {
        PairSpring *ps=new PairSpring;
        ps->whichatom[0]=which_atom(str_relax.atom_pos, tpair(0), inv_cell);
        ps->whichatom[1]=which_atom(str_relax.atom_pos, tpair(1), inv_cell);
        ps->dr=tpair(1)-tpair(0);
        rVector3d ideal_offset=(tpair(1)-tpair(0))-(str_relax.atom_pos(ps->whichatom[1])-str_relax.atom_pos(ps->whichatom[0]));
        ps->cell_shift=inv_cell*ideal_offset;
        ps->forcek=spacegroup.point_op(op)*(fc*(~spacegroup.point_op(op)));
	// test;
	if (simplify!=0) {
	  rVector3d u[3];
	  u[0]=ps->dr/norm(ps->dr);
	  u[1]=find_perpendicular(u[0]);
	  u[2]=u[0]^u[1];
	  rMatrix3d T;
	  for (int ii=0; ii<3; ii++) {
	    for (int jj=0; jj<3; jj++) {
	      T(ii,jj)=u[jj](ii);
	    }
	  }
	  rMatrix3d Tf=(~T)*ps->forcek*T;
	  if (simplify==1) {
	    for (int ii=0; ii<3; ii++) {
	      for (int jj=0; jj<3; jj++) {
		if (ii!=jj) {
		  Tf(ii,jj)=0;
		}
	      }
	    }
	    Tf(1,1)=(Tf(1,1)+Tf(2,2))/2.;
	    Tf(2,2)=Tf(1,1);
	  }
	  else if (simplify==2) {
	    for (int ii=0; ii<3; ii++) {
	      for (int jj=0; jj<ii; jj++) {
		Tf(ii,jj)=(Tf(ii,jj)+Tf(jj,ii))/2.;
	      }
	    }
	  }
	  ps->forcek=T*Tf*(~T);
	}
        (*p_pair) << ps;
      }
    }
  }
  //  cerr << p_pair->get_size() << endl;
}

void read_mass_file(Array<Real> *masses, istream &file, const Array<AutoString> &label) {
  masses->resize(label.get_size());
  for (int i=0; i<masses->get_size(); i++) {(*masses)(i)=-1.;}
  while (skip_delim(file)) {
    AutoString curlab;
    get_string(&curlab,file);
    Real m;
    file >> m;
    int i=index_in_array(label,curlab);
    if (i!=-1) {(*masses)(i)=m;}
  }
  for (int i=0; i<masses->get_size(); i++) {
    if ((*masses)(i)==-1.) {
      cerr << "Mass of " << label(i) << " is not known";
      ERRORQUIT("Aborting");
    }
  }
}

void read_mass_file(Array<Real> *pmasses, char *massfilename, const Structure &rel_str, const Array<AutoString> &label) {
  pmasses->resize(rel_str.atom_type.get_size());
  {
    ifstream massfile;
    Array<Real> atomic_masses;
    if (strlen(massfilename)>0) {
      massfile.open(massfilename);
      if (!massfile) {
        ERRORQUIT("Unable to open atomic masses file");
      }
    }
    else {
      AutoString configfilename(getenv("HOME"));
      configfilename+="/.atat.rc";
      ifstream configfile(configfilename);
      if (!configfile) {
        ERRORQUIT("$HOME/.atat.rc was not found.");
      }
      while (configfile.get()!='=') {};
      skip_delim(configfile," \t");
      AutoString massfilename2;
      get_string(&massfilename2,configfile);
      massfilename2+="/data/masses.in";
      massfile.open(massfilename2);
      if (!massfile) {
	ERRORQUIT("Unable to open atomic masses file");
      }
    }
    read_mass_file(&atomic_masses,massfile,label);
    for (int at=0; at<rel_str.atom_type.get_size(); at++) {
      (*pmasses)(at)=atomic_masses(rel_str.atom_type(at));
    }
  }
}

void make_pow_mat(Array2d<Real> *ppow, const Array<Real> &x, int p) {
 ppow->resize(x.get_size(),p+1);
 for (int i=0; i<x.get_size(); i++) {
   for (int j=0; j<=p; j++) {
     (*ppow)(i,j)=ipow(x(i),j);
   }
 }
}

void fit_best_poly(Array<Real> *pbeta, const Array<Real> &x, const Array<Real> &y, int *pbest_p, int maxp) {
  Real best_cv=MAXFLOAT;
  if (*pbest_p==-1) {
    for (int p=0; p<maxp; p++) {
      Array2d<Real> px;
      make_pow_mat(&px,x,p);
      Real cv=calc_cv(px,y);
      if (cv<best_cv) {
        best_cv=cv;
        *pbest_p=p;
      }
    }
  }
  Array2d<Real> px;
  make_pow_mat(&px,x,*pbest_p);
  calc_ols(pbeta,px,y);
}

void fit_eos(Array<Real> *pb, const Array<Real> &s, const Array<Real> &e, int maxp) {
  if (s.get_size()!=e.get_size()) {
    ERRORQUIT("Error in fit_eos.");
  }
  if (s.get_size()==0) {
    pb->resize(1);
    (*pb)(0)=0.;
  }
  else if (s.get_size()==1) {
    pb->resize(1);
    (*pb)(0)=e(0);
  }
  else if (s.get_size()==2) {
    pb->resize(3);
    (*pb)(0)=MIN(e(0),e(1));
    (*pb)(1)=0.;
    (*pb)(2)=fabs(e(1)-e(0))/sqr(fabs(s(1)-s(0)));
  }
  else {
    int best_p=-1;
    fit_best_poly(pb,s,e,&best_p,maxp);
    if (best_p<2) {best_p=2;}
    fit_best_poly(pb,s,e,&best_p,maxp);
  }
}

Real find_local_minimum_poly(Real *ppos, const Array<Real> &b, Real xmin, Real xmax, int n) {
  Real dx=(xmax-xmin)/(Real)(n-1);
  Real miny=MAXFLOAT;
  for (Real x=xmin; x<=xmax; x+=dx) {
    Real s=0;
    for (int i=0; i<b.get_size(); i++) {
      s+=b(i)*ipow(x,i);
    }
    if (s<miny) {
      miny=s;
      if (ppos) {(*ppos)=x;}
    }
    if (near_zero(xmax-xmin)) break;
  }
  return miny;
}

void add_poly(Array<Real> *pc, const Array<Real> &a, const Array<Real> &b) {
  pc->resize(MAX(a.get_size(), b.get_size()));
  for (int i=0; i<pc->get_size(); i++) {
    (*pc)(i)=(i<a.get_size() ? a(i) : 0.)+(i<b.get_size() ? b(i) : 0.);
  }
}

Real eval_poly(const Array<Real> &a, Real x) {
  Real y=0.;
  for (int i=0; i<a.get_size(); i++) {
    y+=a(i)*ipow(x,i);
  }
  return y;
}

int find_mode_supercell(rMatrix3d *psupercell, const rMatrix3d &cell, const rVector3d &k, Real maxvol) {
  Real cspc=1./norm(k);
  rVector3d a,b,c;
  LatticePointIterator latvect(cell);
  Real maxa=sqrt(maxvol/(sqrt(3.)/2.)/cspc)+zero_tolerance;
  while (1) {
    a=(rVector3d)latvect;
    if (norm(a)>maxa) {return 0;}
    if (near_zero(a*k)) break;
    latvect++;
  }
  Real maxb=sqrt(sqr(maxvol/norm(a)/cspc)+sqr(norm(a)/2.))+zero_tolerance;
  while (1) {
    latvect++;
    b=(rVector3d)latvect;
    if (norm(b)>maxb) {return 0;}
    if (!near_zero(norm(a^b))) {
      if (near_zero(b*k)) break;
    }
  }
  Real maxc=sqrt(sqr(cspc)+sqr(norm(a)/2.+norm(b)/2.))+zero_tolerance;
  latvect.init(cell);
  while (1) {
    c=(rVector3d)latvect;
    if (norm(c)>maxc) {return 0;}
    if (near_zero(c*k-1.)) break;
    latvect++;
  }
  psupercell->set_column(0,a);
  psupercell->set_column(1,b);
  if ((a^b)*c>0) {
    psupercell->set_column(2,c);
  }
  else {
    psupercell->set_column(2,-c);
  }
  *psupercell=find_symmetric_cell(*psupercell);
  return 1;
}

int generate_mode(Array<Structure> *pmodestr, const Structure &str, const rVector3d &k, const Array<Complex> &eigvect, Real maxvol, Real disp_mag) {
  rMatrix3d supercell;
  if (!find_mode_supercell(&supercell, str.cell,k,maxvol)) {
    pmodestr->resize(0);
    return 0;
  }
  pmodestr->resize(2);
  find_all_atom_in_supercell(&((*pmodestr)(0).atom_pos), &((*pmodestr)(0).atom_type),str.atom_pos,str.atom_type,str.cell,supercell);
  (*pmodestr)(0).cell=supercell;
  Array<Array<Real> > disp(2);
  disp(0).resize((*pmodestr)(0).atom_pos.get_size()*3);
  disp(1).resize((*pmodestr)(0).atom_pos.get_size()*3);
  rMatrix3d invcell=!(str.cell);
  for (int at=0; at<(*pmodestr)(0).atom_pos.get_size(); at++) {
    Real phase=2.*M_PI*k*((*pmodestr)(0).atom_pos(at));
    int atincell=which_atom(str.atom_pos,(*pmodestr)(0).atom_pos(at),invcell);
    for (int i=0; i<3; i++) {
      int j=3*atincell+i;
      Complex d=eigvect(j)*Complex(cos(phase),sin(phase));
      disp(0)(j)=real(d);
      disp(1)(j)=imag(d);
    }
  }
  (*pmodestr)(1)=(*pmodestr)(0);
  for (int c=0; c<2; c++) {
    Real m=max(disp(c));
    if (m>0) {
      product(&disp(c),disp(c),disp_mag/m);
    }
    for (int at=0; at<(*pmodestr)(0).atom_pos.get_size(); at++) {
      for (int i=0; i<3; i++) {
	(*pmodestr)(0).atom_pos(at)(i)+=disp(c)(3*at+i);
      }
    }
  }
  return 1;
}

extern char *helpstring;

int main(int argc, char *argv[]) {
  cerr.setf(ios::fixed | ios::showpos);
  cerr.precision(5);
  char *stridealfilename="str.out";
  char *strrelaxfilename="str_relax.out";
  Real enclosed_radius=0.;
  Real radius=0.;
  Real displ_mag=0.2;
  Real max_strain=0.01;
  int norerelax=0;
  int redun=0;
  int nb_vol=2;
  int fitfc=0;
  char *strainfilename="";
  char *massfilename="";
  Real T0=0.;
  Real T1=2000.;
  Real dT=100.;
  Real kppra=1000;
  rVector3d kshift(0.5,0.5,0.5);
  char *kdispfile="";
  Real hplanck=6.6260755e-34/1.60217733e-19;  // h in eV s (not h bar);
  Real kboltzman=1.380658e-23/1.60217733e-19; // k_B in eV/K;
  Real convfk=1.60217733e-19*1e20; // converts eV/A^2 into J/m^2;
  Real mass_unit=1.6605402e-27;    // convert a.u. mass into kg;
  Real press_ang3_to_eV=1e9*1e-30/1.60217733e-19; // convert GPa*Ang^3 to eV;
  Real robust_len=0.;
  Real perfewer=1.;
  Real pressure=0;
  int peratom=0;
  int projdos=0;
  int subE0=0;
  int forceneg=0;
  int find_unstable=0;
  int gen_unstable=0;
  int max_unstable_vol=64;
  int simplify_fc=0;
  int printenergy=0;
  int dohelp=0;
  int sigdig=5;
  int optbin=1000;

  if (file_exists("../Trange.in")) {
    ifstream tfile("../Trange.in");
    T0=0;
    int nT=0;
    tfile >> T1 >> nT;
    dT=(T1-T0)/MAX(1,nT-1);
  }

  AskStruct options[]={
    {"","Fit Stiffness VS Length transferable force constants " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-f","Fit force constants (otherwise, generate pertubations)",BOOLVAL,&fitfc},
    {"-si","Input file defining the ideal structure (default: str.out)",STRINGVAL,&stridealfilename},
    {"-sr","Input file defining the relaxed structure (default: str_relax.out)",STRINGVAL,&strrelaxfilename},
    {"-er","Minimum distance between displaced atoms",REALVAL,&enclosed_radius},
    {"-fr","Force constant range",REALVAL,&radius},
    {"-dr","Displacement of the perturbed atom (default: 0.2)",REALVAL,&displ_mag},
    {"-ms","Strain of the maximum volume sampled (default: 0.01)",REALVAL,&max_strain},
    {"-ns","Number of volume sample (default: 2)",INTVAL,&nb_vol},
    {"-nrr","Do not rerelax structures at each new volume",BOOLVAL,&norerelax},
    {"-ncs","No check for singular matrix in fit",BOOLVAL,&redun},
    {"-sf","Extra strain file",STRINGVAL,&strainfilename},
    {"-m","Input file defining the atomic masses (default: ${atatdir}/data/masses.in)",STRINGVAL,&massfilename},
    {"-T0","Minimum temperature (default: 0) ",REALVAL,&T0},
    {"-T1","Maximum temperature (default: 2000)",REALVAL,&T1},
    {"-dT","Temperature step (default: 100)",REALVAL,&dT},
    {"-P","Pressure (in GPa, default: 0)",REALVAL,&pressure},
    {"-kp","Number of k-points per reciprocal atom (default: 1000)",REALVAL,&kppra},
    {"-sx","k-point shift (along 1st recip lat. vect.)",REALVAL,&kshift(0)},
    {"-sy","k-point shift (along 2nd recip lat. vect.)",REALVAL,&kshift(1)},
    {"-sz","k-point shift (along 3rd recip lat. vect.)",REALVAL,&kshift(2)},
    {"-df","Phonon dispersion curve calculation input file.",STRINGVAL,&kdispfile},
    {"-hp","Planck's constant (default in (eV s))",REALVAL,&hplanck},
    {"-kb","Boltzman's constant (default in eV/K)",REALVAL,&kboltzman},
    {"-cfk","Conversion factor for force constants into energy/dist^2 (default: converts eV/A^2 into J/m^2)",REALVAL,&convfk},
    {"-mu","Mass units (default: converts a.u. mass into kg)",REALVAL,&mass_unit},
    {"-rl","Robust Length algorithm parameter for soft modes (beta)",REALVAL,&robust_len},
    {"-cP","Conversion factor from pressure*volume into eV",REALVAL,&press_ang3_to_eV},
    {"-pa","Output free energy per atom instead of per unit cell",BOOLVAL,&peratom},
    {"-pe","Print energy in 3rd column of fitfc.out",BOOLVAL,&printenergy},
    {"-sc","Correction factor if spectator ion are present (default: 1)",REALVAL,&perfewer},
    {"-me0","Subtract energy at 0K",BOOLVAL,&subE0},
    {"-fn","Force continuation of calculations even if unstable",BOOLVAL,&forceneg},
    {"-fu","Find unstable modes",BOOLVAL,&find_unstable},
    {"-gu","Genenerate unstable modes number n",INTVAL,&gen_unstable},
    {"-mau","Maximum number of atom per supercell for unstable mode generation",INTVAL,&max_unstable_vol},
    {"-sfc","Simplify force constants: 1=streching+bending, 2=symmetric",INTVAL,&simplify_fc},
    {"-apd","Atom-Projected DOS",BOOLVAL,&projdos},
    {"-sig","Number of significant digits printed (default: 5)",INTVAL,&sigdig},
    {"-z","Tolerance for finding symmetry operations (default: 1e-3)",REALVAL,&zero_tolerance},
    {"-h","Display more help",BOOLVAL,&dohelp}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    return 1;
  }
  if (dohelp) {
    cout << helpstring;
    return 1;
  }
  if (fitfc && file_exists("../Trange.in")) {
    cerr << "Read Trange.in file." << endl;
  }

  if (enclosed_radius==0 && !fitfc) {
    ERRORQUIT("You need to specify at least the -er parameter to generate perturbations or the -f option to fit the force constants.");
  }
  if (radius==0 && fitfc) {
    ERRORQUIT("You need to specify at least the -fr parameter with the -f option.");
  }
  int nT=(int)((T1-T0)/dT)+1;

  rMatrix3d extrastrain;
  extrastrain.zero();
  if (strlen(strainfilename)>0) {
    ifstream strainfile(strainfilename);
    strainfile >> extrastrain;
  }
  Structure str;
  Array<AutoString> atom_label;
  rMatrix3d axes;
  {
    Array<Arrayint> site_type_list;
    ifstream file(stridealfilename);
    if (!file) ERRORQUIT("Unable to open ideal structure file.");
    parse_lattice_file(&str.cell, &str.atom_pos, &str.atom_type, &site_type_list, &atom_label, file, &axes);
    if (fabs(det(str.cell))<zero_tolerance) ERRORQUIT("Lattice vectors are coplanar.");
    wrap_inside_cell(&str.atom_pos,str.atom_pos,str.cell);
    fix_atom_type(&str,site_type_list);
    remove_vacancies(&str,atom_label);
    strain_str(&str,str,extrastrain);
  }

  Structure str_relax;
  rMatrix3d axes_relax;
  {
    ifstream file(strrelaxfilename);
    if (!file) ERRORQUIT("Unable to open relaxed structure file.");
    parse_structure_file(&(str_relax.cell), &(str_relax.atom_pos), &(str_relax.atom_type), atom_label, file, &axes_relax);
    if (fabs(det(str_relax.cell))<zero_tolerance) ERRORQUIT("Lattice vectors are coplanar.");
    wrap_inside_cell(&str_relax.atom_pos,str_relax.atom_pos,str_relax.cell);
    reorder_atoms(&str_relax,str);
  }

  if (!fitfc) {
    rMatrix3d ideal_supercell;
    find_smallest_supercell_enclosing_sphere(&ideal_supercell,str.cell,enclosed_radius);
    cerr << "Writing volumes..." << endl;
    for (int iv=0; iv<nb_vol; iv++) {
      Real strain=max_strain*(Real)iv/(Real)MAX(nb_vol-1,1);
      ostrstream volname;
      if (iv==0) {
	volname << "vol_0" << '\0';
      }
      else {
	volname << "vol_" << strain*100 << '\0';
      }
      mkdir(volname.str(),S_IRWXU | S_IRWXG | S_IRWXO);
      chdir_robust(volname.str());
      Structure s_str_relax;
      stretch_str(&s_str_relax,str_relax,strain);
      if (!file_exists("str_relax.out")) {
	cerr << volname.str() << endl;
	{
          ofstream file("str.out");
          write_structure(s_str_relax, atom_label, axes_relax*(1.+strain), file);
	}
        if (norerelax) {
          ofstream file("str_relax.out");
          write_structure(s_str_relax, atom_label, axes_relax*(1.+strain), file);
        }
        else {
          ofstream waitfile("wait");
        }
      }
      chdir_robust("..");
    }
    system("ls vol_*/str_relax.out 2> /dev/null | sed 's+/str_relax.out++g' > vollist.out");
    {
      cerr << "Writing perturbations..." << endl;
      ifstream volfile("vollist.out");
      int atleastone=0;
      while (skip_delim(volfile)) {
	atleastone=1;
        AutoString dirname;
        get_string(&dirname,volfile);
        chdir_robust(dirname);
	cerr << dirname << endl;
        Structure s_str_relax;
        ifstream file("str_relax.out");
	rMatrix3d rel_axes;
        parse_structure_file(&(s_str_relax.cell), &(s_str_relax.atom_pos), &(s_str_relax.atom_type), atom_label, file, &rel_axes);
        wrap_inside_cell(&s_str_relax.atom_pos,s_str_relax.atom_pos,s_str_relax.cell);
        reorder_atoms(&s_str_relax,str_relax);
        SpaceGroup spacegroup;
        spacegroup.cell=s_str_relax.cell;
        find_spacegroup(&spacegroup.point_op,&spacegroup.trans,s_str_relax.cell,s_str_relax.atom_pos,s_str_relax.atom_type);
        Array<int> pert_atom;
        Array<rVector3d> pert_direction;
        Array<int> do_minus;
        find_pertubations(&pert_atom, &pert_direction, &do_minus,s_str_relax,spacegroup);
        for (int p=0; p<pert_atom.get_size(); p++) {
          for (int m=1; m>=(1-2*do_minus(p)); m-=2) {
	    ostrstream pertname;
            pertname << "p" << (m==-1 ? '-' : '+' ) << displ_mag << "_" << enclosed_radius << "_" << p << '\0';
            mkdir(pertname.str(),S_IRWXU | S_IRWXG | S_IRWXO);
            chdir_robust(pertname.str());
	    if (!file_exists("str.out")) {
              cerr << " " << pertname.str() << endl;
	      Structure unpert_str;
	      rMatrix3d supercell=(s_str_relax.cell)*((!str.cell)*ideal_supercell);
	      rVector3d zerovect(0.,0.,0.);
	      calc_pert_str(&unpert_str,supercell,s_str_relax,zerovect,zerovect,0);
	      {
		ofstream unpertstrfile("str_unpert.out");
		unpertstrfile.setf(ios::fixed);
		unpertstrfile.precision(sigdig);
		write_structure(unpert_str, atom_label, rel_axes, unpertstrfile);
	      }
	      
	      Structure ideal_str;
	      calc_pert_str(&ideal_str,ideal_supercell,str,zerovect,zerovect,0.);
	      {
		ofstream idealstrfile("str_ideal.out");
		idealstrfile.setf(ios::fixed);
		idealstrfile.precision(sigdig);
		write_structure(ideal_str, atom_label, axes, idealstrfile);
	      }
	      
	      Structure pert_str;
	      calc_pert_str(&pert_str,supercell,s_str_relax,s_str_relax.atom_pos(pert_atom(p)),(Real)m*displ_mag*pert_direction(p),0);
	      {
		ofstream pertstrfile("str.out");
		pertstrfile.setf(ios::fixed);
		pertstrfile.precision(sigdig);
		write_structure(pert_str, atom_label, rel_axes, pertstrfile);
	      }
	      ofstream waitfile("wait");
	    }
            chdir_robust("..");
          }
        }        
        chdir_robust("..");
      }
      if (!atleastone) {
	cerr << "Next, you need to relax the structures output in the vol_* subdirectories" << endl
	     << "  and rerun fitfc with the same command-line options." << endl;
	cerr << "Exception: If the initial structure was already relaxed and either" << endl
	     << "  (i) you are performing harmonic calculations or" << endl
	     << "  (ii) there are no cell-internal degrees of freedom, rerun fitfc with the -nrr option added." << endl;
      }
    }
    unlink("vollist.out");
  }
  else { // fit forcek instead of generate pert;
    Array<Real> masses;
    read_mass_file(&masses,massfilename,str_relax,atom_label);

    LinkedList<Real> s_elist;
    LinkedList<Real> s_flist;
    LinkedList<Real> energy_list;
    LinkedList<Array<Real> > fvib_list;
    cerr << "Reading..." << endl;
    system("ls vol_*/str_relax.out 2> /dev/null | sed 's+/str_relax.out++g' > vollist.out");
    {
      ifstream volfile("vollist.out");
      while (skip_delim(volfile)) {
        AutoString volname;
        get_string(&volname,volfile);
	cerr << volname << endl;
        chdir_robust(volname);
        ofstream logfile("fitfc.log");
        logfile.setf(ios::fixed | ios::showpos);
        logfile.precision(sigdig);
        Structure s_str_relax;
        ifstream file("str_relax.out");
	rMatrix3d rel_axes;
        parse_structure_file(&(s_str_relax.cell), &(s_str_relax.atom_pos), &(s_str_relax.atom_type), atom_label, file, &rel_axes);
        wrap_inside_cell(&s_str_relax.atom_pos,s_str_relax.atom_pos,s_str_relax.cell);
        reorder_atoms(&s_str_relax,str_relax);
	
	Real strain=pow(det(s_str_relax.cell)/det(str_relax.cell),1./3.)-1.;
	if (file_exists("energy")) {
	  s_elist << new Real(strain);
          ifstream energyfile("energy");
          Real e=MAXFLOAT;
          energyfile >> e;
	  if (e==MAXFLOAT) {ERRORQUIT("Unable to read energy file.");}
          e=e/(Real)(str.atom_pos.get_size())*(Real)(peratom ? 1:str_relax.atom_type.get_size())*perfewer;
          energy_list << new Real(e);
        }
        SpaceGroup spacegroup;
        spacegroup.cell=s_str_relax.cell;
        find_spacegroup(&spacegroup.point_op,&spacegroup.trans,s_str_relax.cell,s_str_relax.atom_pos,s_str_relax.atom_type);
        LinkedList<SpringGeneral> fc_list;
        find_distinct_forcek(&fc_list,radius,str,s_str_relax,spacegroup);

        LinkedList<Structure> pert_str_list;
        LinkedList<Structure> unpert_str_list;
        LinkedList<Array<Real> > forcevector_list;
        int nb_row=0;
        system("ls p*/force.out 2> /dev/null | sed 's+/force.out++g' > pertlist.out");
        ifstream pertfile("pertlist.out");
        while (skip_delim(pertfile)) {
          AutoString dirname;
          get_string(&dirname,pertfile);
	  cerr << " " << dirname << endl;
          chdir_robust(dirname);
          Structure *punpert_str=new Structure;
          {
            ifstream file("str_unpert.out");
            if (!file) {
              cerr << "Unable to open str_unpert.out in directory " << volname << "/" << dirname << endl;
              ERRORQUIT("Aborting");
            }
            parse_structure_file(&(punpert_str->cell),&(punpert_str->atom_pos),&(punpert_str->atom_type),atom_label,file);
            if (punpert_str->atom_pos.get_size()==0) {ERRORQUIT("Unable to read str_unpert.out");}
            unpert_str_list << punpert_str;
          }
          Array<int> copyfrom;
          Array<iVector3d> cellshift;
          Structure *ppert_str=new Structure;
          {
            ifstream file("str_relax.out");
            if (!file) {
              cerr << "Unable to open str_relax.out in directory " << volname << "/" << dirname << endl;
              ERRORQUIT("Aborting");
            }
            parse_structure_file(&(ppert_str->cell),&(ppert_str->atom_pos),&(ppert_str->atom_type),atom_label,file);
            if (ppert_str->atom_pos.get_size()==0) {ERRORQUIT("Unable to read str_relax.out");}
            reorder_atoms(ppert_str,*punpert_str,&copyfrom,&cellshift);
	    // write_structure(*ppert_str,atom_label, axes, cerr, 0);
            pert_str_list << ppert_str;
	  }
	  {
	    ifstream forcefile("force.out");
	    if (!forcefile) {
	      cerr << "Unable to open force.out in directory " << volname << "/" << dirname << endl;
	      ERRORQUIT("Aborting");
	    }
	    Array<Real> *pforce=new Array<Real>;
	    read_force_vector(pforce,forcefile,copyfrom);
	    forcevector_list << pforce;
	    nb_row+=pforce->get_size();
	    Real dE=calc_diff_energy(*punpert_str,*ppert_str,*pforce);
	    if (dE<0) {
	      cerr << "Warning: " << dirname << " is an unstable mode." << endl;
	    }
	    logfile << dirname << " dE= " << dE << endl;
	  }
	  chdir_robust("..");
        }
        unlink("pertlist.out");
	if (nb_row>0) {
	  s_flist << new Real(strain);
	  int nb_col=0;
	  LinkedListIterator<SpringGeneral> i_s(fc_list);
	  for (;i_s; i_s++) {
	    nb_col+=i_s->fc_list.get_size();
	  }

	  Array2d<Real> big_eqnmat(nb_row,nb_col);
	  Array<Real> big_force(nb_row);
	  LinkedListIterator<Structure> i_pert_str(pert_str_list);
	  LinkedListIterator<Structure> i_unpert_str(unpert_str_list);
	  LinkedListIterator<Array<Real> > i_forcevector(forcevector_list);
	  for (int r=0; i_pert_str; i_pert_str++, i_unpert_str++, i_forcevector++) {
	    Array<Array<Real> > eqn;
	    calc_equation_matrix(&eqn,*i_unpert_str,*i_pert_str,fc_list,spacegroup);
	    for (int i=0; i<i_forcevector->get_size(); i++) {
	      big_force(r)=(*i_forcevector)(i);
	      for (int j=0; j<nb_col; j++) {
		big_eqnmat(r,j)=eqn(i)(j);
	      }
	      r++;
	    }
	  }
	  
	  logfile << "Equation Matrix" << endl;
	  logfile << big_eqnmat << endl;
	  logfile << endl;
	  
	  Array<Real> beta;
	  Array<Real> w;
	  calc_ols(&beta,big_eqnmat,big_force,w,!redun);
	  
	  Array<Real> force_hat;
	  product(&force_hat,big_eqnmat,beta);
	  
	  logfile << "actual force   predicted force   difference" << endl;
	  for (int i=0; i<big_force.get_size(); i++) {
	    logfile << big_force(i) << " " << force_hat(i) << " " << (big_force(i)-force_hat(i)) << endl;
	  }
	  logfile << "end" << endl;
	  
	  {
	    ofstream file("fc.out");
	    file.setf(ios::fixed | ios::showpos);
	    file.precision(sigdig);
	    print_fc(file,fc_list,beta,s_str_relax,atom_label,rel_axes);
	  }
	  
	  LinkedList<PairSpring> pairspring;
	  make_forcek(&pairspring,s_str_relax,fc_list,beta,spacegroup,simplify_fc);
	  
	  if (strlen(kdispfile)>0) {
	    chdir_robust("..");
	    ifstream infile(kdispfile);
	    if (!infile) {
	      ERRORQUIT("Unable to open dispersion curve calculation input file.");
	    }
	    chdir_robust(volname);
	    // rMatrix3d recip=(!(~((s_str_relax.cell)*(!str.cell)*rel_axes)));
	    rMatrix3d recip=(!(~((s_str_relax.cell)*(!str.cell)*axes)));
	    // rMatrix3d recip=!(~(rel_axes));
	    ofstream outffile("eigenfreq.out");
	    ofstream outvfile("eigenvect.out");
	    while (skip_delim(infile)) {
	      int nbpts;
	      rVector3d k1,k2;
	      infile >> nbpts >> k1 >> k2;
	      k1=recip*k1;
	      k2=recip*k2;
	      LinkedList<Array<Real> > freq;
	      LinkedList<Array2d<Complex> > eigvect;
	      calc_dispersion_curve(&freq,&eigvect, k1,k2,nbpts,pairspring,masses,convfk/mass_unit);
	      LinkedListIterator<Array<Real> > i_freq(freq);
	      LinkedListIterator<Array2d<Complex> > i_eigvect(eigvect);
	      for (; i_freq; i_freq++, i_eigvect++) {
		for (int i=0; i<i_freq->get_size(); i++) {
		  outffile << (*i_freq)(i) << '\t';
		}
		outffile << endl;
		outvfile << (*i_eigvect);
	      }
	    }
	  }

	  if (find_unstable) {
	    iVector3d kmesh;
	    calc_k_mesh(&kmesh,str.cell,kppra/(Real)(str.atom_pos.get_size()));
	    cerr << "kmesh= " << kmesh << endl;
	    rMatrix3d rec_lat=!(~(s_str_relax.cell));
	    rMatrix3d kcell;
	    for (int i=0; i<3; i++) {
	      kcell.set_column(i,rec_lat.get_column(i)/(Real)kmesh(i));
	    }
	    ofstream modefile("unstable.out");
	    int line=0;
	    MultiDimIterator<iVector3d> ik(kmesh);
	    for (; ik; ik++) {
	      rVector3d k=kcell*(kshift+to_real(ik));
	      k=flip_into_brillouin_1(k,rec_lat);
	      Array<Real> freq;
	      Array2d<Complex> eigvects;
	      calc_normal_modes(&freq,&eigvects, k,pairspring,masses,convfk/mass_unit);
	      for (int i=0; i<freq.get_size(); i++) {
		line++;
		if (freq(i)<0) {
		  Array<Complex> eigvect;
		  extract_column(&eigvect,eigvects,i);
		  Real maxvol=fabs(det(s_str_relax.cell))*(Real)max_unstable_vol/(Real)(s_str_relax.atom_pos.get_size());
		  rMatrix3d supercell;
		  if (find_mode_supercell(&supercell, s_str_relax.cell,k,maxvol)) {
		    int nb_atom=iround((Real)(s_str_relax.atom_pos.get_size())*fabs(det(supercell)/det(s_str_relax.cell)));
		    modefile << "u " << line << " " << nb_atom << " " << (iVector3d &)ik << " " << i << " " << freq(i) << " " << eigvect << endl;
		  }
		  else {
		    modefile << "u " << line << " " << "too_large" << " " << (iVector3d &)ik << " " << i << " " << freq(i) << " " << eigvect << endl;
		  }
		}
	      }
	    }
	  }
	  else if (gen_unstable) {
	    iVector3d kmesh;
	    calc_k_mesh(&kmesh,str.cell,kppra/(Real)(str.atom_pos.get_size()));
	    cerr << "kmesh= " << kmesh << endl;
	    rMatrix3d rec_lat=!(~(s_str_relax.cell));
	    rMatrix3d kcell;
	    for (int i=0; i<3; i++) {
	      kcell.set_column(i,rec_lat.get_column(i)/(Real)kmesh(i));
	    }
	    int line=0;
	    int nbbr=str.atom_pos.get_size()*3;
	    MultiDimIterator<iVector3d> ik(kmesh);
	    for (; ik; ik++) {
	      if (line<abs(gen_unstable) && abs(gen_unstable)<=line+nbbr) {
		rVector3d k=kcell*(kshift+to_real(ik));
		k=flip_into_brillouin_1(k,rec_lat);
		Array<Real> freq;
		Array2d<Complex> eigvects;
		calc_normal_modes(&freq,&eigvects, k,pairspring,masses,convfk/mass_unit);
		for (int i=0; i<freq.get_size(); i++) {
		  line++;
		  if (line==abs(gen_unstable)) {
		    Array<Complex> eigvect;
		    extract_column(&eigvect,eigvects,i);
		    Real maxvol=fabs(det(s_str_relax.cell))*(Real)max_unstable_vol/(Real)(s_str_relax.atom_pos.get_size());
		    ostrstream dirname;
		    dirname << "p_uns_" << displ_mag << "_" << kppra << "_" << gen_unstable << '\0';
		    mkdir(dirname.str(),S_IRWXU | S_IRWXG | S_IRWXO);
		    chdir_robust(dirname.str());
		    
		    Array<Structure> pertstr(2);
		    generate_mode(&pertstr, s_str_relax,k,eigvect,maxvol,displ_mag);
		    {
		      ofstream pertstrfile("str.out");
		      pertstrfile.setf(ios::fixed);
		      pertstrfile.precision(sigdig);
		      write_structure(pertstr((gen_unstable>0 ? 0 : 1)), atom_label, rel_axes, pertstrfile);
		    }
		    Array<Structure> unpertstr(2);
		    generate_mode(&unpertstr, s_str_relax,k,eigvect,maxvol,0.);
		    {
		      ofstream unpertstrfile("str_unpert.out");
		      unpertstrfile.setf(ios::fixed);
		      unpertstrfile.precision(sigdig);
		      write_structure(unpertstr(0), atom_label, rel_axes, unpertstrfile);
		    }
		    
		    ofstream waitfile("wait");
		    chdir_robust("..");
		  }
		}
	      }
	      else {
		line+=nbbr;
	      }
	    }
	  }
	  else if (projdos) {
	    iVector3d kmesh;
	    calc_k_mesh(&kmesh,str.cell,kppra/(Real)(str.atom_pos.get_size()));
	    rMatrix3d kcell=!(~s_str_relax.cell);
	    for (int i=0; i<3; i++) {
	      kcell.set_column(i,kcell.get_column(i)/(Real)kmesh(i));
	    }
	    LinkedList<Real> freqs;
	    LinkedList<Array<Real> > weights;
	    MultiDimIterator<iVector3d> ik(kmesh);
	    logfile << "Phonon frequencies/modes/weights:" << endl;
	    for (; ik; ik++) {
	      rVector3d k=kcell*(kshift+to_real(ik));
	      Array<Real> freq;
	      Array2d<Complex> eigvects;
	      calc_normal_modes(&freq,&eigvects, k,pairspring,masses,convfk/mass_unit);
	      logfile << freq << endl;
	      logfile << eigvects << endl;
	      for (int m=0; m<freq.get_size(); m++) {
		freqs << new Real(freq(m));
		Array<Complex> eigvect;
		extract_column(&eigvect,eigvects,m);
		Array<Real> weight(str.atom_type.get_size());
		Real s=0;
		for (int at=0; at<str.atom_type.get_size(); at++) {
		  weight(at)=sqr(abs(eigvect(3*at+0)))+sqr(abs(eigvect(3*at+1)))+sqr(abs(eigvect(3*at+2)));
		  s+=weight(at);
		}
		product(&weight,weight,3./(s*(Real)(str.atom_pos.get_size())));
		logfile << weight << endl;
		weights << new Array<Real>(weight);
	      }
	    }
	    logfile << "end" << endl << endl;
	    Real xmin,xmax;
	    Array<Array<Real> > pvdos(str.atom_type.get_size());
	    for (int at=0; at<str.atom_type.get_size(); at++) {
	      Array<Real> freqa;
	      LinkedList_to_Array(&freqa,freqs);
	      Array<Real> weighta(freqa.get_size());
	      LinkedListIterator<Array<Real> > it(weights);
	      for (int i=0; it; i++,it++) {
		weighta(i)=(*it)(at);
	      }
	      smooth_density(&xmin,&xmax,&(pvdos(at)),freqa,weighta);
	    }
	    Real df=(xmax-xmin)/(Real)(pvdos(0).get_size()-1);
	    Real w=(1./(kmesh(0)*kmesh(1)*kmesh(2)))*(Real)(peratom ? 1:str_relax.atom_type.get_size())*perfewer;
	    {
	      ofstream vdosfile("pvdos.out");
	      for (int i=0; i<pvdos(0).get_size(); i++) {
		vdosfile << xmin+df*(Real)i;
		for (int at=0; at<str.atom_type.get_size(); at++) {
		  vdosfile << "\t" << w*pvdos(at)(i);
		}
		vdosfile << endl;
	      }
	    }
	  }
	  else {
	    iVector3d kmesh;
	    calc_k_mesh(&kmesh,str.cell,kppra/(Real)(str.atom_pos.get_size()));
	    LinkedList<Real> freq,weight;
	    if (!list_phonon_freq(&freq, &weight, kmesh, kshift, pairspring, s_str_relax.cell,masses,convfk/mass_unit)) {
	      cerr << "Unstable modes found." << endl;
	      if (!forceneg && robust_len==0.) {
		ERRORQUIT("Aborting.");
	      }
	    }
	  
	    LinkedListIterator<Real> ifreq(freq);
	    logfile << "Phonon frequencies:" << endl;
	    for ( ; ifreq; ifreq++) {
	      logfile << *ifreq << endl;
	    }
	    logfile << "end" << endl << endl;
	    {
	      ofstream vdosfile("vdos.out");
	      Array<Real> freqa;
	      LinkedList_to_Array(&freqa,freq);
	      Real xmin,xmax;
	      Array<Real> vdos;
	      smooth_density(&xmin,&xmax,&vdos,freqa);
	      Real df=(xmax-xmin)/(Real)(vdos.get_size()-1);
	      LinkedListIterator<Real> iw(weight);
	      Real w=(*iw)*(Real)(peratom ? 1:str_relax.atom_type.get_size())*perfewer;
	      for (int i=0; i<vdos.get_size(); i++) {
		vdosfile << xmin+df*(Real)i << "\t" << w*vdos(i) << endl;
	      }
	      vdosfile << endl;
	    }
	    {
	      ofstream file("svib_ht");
	      file << calc_vib_entropy(freq,weight)*(Real)(peratom ? 1:str_relax.atom_type.get_size())*perfewer << endl;
	    }
	    Array<Real> fvib(nT);
	    {
	      Real energy=MAXFLOAT;
	      if (file_exists("energy")) {
		ifstream energyfile("energy");
		energyfile >> energy;
	      }
	      if (energy==MAXFLOAT) {ERRORQUIT("Unable to read energy file.");}
	      energy=energy/(Real)(str_relax.atom_pos.get_size()) *(Real)(peratom ? 1:str_relax.atom_type.get_size())*perfewer;
	      ofstream file("fvib");
	      for (int it=0; it<nT; it++) {
		Real T=T0+(Real)it*dT;
		if (robust_len>0) {
		  Real geomass=1.;
		  for (int at=0; at<masses.get_size(); at++) {geomass*=masses(at);}
		  geomass=mass_unit*pow(geomass,1./(Real)(masses.get_size()));
		  Real L=pow(det(s_str_relax.cell)/(Real)(s_str_relax.atom_type.get_size())*exp(1.),1./3.)*sqrt(geomass/convfk)/2.;
		  fvib(it)=calc_vib_free_energy_robust(freq,weight,kboltzman*T,hplanck,robust_len*L);
		}
		else {
		  fvib(it)=calc_vib_free_energy(freq,weight,kboltzman*T,hplanck);
		}
                fvib(it)+=pressure*det(s_str_relax.cell)/(Real)(str_relax.atom_type.get_size())*press_ang3_to_eV;
		fvib(it)*=(Real)(peratom ? 1:str_relax.atom_type.get_size())*perfewer;
		file << (subE0 ? 0 : energy)+fvib(it) << endl;;
	      }
	    }
	    fvib_list << new Array<Real>(fvib);
	  }
	}
	chdir_robust("..");
      }
    }
    unlink("vollist.out");
    if (strlen(kdispfile)==0 && gen_unstable==0 && !projdos) {
      ofstream outfile("fitfc.out");
      ofstream fvibfile("fvib");
      Array<Real> energy_a;
      Array<Array<Real> > fvib_a;
      Array<Real> s_e_a;
      Array<Real> s_f_a;
      LinkedList_to_Array(&energy_a,energy_list);
      LinkedList_to_Array(&fvib_a,fvib_list);
      LinkedList_to_Array(&s_e_a,s_elist);
      LinkedList_to_Array(&s_f_a,s_flist);
      if (fvib_a.get_size()==0) {
	ERRORQUIT("No vol_*/p*/force.out files found or no vol_*/str_relax.out files found.");
      }
      if (s_e_a.get_size()==0) {
	cerr << "Warning: no vol_*/energy file found. Consequently, T=0 static energy will not be added to output free energies." << endl;
	s_e_a.resize(1);
	energy_a.resize(1);
	s_e_a(0)=0.;
	energy_a(0)=0.;
      }
      Real smin=MIN(min(s_e_a),min(s_f_a));
      Real smax=MAX(max(s_e_a),max(s_f_a));
      Array<Real> eos;
      if (!file_exists("user_eos.in")) {
	fit_eos(&eos,s_e_a,energy_a,4);
      }
      else {
	ifstream eosfile("user_eos.in");
	LinkedList<Real> epoly;
	smin=0;
	eosfile >> smax;
	while (skip_delim(eosfile)) {
	  Real c;
	  eosfile >> c;
	  // c=c/(Real)(str.atom_pos.get_size())*(Real)(peratom ? 1:str_relax.atom_type.get_size())*perfewer;
	  epoly << new Real(c);
	}
	LinkedList_to_Array(&eos,epoly);
      }

      {
	ofstream eosfile("eos0.out");
	for (int j=0; j<energy_a.get_size(); j++ ) {
	  eosfile << s_e_a(j) << " " << energy_a(j) << endl;
	}
	ofstream eosgnufile("eos0.gnu");
	eosgnufile << "plot 'eos0.out',";
	for (int j=0; j<eos.get_size(); j++) {
	  if (j!=0) {eosgnufile << "+";}
	  eosgnufile << eos(j) << "*x**(" << j << ")";
	}
	eosgnufile << endl;
	eosgnufile << "pause -1" << endl;
      }

      Real E0=find_local_minimum_poly(NULL,eos,smin,smax,optbin);
      Array<Real> f_vs_T(nT);
      ofstream tfvibfile("tfvib.gnu");
      ofstream tfvibdatafile("tfvib.out");
      ofstream teosfile("teos.gnu");
      for (int it=0; it<nT; it++) {
	Real T=(T0+(Real)it*dT);
        Array<Real> f_vs_s(s_f_a.get_size());
        for (int iv=0; iv<s_f_a.get_size(); iv++) {
          f_vs_s(iv)=fvib_a(iv)(it);
	  tfvibdatafile << T << " " << s_f_a(iv) << " " << f_vs_s(iv) << endl;
        }
/*
	Array2d<Real> ps;
	make_pow_mat(&ps,s_f_a,s_f_a.get_size()-1);
	Array<Real> eos_vib;
	calc_ols(&eos_vib,ps,f_vs_s);
*/
        Array<Real> eos_vib;
        int bestp=-1;
	if (s_f_a.get_size()<=2) {bestp=s_f_a.get_size()-1;}
	fit_best_poly(&eos_vib, s_f_a, f_vs_s, &bestp, s_f_a.get_size()-1);

	tfvibfile << "plot [" << smin << ":" << smax << "] 'tfvib.out' u ($2):(($1==" << T << " ? $3 : 1./0.)), ";
	for (int j=0; j<eos_vib.get_size(); j++) {
	  if (j!=0) {tfvibfile << "+";}
	  tfvibfile << eos_vib(j) << "*x**(" << j << ")";
	}
	tfvibfile << endl;
	tfvibfile << "pause -1" << endl;

	Array<Real> poly;
	add_poly(&poly,eos,eos_vib);
	teosfile << "plot [" << smin << ":" << smax << "]";
	for (int j=0; j<poly.get_size(); j++) {
	  if (j!=0) {teosfile << "+";}
	  teosfile << poly(j) << "*x**(" << j << ")";
	}
	teosfile << endl;
	teosfile << "pause -1" << endl;
        Real pos=0;
        f_vs_T(it)=find_local_minimum_poly(&pos,poly,smin,smax,optbin);
	Real E=eval_poly(eos,pos);
	if (subE0) {f_vs_T(it)-=E0; E-=E0;}
        outfile << (T0+(Real)it*dT) << " " << f_vs_T(it);
	if (printenergy) {outfile << " " << E;}
	outfile << " " << pos << endl;
        fvibfile << f_vs_T(it) << endl;
      }
      {
	ofstream sfile("svib");
	sfile << 0. << endl;
        for (int it=1; it<nT-1; it++) {
	  sfile << -(f_vs_T(it+1)-f_vs_T(it-1))/(2*dT) << endl;
	}
	sfile << -(f_vs_T(nT-1)-f_vs_T(nT-2))/dT << endl;
      }
    }
  }
}
