#include <fstream.h>
#include <sys/stat.h>
#include "phonlib.h"
#include "getvalue.h"
#include "parse.h"
#include "version.h"
#include "lattype.h"
#include "arraylist.h"

void calc_force_vector(Array<Real> *forces, const Structure &unpert_str, const Structure &pert_str, const LinkedList<PairSpring> &springlist) {
  forces->resize(unpert_str.atom_pos.get_size()*3);
  zero_array(forces);
  LinkedListIterator<PairSpring> ip(springlist);
  for (; ip; ip++) {
    rMatrix3d fk=ip->forcek;
    for (int wend=0; wend<2; wend++) {
      rVector3d dr=(pert_str.atom_pos(ip->whichatom[wend])-unpert_str.atom_pos(ip->whichatom[wend]));
//if (norm(dr)>zero_tolerance) {
//  cerr << dr << endl;
//  cerr << fk << endl;
//}
      rVector3d f=fk*dr;
      add_rVector3d_to_ArrayReal(forces,ip->whichatom[1-wend],f);
      add_rVector3d_to_ArrayReal(forces,ip->whichatom[  wend],-f);
      fk=~fk;
    }
  }
}

class SameEndAtom {
 public:
  int operator () (const SpringSvsL &a, const SpringSvsL &b) const {
    if ((a.atom_type[0]==b.atom_type[0] && a.atom_type[1]==b.atom_type[1])
        || (a.atom_type[0]==b.atom_type[1] && a.atom_type[1]==b.atom_type[0])) return 1;
    return 0; 
  }
};

void find_pertubations(Array<int> *pwhich_atom, Array<rVector3d> *pdirection,
		       const Structure &str,
		       const SpaceGroup &spacegroup) {
  LinkedList<int> which_atom_list;
  LinkedList<rVector3d> dir_list;
  LinkedList<Array<rVector3d> > pointlist;
  for (int at=0; at<str.atom_pos.get_size(); at++) {
    Array<rVector3d> p(1);
    p(0)=str.atom_pos(at);
    if (add_unique(&pointlist,p,spacegroup)) {
      Array<rMatrix3d> atom_point_group;
      find_atom_point_group(&atom_point_group,str.atom_pos(at),str.cell,spacegroup.point_op,spacegroup.trans);
      Array<rVector3d> special_dir;
      find_special_direction(&special_dir, NULL, atom_point_group);
      for (int i=0; i<special_dir.get_size(); i++) {
        which_atom_list << new int(at);
        dir_list << new rVector3d(special_dir(i));
      }
    }
  }
  LinkedList_to_Array(pwhich_atom,which_atom_list);
  LinkedList_to_Array(pdirection,dir_list);
}

void calc_equation_matrix(Array<Array<Real> > *peqn_mat, const Structure &ideal_str, const Structure &unpert_str, const Structure &pert_str, LinkedList<SpringSvsL> *pspring_list, const LinkedList<rMatrix3d> &dirdep_mat, const Array<Real> &conc, const MultiDimPoly &multipoly, Real maxfklen=MAXFLOAT) {
  int nbcol=0;
  LinkedListIterator<SpringSvsL> i_s(*pspring_list);
  for (;i_s; i_s++) {
    for (int sb=0; sb<i_s->forcek.get_size(); sb++) {
      nbcol+=i_s->forcek(sb).get_size();
      zero_array(&(i_s->forcek(sb)));
   }
  }
  int nbrow=pert_str.atom_pos.get_size()*3;
  resize(peqn_mat,nbrow,nbcol);
  int c=0;
  LinkedList<Array<rVector3d> > pairbuf;

  i_s.init(*pspring_list);
  for (;i_s; i_s++) {
    for (int sb=0; sb<i_s->forcek.get_size(); sb++) {
      for (int np=0; np<i_s->forcek(sb).get_size(); np++) {
	i_s->forcek(sb)(np)=1.;
	LinkedList<PairSpring> pair_list;
//	make_nn_forcek(&pair_list, ideal_str, unpert_str, *pspring_list, dirdep_mat, conc, multipoly);
	make_nn_forcek(&pair_list, ideal_str, unpert_str, *pspring_list, dirdep_mat, conc, multipoly, &pairbuf);
	{
	  LinkedListIterator<PairSpring> ifk(pair_list);
	  while (ifk) {
	    if (norm(ifk->dr)>maxfklen) {
	      delete pair_list.detach(ifk);
	    }
	    else {
	      ifk++;
	    }
	  }
	}
	Array<Real> force;
	calc_force_vector(&force, unpert_str, pert_str, pair_list);
        for (int i=0; i<nbrow; i++) {
          (*peqn_mat)(i)(c)=force(i);
        }
        c++;
	i_s->forcek(sb)(np)=0.;
      }
    }
  }
}

void write_spring_file(ostream &file, const LinkedList<SpringSvsL> &springlist, const Array<AutoString> &label, int maxlpow, int dirdep, int concpow) {
  if (dirdep!=0 || concpow!=0) {
    file << "maxlpow " << maxlpow << endl;
    file << "dirdep " << dirdep << endl;
    file << "concpow " << concpow << endl;
  }
  LinkedListIterator<SpringSvsL> i_s(springlist);
  for (;i_s; i_s++) {
    file << label(i_s->atom_type[0]) << " " << label(i_s->atom_type[1]) << endl;
    for (int sb=0; sb<i_s->forcek.get_size(); sb++) {
      file << i_s->forcek(sb);
    }
  }
}

class BondLengthData {
public:
  int atom_type[2];
  Real minl,maxl;
  BondLengthData(void) {atom_type[0]=0; atom_type[1]=0; minl=MAXFLOAT; maxl=-MAXFLOAT;}
  BondLengthData(const BondLengthData &bld) {
    atom_type[0]=bld.atom_type[0];
    atom_type[1]=bld.atom_type[1];
    minl=bld.minl;
    maxl=bld.maxl;
  }
};

void find_nn_range(LinkedList<BondLengthData> *bnd_len_data, const Structure &str_ideal, const Structure &str_relax) {
  rMatrix3d inv_cell=!str_ideal.cell;
  for (int at=0; at<str_ideal.atom_pos.get_size(); at++) {
    Real r=find_1nn_radius(str_ideal,at)+zero_tolerance;
    AtomPairIterator pair(str_ideal.cell,str_ideal.atom_pos(at),str_ideal.atom_pos);
    while (get_length_quick(pair)<=r) {
      int whichatom[2];
      BondLengthData bld;
      whichatom[0]=which_atom(str_ideal.atom_pos, pair(0), inv_cell);
      whichatom[1]=which_atom(str_ideal.atom_pos, pair(1), inv_cell);
      bld.atom_type[0]=str_ideal.atom_type(whichatom[0]);
      bld.atom_type[1]=str_ideal.atom_type(whichatom[1]);
      rVector3d ideal_offset=(pair(1)-pair(0))-(str_ideal.atom_pos(whichatom[1])-str_ideal.atom_pos(whichatom[0]));
      rVector3d cell_shift=inv_cell*ideal_offset;
      rVector3d dr=str_relax.atom_pos(whichatom[1])+str_relax.cell*cell_shift-str_relax.atom_pos(whichatom[0]);
      Real l=norm(dr);
      bld.minl=l;
      bld.maxl=l;
      LinkedListIterator<BondLengthData> bi(*bnd_len_data);
      for ( ; bi; bi++) {
        if (bld.atom_type[0]==bi->atom_type[0] && bld.atom_type[1]==bi->atom_type[1]) break;
	if (bld.atom_type[1]==bi->atom_type[0] && bld.atom_type[0]==bi->atom_type[1]) break;
      }
      if (!bi) {
        (*bnd_len_data) << new BondLengthData(bld);
      }
      else {
        bi->minl=MIN(bi->minl,l);
        bi->maxl=MAX(bi->maxl,l);
      }
      pair++;
    }
  }
}

class BondForcekData {
public:
  int atom_type[2];
  Real len;
  Array<Real> forcek;
  
  BondForcekData(void): forcek(0) {}
  BondForcekData(const BondForcekData &bfd) {
    atom_type[0]=bfd.atom_type[0];
    atom_type[1]=bfd.atom_type[1];
    len=bfd.len;
    forcek=bfd.forcek;
  }
};

void list_nn_forcek(LinkedList<BondForcekData> *p_forcek_data, const Structure &str_ideal, const Structure &str_unpert, const Structure &str_pert, const Array<Real> &forces) {
  rMatrix3d inv_cell=!str_ideal.cell;
  int pert_at;
  Real maxmove=0.;
  for (int at=0; at<str_unpert.atom_pos.get_size(); at++) {
    Real r=norm(str_pert.atom_pos(at)-str_unpert.atom_pos(at));
    if (r>maxmove) {
      maxmove=r;
      pert_at=at;
    }
  }
  rVector3d d=str_pert.atom_pos(pert_at)-str_unpert.atom_pos(pert_at);
  //  cerr << pert_at << endl;
  //  cerr << d << endl;
  Real r=find_1nn_radius(str_ideal,pert_at)+zero_tolerance;
  AtomPairIterator pair(str_ideal.cell,str_ideal.atom_pos(pert_at),str_ideal.atom_pos);
  while (get_length_quick(pair)<=r) {
    int whichatom[2];
    BondForcekData bfd;
    whichatom[0]=which_atom(str_ideal.atom_pos, pair(0), inv_cell);
    whichatom[1]=which_atom(str_ideal.atom_pos, pair(1), inv_cell);
    bfd.atom_type[0]=str_ideal.atom_type(whichatom[0]);
    bfd.atom_type[1]=str_ideal.atom_type(whichatom[1]);
    if (bfd.atom_type[0]>bfd.atom_type[1]) {
      swap(&(bfd.atom_type[0]),&(bfd.atom_type[1]));
    }
    rVector3d ideal_offset=(pair(1)-pair(0))-(str_ideal.atom_pos(whichatom[1])-str_ideal.atom_pos(whichatom[0]));
    rVector3d cell_shift=inv_cell*ideal_offset;
    // cerr << cell_shift << endl;
    rVector3d dr=str_unpert.atom_pos(whichatom[1])+str_unpert.cell*cell_shift-str_unpert.atom_pos(whichatom[0]);
    Real l=norm(dr);
    bfd.len=l;
    bfd.forcek.resize(4);
    zero_array(&(bfd.forcek));
    LinkedListIterator<BondForcekData> i_bnd(*p_forcek_data);
    for (; i_bnd; i_bnd++) {
      if (i_bnd->atom_type[0]==bfd.atom_type[0] && i_bnd->atom_type[1]==bfd.atom_type[1] && near_zero(i_bnd->len-bfd.len)) break;
    }
    if (!i_bnd) {
      (*p_forcek_data) << new BondForcekData(bfd);
    }
    Real perp_tresh=0.1;
    rVector3d u=dr/l;
    rVector3d f;
    for (int i=0; i<3; i++) {f(i)=forces(3*whichatom[1]+i);}
    Real x=(u*d);
    Real y=(u*f);
    if (2.*fabs(x)>perp_tresh*norm(d)) {
      //cerr << bfd.atom_type[0] << " " << bfd.atom_type[1] << " " << x << " " << y << endl;
      i_bnd->forcek(0)+=x*x;
      i_bnd->forcek(1)+=x*y;
    }
    rVector3d v=d-u*(u*d);
    x=norm(v);
    v.normalize();
    y=v*f; //(v*(f-u*(u*f)));
    if (2.*fabs(d*v)>perp_tresh*norm(d)) {
      i_bnd->forcek(2)+=x*x;
      i_bnd->forcek(3)+=x*y;
    }
    pair++;
  }
}

void finish_list_nn_forcek(LinkedList<BondForcekData> *p_forcek_data) {
  LinkedListIterator<BondForcekData> i_bnd(*p_forcek_data);
  for (; i_bnd; i_bnd++) {
    Array<Real> fk(2);
    for (int j=0; j<2; j++) {
      //cerr << j << " " << i_bnd->atom_type[0] << " " << i_bnd->atom_type[1] << " " << i_bnd->forcek(2*j+1) << " " << i_bnd->forcek(2*j) << endl;
      fk(j)=i_bnd->forcek(2*j+1)/i_bnd->forcek(2*j);
    }
    i_bnd->forcek=fk;
  }  
}

void list_nn_forcek_new(LinkedList<BondForcekData> *p_forcek_data, const Structure &str_ideal, const Structure &str_unpert, const Structure &str_pert) {
  rMatrix3d inv_cell=!str_ideal.cell;
  int pert_at;
  Real maxmove=0.;
  for (int at=0; at<str_unpert.atom_pos.get_size(); at++) {
    Real r=norm(str_pert.atom_pos(at)-str_unpert.atom_pos(at));
    if (r>maxmove) {
      maxmove=r;
      pert_at=at;
    }
  }
  rVector3d d=str_pert.atom_pos(pert_at)-str_unpert.atom_pos(pert_at);
  Real r=find_1nn_radius(str_ideal,pert_at)+zero_tolerance;
  AtomPairIterator pair(str_ideal.cell,str_ideal.atom_pos(pert_at),str_ideal.atom_pos);
  while (get_length_quick(pair)<=r) {
    int whichatom[2];
    BondForcekData bfd;
    whichatom[0]=which_atom(str_ideal.atom_pos, pair(0), inv_cell);
    whichatom[1]=which_atom(str_ideal.atom_pos, pair(1), inv_cell);
    bfd.atom_type[0]=str_ideal.atom_type(whichatom[0]);
    bfd.atom_type[1]=str_ideal.atom_type(whichatom[1]);
    if (bfd.atom_type[0]>bfd.atom_type[1]) {
      swap(&(bfd.atom_type[0]),&(bfd.atom_type[1]));
    }
    rVector3d ideal_offset=(pair(1)-pair(0))-(str_ideal.atom_pos(whichatom[1])-str_ideal.atom_pos(whichatom[0]));
    rVector3d cell_shift=inv_cell*ideal_offset;
    rVector3d dr=str_unpert.atom_pos(whichatom[1])+str_unpert.cell*cell_shift-str_unpert.atom_pos(whichatom[0]);
    bfd.len=norm(dr);
    bfd.forcek.resize(2);
    zero_array(&(bfd.forcek));
    LinkedListIterator<BondForcekData> i_bnd(*p_forcek_data);
    for (; i_bnd; i_bnd++) {
      if (i_bnd->atom_type[0]==bfd.atom_type[0] && i_bnd->atom_type[1]==bfd.atom_type[1] && near_zero(i_bnd->len-bfd.len)) break;
    }
    if (!i_bnd) {
      (*p_forcek_data) << new BondForcekData(bfd);
    }
    pair++;
  }
}

void calc_equation_forcek(Array2d<Real> *eqn, LinkedList<BondForcekData> *p_forcek_data, const Structure &str_ideal, const Structure &str_unpert, const Structure &str_pert, Real tol) {
  eqn->resize(iVector2d(3*str_ideal.atom_pos.get_size(),p_forcek_data->get_size()*2));
  zero_array(eqn);

  rMatrix3d inv_cell=!str_ideal.cell;
  int pert_at;
  Real maxmove=0.;
  for (int at=0; at<str_unpert.atom_pos.get_size(); at++) {
    Real r=norm(str_pert.atom_pos(at)-str_unpert.atom_pos(at));
    if (r>maxmove) {
      maxmove=r;
      pert_at=at;
    }
  }
  rVector3d d=str_pert.atom_pos(pert_at)-str_unpert.atom_pos(pert_at);
  Real r=find_1nn_radius(str_ideal,pert_at)+zero_tolerance;
  AtomPairIterator pair(str_ideal.cell,str_ideal.atom_pos(pert_at),str_ideal.atom_pos);
  while (get_length_quick(pair)<=r) {
    int whichatom[2];
    BondForcekData bfd;
    whichatom[0]=which_atom(str_ideal.atom_pos, pair(0), inv_cell);
    whichatom[1]=which_atom(str_ideal.atom_pos, pair(1), inv_cell);
    bfd.atom_type[0]=str_ideal.atom_type(whichatom[0]);
    bfd.atom_type[1]=str_ideal.atom_type(whichatom[1]);
    if (bfd.atom_type[0]>bfd.atom_type[1]) {
      swap(&(bfd.atom_type[0]),&(bfd.atom_type[1]));
    }
    rVector3d ideal_offset=(pair(1)-pair(0))-(str_ideal.atom_pos(whichatom[1])-str_ideal.atom_pos(whichatom[0]));
    rVector3d cell_shift=inv_cell*ideal_offset;
    rVector3d dr=str_unpert.atom_pos(whichatom[1])+str_unpert.cell*cell_shift-str_unpert.atom_pos(whichatom[0]);
    bfd.len=norm(dr);
    LinkedListIterator<BondForcekData> i_bnd(*p_forcek_data);
    int idx=0;
    for (; i_bnd; i_bnd++, idx++) {
      if (i_bnd->atom_type[0]==bfd.atom_type[0] && i_bnd->atom_type[1]==bfd.atom_type[1] && near_zero(i_bnd->len-bfd.len)) break;
    }

    rVector3d u=dr;
    u.normalize();
    u=u*(u*d);
    rVector3d v=d-u;

    if (norm(u)>tol) {
      for (int i=0; i<3; i++) {
	(*eqn)(whichatom[0]*3+i,idx*2)  +=-u(i);
	(*eqn)(whichatom[1]*3+i,idx*2)  += u(i);
      }
    }
    if (norm(v)>tol) {
     for (int i=0; i<3; i++) {
	(*eqn)(whichatom[0]*3+i,idx*2+1)+=-v(i);
	(*eqn)(whichatom[1]*3+i,idx*2+1)+= v(i);
      }
    }
    pair++;
  }
}

void finish_list_nn_forcek(LinkedList<BondForcekData> *p_forcek_data, const LinkedList<Structure> &ideal_str_list,
			   const LinkedList<Structure> &unpert_str_list,
			   const LinkedList<Structure> &pert_str_list, const Array<Real> &forces, Real tol) {
  Array2d<Real> big_eqn(iVector2d(forces.get_size(),p_forcek_data->get_size()*2));
  zero_array(&big_eqn);
  LinkedListIterator<Structure> i_ideal_str(ideal_str_list);
  LinkedListIterator<Structure> i_unpert_str(unpert_str_list);
  LinkedListIterator<Structure> i_pert_str(pert_str_list);
  int r=0;
  for (; i_ideal_str; i_ideal_str++, i_unpert_str++, i_pert_str++) {
    Array2d<Real> eqn;
    calc_equation_forcek(&eqn,p_forcek_data,*i_ideal_str,*i_unpert_str,*i_pert_str,zero_tolerance);
    for (int j=0; j<eqn.get_size()(0); j++, r++) {
      for (int i=0; i<eqn.get_size()(1); i++) {
	big_eqn(r,i)=eqn(j,i);
      }
    }
  }
  Array<Real> beta;
  Array<Real> w;
  calc_ols(&beta,big_eqn,forces,w,0);
  Array2d<Real> varbeta;
  inner_product(&varbeta,big_eqn,big_eqn);
  //invert_matrix(&varbeta,varbeta);
  calc_ols_var(&varbeta,big_eqn,forces);
  LinkedListIterator<BondForcekData> ibnd(*p_forcek_data);
  for (int i=0; ibnd; ibnd++, i+=2) {
    if (varbeta(i,i)>tol) {
      ibnd->forcek(0)=0;
    }
    else {
      ibnd->forcek(0)=beta(i);
    }
    if (varbeta(i+1,i+1)>tol) {
      ibnd->forcek(1)=0;
    }
    else {
      ibnd->forcek(1)=beta(i+1);
    }
  }
  //  cerr << big_eqn << endl;
  //  cerr << forces << endl;
}

extern char *helpstring;

int main(int argc, char *argv[]) {
  //  cerr.setf(ios::fixed | ios::showpos);
  cerr.precision(5);
  char *latfilename="lat.in";
  char *strnamefilename="strname.in";
  Real enclosed_radius=0.;
  Real displ_mag=0.2;
  Real max_strain=0.01;
  int nb_vol=2;
  int fitfc=0;
  int redun=0;
  int maxlpow=1;
  int dirdep=0;
  int concpow=0;
  Real maxfklen=MAXFLOAT;
  char *strainfilename="";
  Real fkeqntol=2.0;
  int dohelp=0;
  int sigdig=5;
  int olddirstruct=0;
  AskStruct options[]={
    {"","Fit Stiffness VS Length transferable force constants " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-f","Fit force constants (otherwise, generate pertubations)",BOOLVAL,&fitfc},
    {"-l","Input file defining the lattice (default: lat.in)",STRINGVAL,&latfilename},
    {"-dn","Input file listing the directories containing the structures used to calculate force constants (default: strname.in)",STRINGVAL,&strnamefilename},
    {"-er","Minimum distance between displaced atoms",REALVAL,&enclosed_radius},
    {"-dr","Displacement of the perturbed atom (default: 0.2)",REALVAL,&displ_mag},
    {"-ms","Strain of the maximum volume sampled (default: 0.01)",REALVAL,&max_strain},
    {"-ns","Number of volume sample (default: 2)",INTVAL,&nb_vol},
    {"-op","Order of the polynomial used to fit stiffness vs length",INTVAL,&maxlpow},
    {"-dd","Use direction-dependent force constants",BOOLVAL,&dirdep},
    {"-pc","Maximum power of composition used in fit (default: no composition dependence)",INTVAL,&concpow},
    {"-msl","Maximum spring length",REALVAL,&maxfklen},
    {"-sf","Extra strain file",STRINGVAL,&strainfilename},
    {"-eqt","Tolerance for excluding force constants in plots.",REALVAL,&fkeqntol},
    {"-sig","Number of significant digits printed (default: 5)",INTVAL,&sigdig},
    {"-z","Tolerance for finding symmetry operations (default: 1e-3)",REALVAL,&zero_tolerance},
    {"-ncs","No check for singular matrix in fit",BOOLVAL,&redun},
    {"-od","Use older 1-level only directory tree",BOOLVAL,&olddirstruct},
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

  if (enclosed_radius==0 && !fitfc) {
    ERRORQUIT("You need to specify at least the -er parameter to generator perturbations or the -f option to fit the force constants.");
  }

  rMatrix3d extrastrain;
  extrastrain.identity();
  if (strlen(strainfilename)>0) {
    ifstream strainfile(strainfilename);
    strainfile >> extrastrain;
  }
  Structure lat;
  Array<Arrayint> site_type_list;
  Array<AutoString> atom_label;
  rMatrix3d axes;
  {
    ifstream file(latfilename);
    if (!file) ERRORQUIT("Unable to open lattice file.");
    parse_lattice_file(&lat.cell, &lat.atom_pos, &lat.atom_type, &site_type_list, &atom_label, file, &axes);
    if (fabs(det(lat.cell))<zero_tolerance) ERRORQUIT("Lattice vectors are coplanar.");
    wrap_inside_cell(&lat.atom_pos,lat.atom_pos,lat.cell);
  }

  LinkedList<SpringSvsL> springsvsl_list;
  rMatrix3d ilatcell=!(lat.cell);
  for (int at=0; at< lat.atom_type.get_size(); at++) {
    Real r1nn=find_1nn_radius(lat, at)+zero_tolerance;
    AtomPairIterator pair(lat.cell,lat.atom_pos(at),lat.atom_pos);
    for (; norm(pair(1)-pair(0))<=r1nn; pair++) {
      Array<int> nb_atom_type(2);
      Array<int> index_site_type_list(2);
      for (int endpt=0; endpt<2; endpt++) {
        index_site_type_list(endpt)=lat.atom_type(which_atom(lat.atom_pos,pair(endpt),ilatcell));
        nb_atom_type(endpt)=site_type_list(index_site_type_list(endpt)).get_size();
      }
      MultiDimIterator<Array<int> > pairtype(nb_atom_type);
      for (; pairtype; pairtype++) {
        SpringSvsL spr;
        for (int endpt=0; endpt<2; endpt++) {
          spr.atom_type[endpt]=site_type_list(index_site_type_list(endpt))(((Array<int> &)pairtype)(endpt));
        }
	if (spr.atom_type[0]>spr.atom_type[1]) {
	  swap(&(spr.atom_type[0]),&(spr.atom_type[1]));
	}
	if (strcmp(atom_label(spr.atom_type[0]),"Vac")!=0 && strcmp(atom_label(spr.atom_type[0]),"Vac")!=0) {
	  add_unique(&springsvsl_list,spr,SameEndAtom());
	}
      }
    }
  }

  LinkedList<rMatrix3d> dirdep_mat;
  if (!dirdep) {
    rMatrix3d *pid=new rMatrix3d;
    pid->identity();
    dirdep_mat << pid;
  }
  else {
    SpaceGroup spacegroup;
    find_spacegroup(&spacegroup.point_op,&spacegroup.trans,lat.cell,lat.atom_pos,lat.atom_type);
    Array<rMatrix3d> pointgroup;
    pointgroup_from_spacegroup(&pointgroup, spacegroup.point_op);
    find_sym_const_tensor(&dirdep_mat, pointgroup,1);
  }

  Array2d<Real> redtonored;
  {
    Array2d<Real> m;
    Array<Real> v;
    calc_noredtored(&m,&v,lat,site_type_list);
    calc_redtonored(&redtonored, m,v);
  }

  int ddms=dirdep_mat.get_size();
  Array<Real> maxpows(1+ddms+redtonored.get_size()(0));
  Array<Real> dlineq(maxpows.get_size()),duineq(maxpows.get_size());
  Array<Real> cineq(maxpows.get_size());
  zero_array(&dlineq);
  zero_array(&duineq);
  zero_array(&cineq);
  maxpows(0)=maxlpow;
  for (int i=0; i<dirdep_mat.get_size(); i++) {
    maxpows(1+i)=1;
    dlineq(1+i)=-1;
    duineq(1+i)=1;
  }
  for (int i=0; i<redtonored.get_size()(0); i++) {
    maxpows(1+ddms+i)=concpow;
    cineq(1+ddms+i)=1;
  }
  LinkedList<LinearInequality> linineq;
  linineq << new LinearInequality(dlineq,-1) << new LinearInequality(duineq,1) << new LinearInequality(cineq,concpow);
  MultiDimPoly multipoly(maxpows,linineq);

  if (!fitfc) {
    LinkedList<Structure> str_list;
    LinkedList<Structure> rel_str_list;
    LinkedList<AutoString> dir_list;
    {
      ifstream strnamefile(strnamefilename);
      if (!strnamefile) {ERRORQUIT("Unable to open structure name file");}
      while (skip_delim(strnamefile)) {
        AutoString dirname;
        get_string(&dirname,strnamefile);
        dir_list << new AutoString(dirname);
        Structure str;
        {
          ostrstream ostrfilename;
          ostrfilename << dirname << '/' << "str.out" << '\0';
          ifstream strfile(ostrfilename.str());
          if (!strfile) {
            cerr << "Unable to open " << ostrfilename.str() << endl;
            ERRORQUIT("Aborting");
          }
          if (!parse_structure_file(&str.cell,&str.atom_pos,&str.atom_type,atom_label,strfile)) {
	    cerr << "Error in file " << dirname << '/' << "str.out" << endl;
	    ERRORQUIT("Aborting");
	  }
          wrap_inside_cell(&str.atom_pos,str.atom_pos,str.cell);
          str_list << new Structure(str);
	  remove_vacancies(&str,atom_label);
        }
        Structure rel_str;
        {
          ostrstream rel_strfilename;
          rel_strfilename << dirname << '/' << "str_relax.out" << '\0';
          ifstream rel_strfile(rel_strfilename.str());
          if (!rel_strfile) {
            cerr << "Unable to open " << rel_strfilename.str() << endl;
            ERRORQUIT("Aborting");
          }
          if (!parse_structure_file(&rel_str.cell,&rel_str.atom_pos,&rel_str.atom_type,atom_label,rel_strfile)) {
	    cerr << "Error in file " << dirname << '/' << "str_relax.out" << endl;
	    ERRORQUIT("Aborting");
	  }
          wrap_inside_cell(&rel_str.atom_pos,rel_str.atom_pos,rel_str.cell);
          reorder_atoms(&rel_str,str);
          rel_str_list << new Structure(rel_str);
        }
      }
    }
  
    LinkedListIterator<AutoString> i_dir(dir_list);
    LinkedListIterator<Structure> i_str(str_list);
    LinkedListIterator<Structure> i_rel_str(rel_str_list);
    cerr << "Generating perturbations..." << endl;
    for ( ; i_str; i_dir++, i_str++, i_rel_str++) {
      char thecwd[MAX_LINE_LEN];
      if (!getcwd(thecwd,MAX_LINE_LEN)) {ERRORQUIT("Pathname too long, Increase MAX_LINE_LEN.");}
      chdir_robust(*i_dir);
      cerr << *i_dir << endl;
      SpaceGroup str_spacegroup;
      str_spacegroup.cell=i_rel_str->cell;
      find_spacegroup(&str_spacegroup.point_op,&str_spacegroup.trans,i_rel_str->cell,i_rel_str->atom_pos,i_rel_str->atom_type);
      Array<int> pert_atom;
      Array<rVector3d> pert_direction;
      find_pertubations(&pert_atom, &pert_direction,*i_rel_str,str_spacegroup);
      for (int iv=0; iv<nb_vol; iv++) {
        Real strain=max_strain*(Real)iv/(Real)MAX(nb_vol-1,1);
	if (!olddirstruct) {
	  ostrstream volname;
	  volname << "vol_" << strain*100 << '\0';
	  cerr << "  " << volname.str() << endl;
	  mkdir(volname.str(),S_IRWXU | S_IRWXG | S_IRWXO);
	  chdir_robust(volname.str());
	  Structure s_str_relax;
	  stretch_str(&s_str_relax,*i_rel_str,strain);
	  rMatrix3d newaxes=(1.+strain)*(i_rel_str->cell)*(!i_str->cell)*axes;
	  {
	    ofstream file("str_relax.out");
	    write_structure(s_str_relax, atom_label, newaxes, file);
	  }
	}
        for (int p=0; p<pert_atom.get_size(); p++) {
	  ostrstream pertname;
	  if (olddirstruct) {
	    pertname << displ_mag << "_" << strain*100 << "_" << p << '\0';
	  }
	  else {
            pertname << "p+" << displ_mag << "_" << enclosed_radius << "_" << p << '\0';
	    cerr << "    " << pertname.str() << endl;
	  }
	  mkdir(pertname.str(),S_IRWXU | S_IRWXG | S_IRWXO);
	  chdir_robust(pertname.str());
          rMatrix3d ideal_supercell;
          find_smallest_supercell_enclosing_sphere(&ideal_supercell,i_str->cell,enclosed_radius);
          rVector3d zerovect(0.,0.,0.);
          Structure ideal_str;
          calc_pert_str(&ideal_str,ideal_supercell,*i_str,zerovect,zerovect,0.);
          {
            ofstream idealstrfile("str_ideal.out");
	    idealstrfile.setf(ios::fixed);
	    idealstrfile.precision(sigdig);
            write_structure(ideal_str, atom_label, axes, idealstrfile);
          }
          Structure unpert_str;
          rMatrix3d supercell=(i_rel_str->cell)*((!i_str->cell)*ideal_supercell);
          calc_pert_str(&unpert_str,supercell,*i_rel_str,zerovect,zerovect,strain);
	  rMatrix3d newaxes=(1.+strain)*(i_rel_str->cell)*(!i_str->cell)*axes;
          {
            ofstream unpertstrfile("str_unpert.out");
	    unpertstrfile.setf(ios::fixed);
	    unpertstrfile.precision(sigdig);
            write_structure(unpert_str, atom_label, newaxes, unpertstrfile);
          }
          Structure pert_str;
          calc_pert_str(&pert_str,supercell,*i_rel_str,i_rel_str->atom_pos(pert_atom(p)),displ_mag*pert_direction(p),strain);
          {
            ofstream pertstrfile("str.out");
	    pertstrfile.setf(ios::fixed);
	    pertstrfile.precision(sigdig);
            write_structure(pert_str, atom_label, newaxes, pertstrfile);
          }
          if (!file_exists("force.out")) {
            ofstream waitfile("wait");
          }
          chdir("..");
        }
	if (!olddirstruct) {
	  chdir("..");
	}
      }
      chdir(thecwd);
    }
  }
  else { // fit forcek instead of generate pert;
    LinkedList<Structure> ideal_str_list;
    LinkedList<Structure> pert_str_list;
    LinkedList<Structure> unpert_str_list;
    LinkedList<Array<Real> > forcevector_list;
    int nb_row=0;
    system("ls */force.out */*/force.out */*/*/force.out 2> /dev/null | sed 's+/force.out++g' > pertlist.out");
    ifstream pertfile("pertlist.out");
    if (!pertfile) {ERRORQUIT("Unable to open pertlist.out.");}
    int strcnt=0;
    while (skip_delim(pertfile)) {
      AutoString dirname;
      get_string(&dirname,pertfile);
      char thecwd[MAX_LINE_LEN];
      if (!getcwd(thecwd,MAX_LINE_LEN)) {ERRORQUIT("Pathname too long, Increase MAX_LINE_LEN.");}
      chdir_robust(dirname);
      if (file_exists("str_unpert.out")) {
	cerr << strcnt << ") Reading " << dirname << endl;
	strcnt++;
	Array<int> copyfrom;
	Array<iVector3d> cellshift;
	Structure *pideal_str=new Structure;
	{
	  ifstream file("str_ideal.out");
	  if (!file) {
	    cerr << "Unable to open str_ideal.out in directory " << dirname << endl;
	    ERRORQUIT("Aborting");
	  }
	  parse_structure_file(&(pideal_str->cell),&(pideal_str->atom_pos),&(pideal_str->atom_type),atom_label,file);
	  if (pideal_str->atom_pos.get_size()==0) {ERRORQUIT("Unable to read str_ideal.out");}
	  pideal_str->cell=extrastrain*(pideal_str->cell);
	  for (int at=0; at<pideal_str->atom_pos(); at++) {
	    pideal_str->atom_pos(at)=extrastrain*(pideal_str->atom_pos(at));
	  }
	  ideal_str_list << pideal_str;
	}
	Structure *punpert_str=new Structure;
	{
	  ifstream file("str_unpert.out");
	  if (!file) {
	    cerr << "Unable to open str_unpert.out in directory " << dirname << endl;
	    ERRORQUIT("Aborting");
	  }
	  parse_structure_file(&(punpert_str->cell),&(punpert_str->atom_pos),&(punpert_str->atom_type),atom_label,file);
	  if (punpert_str->atom_pos.get_size()==0) {ERRORQUIT("Unable to read str_unpert.out");}
	  Structure ideal_str_novac=*pideal_str;
	  remove_vacancies(&ideal_str_novac,atom_label);
	  reorder_atoms(punpert_str,ideal_str_novac);
	  unpert_str_list << punpert_str;
	}
	Structure *ppert_str=new Structure;
	{
	  ifstream file("str_relax.out");
	  if (!file) {
	    cerr << "Unable to open str_relax.out in directory " << dirname << endl;
	    ERRORQUIT("Aborting");
	  }
	  parse_structure_file(&(ppert_str->cell),&(ppert_str->atom_pos),&(ppert_str->atom_type),atom_label,file);
	  if (ppert_str->atom_pos.get_size()==0) {ERRORQUIT("Unable to read str_relax.out");}
	  reorder_atoms(ppert_str,*punpert_str,&copyfrom,&cellshift);
	  pert_str_list << ppert_str;
	}
	{
	  ifstream forcefile("force.out");
	  if (!forcefile) {
	    cerr << "Unable to open force.out in directory " << dirname << endl;
	    ERRORQUIT("Aborting");
	  }
	  Array<Real> *pforce=new Array<Real>;
	  read_force_vector(pforce,forcefile,copyfrom);
	  {
	    ofstream tmpf("tmp3.out");
	    for (int i=0; i<pforce->get_size()/3; i++) {
	      tmpf << (*pforce)(3*i+0) << " " << (*pforce)(3*i+1) << " "<< (*pforce)(3*i+2) << endl;
	    }
	  }
	  forcevector_list << pforce;
	  nb_row+=pforce->get_size();
	}
      }
      chdir(thecwd);
    }

    int nb_col=0;
    {
      LinkedListIterator<SpringSvsL> i_s(springsvsl_list);
      for (;i_s; i_s++) {
	i_s->forcek.resize(2);
	i_s->forcek(0).resize(multipoly.get_dim_param());
	i_s->forcek(1).resize(multipoly.get_dim_param());
	nb_col+=i_s->forcek(0).get_size() + i_s->forcek(1).get_size();
      }
    }
    
    Array2d<Real> big_eqnmat(nb_row,nb_col);
    Array<Real> big_force(nb_row);
    LinkedListIterator<Structure> i_ideal_str(ideal_str_list);
    LinkedListIterator<Structure> i_pert_str(pert_str_list);
    LinkedListIterator<Structure> i_unpert_str(unpert_str_list);
    LinkedListIterator<Array<Real> > i_forcevector(forcevector_list);
    LinkedList<BondForcekData> raw_forcek_data;
    strcnt=0;
    for (int r=0; i_ideal_str; i_ideal_str++, i_pert_str++, i_unpert_str++, i_forcevector++, strcnt++) {
      //cerr << "r= " << r << endl;
      Array<Real> conc,red_conc;
      //      if (!calc_concentration(&red_conc, lat,site_type_list,*i_ideal_str)) {
      //	cerr << "Structure " << strcnt << endl;
      //	ERRORQUIT("Problem calculating concentration");
      //      }
      calc_concentration(&red_conc, i_ideal_str->atom_type,atom_label.get_size());
      product(&conc,redtonored,red_conc);
      remove_vacancies(&(*i_ideal_str),atom_label);
      Array<Array<Real> > eqn;
      calc_equation_matrix(&eqn,*i_ideal_str,*i_unpert_str,*i_pert_str,&springsvsl_list,dirdep_mat,conc,multipoly,maxfklen);
      list_nn_forcek_new(&raw_forcek_data, *i_ideal_str, *i_unpert_str, *i_pert_str);
      // cerr << raw_forcek_data.get_size() << endl;
      for (int i=0; i<i_forcevector->get_size(); i++) {
        big_force(r)=(*i_forcevector)(i);
        for (int j=0; j<nb_col; j++) {
          big_eqnmat(r,j)=eqn(i)(j);
        }
        r++;
      }
    }
    //cerr << "fin1" << endl;
    finish_list_nn_forcek(&raw_forcek_data,ideal_str_list,unpert_str_list,pert_str_list,big_force,sqr(fkeqntol));
    Array<Real> beta;
    Array<Real> w;
    //cerr << "ols" << endl;
    calc_ols(&beta,big_eqnmat,big_force,w,!redun);

    Array<Real> force_hat;
    //cerr << "prod" << endl;
    product(&force_hat,big_eqnmat,beta);
    //cerr << "write" << endl;
    {
      ofstream logfile("fitsvsl.log");
      logfile.setf(ios::fixed | ios::showpos);
      logfile.precision(sigdig);
      logfile << "Equation Matrix" << endl;
      logfile << big_eqnmat << endl;
      logfile << endl;
      logfile << "actual force   predicted force   difference" << endl;
      for (int i=0; i<big_force.get_size(); i++) {
	logfile << big_force(i) << " " << force_hat(i) << " " << (big_force(i)-force_hat(i)) << endl;
      }
    }

    LinkedListIterator<SpringSvsL> i_s(springsvsl_list);
    for (int c=0;i_s; i_s++) {
      for (int sb=0; sb<i_s->forcek.get_size(); sb++) {
        for (int p=0; p<i_s->forcek(sb).get_size(); p++) {
          i_s->forcek(sb)(p)=beta(c);
          c++;
        }
      }
    }
    {
      ofstream spr_file("slspring.out");
      spr_file.setf(ios::fixed);
      spr_file.precision(sigdig);
      write_spring_file(spr_file, springsvsl_list, atom_label, maxlpow, dirdep, concpow);
    }

    {
      ofstream file("fitsvsl.gnu");
      LinkedList<BondLengthData> bld;
      LinkedListIterator<Structure> i_ideal_str(ideal_str_list);
      LinkedListIterator<Structure> i_unpert_str(unpert_str_list);
      for ( ; i_ideal_str; i_ideal_str++, i_unpert_str++) {
	find_nn_range(&bld, *i_ideal_str, *i_unpert_str);
      }
      LinkedListIterator<SpringSvsL> i_s(springsvsl_list);
      for (;i_s; i_s++) {
	LinkedListIterator<BondLengthData> i_bld(bld);
	for (;i_bld; i_bld++) {
	  if (i_s->atom_type[0]==i_bld->atom_type[0] && i_s->atom_type[1]==i_bld->atom_type[1] ) break;
	  if (i_s->atom_type[1]==i_bld->atom_type[0] && i_s->atom_type[0]==i_bld->atom_type[1] ) break;
	}
	if (i_bld) {
	  ostrstream line;
	  line << "f_" << atom_label(i_s->atom_type[0]) << "-" << atom_label(i_s->atom_type[1]) << ".dat" << '\0';
	  unlink(line.str());
	  file << "set title '" << atom_label(i_s->atom_type[0]) << "-" << atom_label(i_s->atom_type[1]) << "'" << endl;
	  file << "plot [" << i_bld->minl << ":" << i_bld->maxl << "] ";
	  for (int sb=0; sb<i_s->forcek.get_size(); sb++) {
	    if (dirdep==0 && concpow==0) {
	      for (int p=0; p<i_s->forcek(sb).get_size(); p++) {
		file << "(" << i_s->forcek(sb)(p) << ")";
		switch(p) {
		case 0:
		  break;
		case 1:
		  file << "*x";
		  break;
		default:
		  file << "*x**" << p;
		}
		if (p<i_s->forcek(sb).get_size()-1) file << "+";
	      }
	      file << " t '" << (sb==0 ? 's' : 'b') << "', ";
	    }
	    file << "'" << line.str() << "' u 1:" << (sb==0 ? '2' : '3');
	    if (sb<i_s->forcek.get_size()-1) file << ",";
	  }
	}
	file << endl;
	file << "pause -1" << endl;
      }
    }
    {
      LinkedListIterator<BondForcekData> i_bfd(raw_forcek_data);
      for (; i_bfd; i_bfd++) {
	ostrstream line;
	line << "f_" << atom_label(i_bfd->atom_type[0]) << "-" << atom_label(i_bfd->atom_type[1]) << ".dat" << '\0';
	{
	  ofstream datfile(line.str(),ios::app);
	  datfile << i_bfd->len << " " << i_bfd->forcek(0) << " " << i_bfd->forcek(1) << endl;
	}
      }
    }
  }
}
