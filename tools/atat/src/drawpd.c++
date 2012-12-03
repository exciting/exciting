#include "drawpd.h"

MF_ThermoData::MF_ThermoData(void): PhaseThermoData() {}

MF_ThermoData::MF_ThermoData(  const Structure &lattice, const SpaceGroup &space_group,
                  const Structure &str,
                  const LinkedList<ArrayrVector3d> &cluster_list, const LinkedList<Real> &eci_list): PhaseThermoData() {
    init(lattice,space_group,str,cluster_list,eci_list);
  }

void MF_ThermoData::init(const Structure &lattice, const SpaceGroup &space_group,
                  const Structure &str,
                  const LinkedList<ArrayrVector3d> &cluster_list, const LinkedList<Real> &eci_list) {
  cur_phase=this;

  volume=str.atom_pos.get_size();
  max_flip=2.0;

  E_ref=find_empty_eci(cluster_list,eci_list)/(Real)(lattice.atom_pos.get_size());
  ideal_spin.resize(str.atom_type.get_size());
  for (int i=0; i<ideal_spin.get_size(); i++) {ideal_spin(i)=(Real)str.atom_type(i);}
  clusters.resize(str.atom_pos.get_size());
  eci.resize(str.atom_pos.get_size());
  rMatrix3d inv_cell=!str.cell;
  for (int s=0; s<str.atom_pos.get_size(); s++) {
    Array<ArrayrVector3d> o_cluster;
    Array<Real> o_eci;
    find_clusters_overlapping_site(&o_cluster,&o_eci,
                  str.atom_pos(s),lattice,space_group,cluster_list,eci_list,1);
    clusters(s).resize(o_cluster.get_size());
    eci(s).resize(o_eci.get_size());
    for (int c=0; c<o_cluster.get_size(); c++) {
      clusters(s)(c).resize(o_cluster(c).get_size()-1);
      eci(s)(c)=o_eci(c);
      for (int i=1; i<o_cluster(c).get_size(); i++) {
        clusters(s)(c)(i-1)=which_atom(str.atom_pos,o_cluster(c)(i),inv_cell);
      }
    }
  }
}

  void MF_ThermoData::set_T_mu(Real T, Real mu) {
    Array<Real> spin=ideal_spin;
    Array<Real> dE(spin.get_size());
    Real max_diff=MAXFLOAT;
    int iter=0;
    while (max_diff>spin_precision && iter<10000) {
      iter++;
      for (int s=0; s<spin.get_size(); s++) {
        dE(s)=-mu;
        for (int c=0; c<eci(s).get_size(); c++) {
          Real spin_product=1.;
          for (int i=0; i<clusters(s)(c).get_size(); i++) {
            spin_product*=spin(clusters(s)(c)(i));
          }
          dE(s)+=eci(s)(c)*spin_product;
        }
      }
      max_diff=0.;
      for (int s=0; s<spin.get_size(); s++) {
        Real new_spin;
	if (-2.*dE(s) > 1e2*T) {
	  new_spin=spin(s);
	}
	else {
	  new_spin=(-1.+1.*exp(-2.*dE(s)/T))/(1.+exp(-2.*dE(s)/T));
	}
        new_spin=(new_spin+spin(s))/2.;
        max_diff=MAX(max_diff,fabs(new_spin-spin(s)));
        spin(s)=new_spin;
      }
    }
    Real flip=0.;
    for (int i=0; i<spin.get_size(); i++) {
      flip=MAX(flip,fabs(ideal_spin(i)-spin(i)));
    }
    if (flip>max_flip) {
      cur_x=0;
      for (int s=0; s<spin.get_size(); s++) {
        cur_x+=spin(s);
      }
      cur_x/=volume;
      cur_phi=MAXFLOAT;
    }
    else {
      cur_x=0.;
      cur_phi=0.;
      for (int s=0; s<spin.get_size(); s++) {
        cur_x+=spin(s);
	cur_phi+=-mu*spin(s);
	for (int c=0; c<clusters(s).get_size(); c++) {
	  Real prod=spin(s);
	  for (int i=0; i<clusters(s)(c).get_size(); i++) {
	    prod*=spin(clusters(s)(c)(i));
	  }
	  cur_phi+=eci(s)(c)*prod/(clusters(s)(c).get_size()+1);
	}
      }
      if (T>0) {
	for (int s=0; s<spin.get_size(); s++) {
	  cur_phi+=-T*log(exp(-(1-spin(s))*dE(s)/T)+exp(-(-1-spin(s))*dE(s)/T));
	}
      }
      cur_x/=volume;
      cur_phi=cur_phi/volume+E_ref;
    }
  }

  void MF_ThermoData::read(istream &file) {
    file >> volume;
    file >> ideal_spin;
    file >> clusters;
    file >> eci;
    file >> max_flip;
    file >> E_ref;
  }

  void MF_ThermoData::write(ostream &file) {
    file << object_label() << endl;
    file << volume << endl;
    file << ideal_spin << endl;
    file << clusters << endl;
    file << eci << endl;
    file << max_flip << endl;
    file << E_ref << endl;
  }

Real MF_ThermoData::spin_precision=1e-3;

LTE_ThermoData::LTE_ThermoData(void): PhaseThermoData() {}

LTE_ThermoData::LTE_ThermoData(const Structure &lattice,
			       const SpaceGroup &space_group,
			       const Structure &str,
			       const LinkedList<ArrayrVector3d> &cluster_list,
			       const LinkedList<Real> &eci_list): PhaseThermoData() {
    init(lattice,space_group,str,cluster_list,eci_list);
}

void LTE_ThermoData::init(const Structure &lattice,
		     const SpaceGroup &space_group,
		     const Structure &str,
		     const LinkedList<ArrayrVector3d> &cluster_list,
		     const LinkedList<Real> &eci_list) {
  cur_phase=this;

  x0=0.;
  volume=str.atom_pos.get_size();
  E0=volume*find_empty_eci(cluster_list,eci_list)/(Real)(lattice.atom_pos.get_size());

  dE.resize(str.atom_pos.get_size());
  dx.resize(str.atom_pos.get_size());
  zero_array(&(dE));
  zero_array(&(dx));

  rMatrix3d inv_cell=!str.cell;
  for (int s=0; s<str.atom_pos.get_size(); s++) {
    dx(s)=-2*str.atom_type(s);
    x0+=str.atom_type(s);
    Array<ArrayrVector3d> o_cluster;
    Array<Real> o_eci;
    find_clusters_overlapping_site(&o_cluster,&o_eci,
                  str.atom_pos(s),lattice,space_group,cluster_list,eci_list,1);
    for (int c=0; c<o_cluster.get_size(); c++) {
      int prod=str.atom_type(s);
      for (int i=1; i<o_cluster(c).get_size(); i++) {
        int atom=which_atom(str.atom_pos,o_cluster(c)(i),inv_cell);
        prod*=str.atom_type(atom);
      }
      dE(s)+=-2.*(Real)prod*o_eci(c);
      E0+=(Real)prod*o_eci(c)/(Real)(o_cluster(c).get_size());
    }
  }
}

  void LTE_ThermoData::set_T_mu(Real T, Real mu) {
    if (T==0) {
      cur_phi=(E0-mu*x0)/volume;
      cur_x=x0/volume;
      cur_E=E0/volume;
    }
    else {
      cur_phi=MAXFLOAT;
      cur_x=x0/volume;
      cur_E=E0/volume;
      Real num=x0;
      Real numE=0.;
      Real den=0.;
      for (int i=0; i<dE.get_size(); i++) {
	Real dF=dE(i)-mu*dx(i);
	if (T>dF) return;
	if (dF/T < 20) {
	  Real b=exp(-dF/T);
	  num+=dx(i)*b;
	  numE+=dF*b;
	  den+=b;
	}
      }
      cur_phi=(E0-mu*x0-T*den)/volume;
      cur_x=num/volume;
      cur_E=(E0-mu*x0+numE)/volume;
    }
  }

  Real LTE_ThermoData::E(Real T, Real mu) {
    set_T_mu(T,mu);
    return cur_E;
  }

  void LTE_ThermoData::read(istream &file) {
    file >> E0;
    file >> x0;
    file >> volume;
    file >> dE;
    file >> dx;
  }

  void LTE_ThermoData::write(ostream &file) {
    file << object_label() << endl;
    file << E0 << endl;
    file << x0 << endl;
    file << volume << endl;
    file << dE << endl;
    file << dx << endl;
  }

  Real HTE_ThermoData::F_of_x(Real T, Real x) {
    Real e=0;
    for (int i=0; i<poly_x.get_size(); i++) {
      e+=poly_x(i)*ipow(x,i);
    }
    if (x<-1.+zero_tolerance || x>1.-zero_tolerance) {
      return e;
    }
    return e+T*((1+x)*log(1+x)+(1-x)*log(1-x)-2.*log(2.))/2.;
  }

  Real HTE_ThermoData::mu_of_x(Real T, Real x) {
    if ((1.-fabs(x))<zero_tolerance) return (x<0 ? -MAXFLOAT : MAXFLOAT);
    Real de=0;
    for (int i=1; i<poly_x.get_size(); i++) {
      de+=poly_x(i)*ipow(x,i-1)*(Real)i;
    }
    return de+T*log((1+x)/(1-x))/2.;
  }

  Real HTE_ThermoData::x_of_mu(Real T, Real mu) {
    Real x=0.;
    Real oldx=0.;
    Real side=( mu_of_x(T,x)>mu ? -1 : 1);
    while ( side*(mu-mu_of_x(T,x)) > 0. ) {
      oldx=x;
      x=(side+x)/2.;
    }
    Real min_x=MIN(x,oldx);
    Real max_x=MAX(x,oldx);
    while ((max_x-min_x)>dx) {
      Real x=(max_x+min_x)/2.;
      if (mu_of_x(T,x)>mu) {
        max_x=x;
      }
      else {
        min_x=x;
      }
    }
    return (max_x+min_x)/2.;
  }

HTE_ThermoData::HTE_ThermoData(void): PhaseThermoData() {cur_phase=this;}

HTE_ThermoData::HTE_ThermoData(const Structure &lattice,
			       const SpaceGroup &space_group,
			       const Structure &str,
			       const LinkedList<ArrayrVector3d> &cluster_list,
			       const LinkedList<Real> &eci_list): PhaseThermoData() {
    init(lattice,space_group,str,cluster_list,eci_list);
}

void HTE_ThermoData::init(const Structure &lattice,
		     const SpaceGroup &space_group,
		     const Structure &str,
		     const LinkedList<ArrayrVector3d> &cluster_list,
		     const LinkedList<Real> &eci_list) {
  simple_init(lattice.atom_pos.get_size(),space_group,cluster_list,eci_list);
}

void HTE_ThermoData::simple_init(int atom_per_unit_cell,
			    const SpaceGroup &space_group,
			    const LinkedList<ArrayrVector3d> &cluster_list,
			    const LinkedList<Real> &eci_list) {
    cur_phase=this;
    int max_size=0;
    LinkedListIterator<ArrayrVector3d> c(cluster_list);
    for (; c; c++) {
      max_size=MAX(max_size,c->get_size());
    }
    LinkedListIterator<Real> e(eci_list);
    poly_x.resize(max_size+1);
    zero_array(&poly_x);
    c.init(cluster_list);
    for (; c; c++,e++) {
      poly_x(c->get_size())+=(*e)*(Real)calc_multiplicity(*c,space_group.cell,space_group.point_op,space_group.trans)/(Real)atom_per_unit_cell;
    }
  }

  void HTE_ThermoData::set_T_mu(Real T, Real mu) {
    cur_x=x_of_mu(T,mu);
    cur_phi=F_of_x(T,cur_x)-mu*cur_x;
  }

  void HTE_ThermoData::read(istream &file) {
    file >> poly_x;
  }

  void HTE_ThermoData::write(ostream &file) {
    file << object_label() << endl;
    file << poly_x;
  }

Real HTE_ThermoData::dx=1e-3;

  SampledThermoData::SampledThermoData(void): PhaseThermoData() {}

  SampledThermoData::SampledThermoData(Real _T0, Real _dT, const Array<Real> &_mu0,
                     Real _dmu, const Array<ArrayReal> &_xs,
                     const Array<ArrayReal> &_phis):
    PhaseThermoData(), T0(_T0), dT(_dT), mu0(_mu0), dmu(_dmu), xs(_xs), phis(_phis) {
    cur_phase=this;
  }

  void SampledThermoData::set_T_mu(Real T, Real mu) {
    cur_phi=MAXFLOAT;
    cur_x=0.;
    int iT=(int)floor((T-T0)/dT);
    if (iT<0 || iT>=mu0.get_size()-1) return;
    Real fT=(T-T0-(Real)iT*dT)/dT;
    int imu[2];
    for (int l=0; l<2; l++) {
      imu[l]=(int)floor((mu-mu0(iT+l))/dmu);
      if (imu[l]<0) {
        cur_x=xs(iT+l)(0);
        return;
      }
      if (imu[l]>=xs(iT+l).get_size()-1) {
        cur_x=xs(iT+l)(xs(iT+l).get_size()-1);
        return;
      }
    }
    Real fmu=(mu-(mu0(iT)+imu[0]*dmu))/dmu;
    cur_x  =(1-fT)*((1-fmu)*xs(iT  )(imu[0])+fmu*xs(iT  )(imu[0]+1))
            +  fT *((1-fmu)*xs(iT+1)(imu[1])+fmu*xs(iT+1)(imu[1]+1));
    cur_phi=(1-fT)*((1-fmu)*phis(iT  )(imu[0])+fmu*phis(iT  )(imu[0]+1))
            +  fT *((1-fmu)*phis(iT+1)(imu[1])+fmu*phis(iT+1)(imu[1]+1));
  }

  void SampledThermoData::read(istream &file) {
    file >> T0;
    file >> dT;
    file >> mu0;
    file >> dmu;
    file >> xs;
    file >> phis;
  }

  void SampledThermoData::write(ostream &file) {
    file << object_label() << endl;
    file << T0 << endl;
    file << dT << endl;
    file << mu0 << endl;
    file << dmu << endl;
    file << xs << endl;
    file << phis << endl;
  }

  PhaseThermoData* PhaseListT::find_phase(Real T) {
    LinkedListIterator<PhaseThermoData> i(list);
    LinkedListIterator<Real> j(T_limit);
    while (j && (*j)<T) {i++; j++;}
    return i;
  }

  PhaseListT::PhaseListT(void): PhaseThermoData(), list() {cur_phase=this;}

  PhaseThermoData *PhaseListT::add_phase(PhaseThermoData *phase, Real Tmax) {
    list << phase;
    T_limit << new Real(Tmax);
    return phase;
  }

  void PhaseListT::set_T_mu(Real T, Real mu) {
    PhaseThermoData *phase=find_phase(T);
    if (!phase) {
      cur_phi=MAXFLOAT;
      cur_x=0;
      return;
    }
    phase->set_T_mu(T,mu);
    cur_phi=phase->phi();
    cur_x=phase->x();
  }

  void PhaseListT::read(istream &file) {
    int size;
    file >> size;
    for (int i=0; i<size; i++) {
      Real T;
      file >> T;
      T_limit << new Real(T);
    }
    for (int i=0; i<size; i++) {
      PhaseThermoData *phase=read_phase(file);
      list << phase;
    }
  }

  void PhaseListT::write(ostream &file) {
    file << object_label() << endl;
    file << list.get_size() << endl;
    LinkedListIterator<Real> iT(T_limit);
    for ( ; iT; iT++) {
      file << *iT << endl;
    }
    LinkedListIterator<PhaseThermoData> il(list);
    for ( ; il; il++) {
      il->write(file);
    }
  }

  PhaseThermoData* PhaseList::find_phase(Real T, Real mu) {
    LinkedListIterator<PhaseThermoData> i(list);
    Real min_phi=MAXFLOAT;
    PhaseThermoData *phase=NULL;
    for (; i; i++) {
      if (i->phi(T,mu)<min_phi) {
        min_phi=i->phi();
        phase=i;
      }
    }
    return phase;
  }

  PhaseList::PhaseList(void): PhaseThermoData(), list() {cur_phase=this;}

  PhaseThermoData *PhaseList::add_phase(PhaseThermoData *phase) {
    list << phase;
    return phase;
  }

  void PhaseList::set_T_mu(Real T, Real mu) {
    PhaseThermoData *phase=find_phase(T,mu);
    phase->set_T_mu(T,mu);
    cur_phi=phase->phi();
    cur_x=phase->x();
  }

  void PhaseList::read(istream &file) {
    int size;
    file >> size;
    for (int i=0; i<size; i++) {
      PhaseThermoData *phase=read_phase(file);
      list << phase;
    }
  }

  void PhaseList::write(ostream &file) {
    file << object_label() << endl;
    file << list.get_size() << endl;
    LinkedListIterator<PhaseThermoData> il(list);
    for ( ; il; il++) {
      il->write(file);
    }
  }

  MiscibilityGap::MiscibilityGap(void): PhaseThermoData() {hi=NULL; lo[0]=NULL; lo[1]=NULL;}

  MiscibilityGap::MiscibilityGap(PhaseThermoData *_hi, PhaseThermoData *_lo0, PhaseThermoData *_lo1): PhaseThermoData() {
    hi=_hi;
    lo[0]=_lo0;
    lo[1]=_lo1;
  }

  MiscibilityGap::~MiscibilityGap() {delete hi; delete lo[0]; delete lo[1];}

  void MiscibilityGap::set_T_mu(Real T, Real mu) {
    if (lo[0]->phi(T,mu) < lo[1]->phi(T,mu)) {
      if (hi->phi(T,mu) < lo[0]->phi(T,mu))
	cur_phase=hi;
      else
        cur_phase=lo[0];
    }
    else {
      if (hi->phi(T,mu) < lo[1]->phi(T,mu))
	cur_phase=hi;
      else
        cur_phase=lo[1];
    }
    cur_phi=cur_phase->phi();
    cur_x=cur_phase->x();
  }

  void MiscibilityGap::read(istream &file) {
    hi=read_phase(file);
    lo[0]=read_phase(file);
    lo[1]=read_phase(file);
  }

  void MiscibilityGap::write(ostream &file) {
    file << object_label() << endl;
    hi->write(file);
    lo[0]->write(file);
    lo[1]->write(file);
  }

  PhaseObjectList::PhaseObjectList(void): list() {
    list << new MF_ThermoData() << new LTE_ThermoData() 
         << new HTE_ThermoData() << new SampledThermoData()
         << new PhaseListT()  << new PhaseList()
         << new MiscibilityGap();
  }

PhaseObjectList phaseobjectlist;

PhaseThermoData *read_phase(istream &file) {
  char buf[MAX_LINE_LEN];
  buf[0]=0;
  do {
    file.get(buf,MAX_LINE_LEN-1);
    file.get();
  } while (strlen(buf)==0);
  if (file.eof()) return NULL;
  LinkedListIterator<PhaseThermoData> i(phaseobjectlist.list);
  for ( ;i; i++) {
    if (strcmp(buf,i->object_label())==0) {
      PhaseThermoData *phase=i->new_object();
      phase->read(file);
      return phase;
    }
  }
  return NULL;
}

Real tiny_phi=1e-4;

void recurse_phase_boundaries(LinkedList<PhaseBoundary> *pb_list, const LinkedList<PhaseThermoData> &phases, Real T, Real mu0, Real mu1, Real dx) {
  LinkedListIterator<PhaseThermoData> a(phases);
  PhaseThermoData *phase[2]={NULL,NULL};
  Real phi[2]={MAXFLOAT,MAXFLOAT};
  for ( ; a; a++) {
    if (a->phi(T,mu0)<phi[0]-tiny_phi) {
      phase[0]=a;
      phi[0]=a->phi(T,mu0);
    }
    if (a->phi(T,mu1)<phi[1]-tiny_phi) {
      phase[1]=a;
      phi[1]=a->phi(T,mu1);
    }
  }
  if (phase[0]!=phase[1] ||
      (phase[0]==phase[1] && phase[0]->get_sub_phase(T,mu0)!=phase[1]->get_sub_phase(T,mu1))
     ) {
    Real mu;
    Real phi0=phase[0]->phi(T,mu0);
    Real phi1=phase[1]->phi(T,mu1);
    Real x0=phase[0]->x(T,mu0);
    Real x1=phase[1]->x(T,mu1);
    //    cout << phi0 << " " << x0 << " " << phi1 << " " << x1 << endl;
    mu=((phi1+x1*mu1)-(phi0+x0*mu0))/(x1-x0);
    //    cout << mu0 << " " << mu << " " << mu1 << endl;
    if (mu<mu0+0.1*(mu1-mu0) || mu>mu1-0.1*(mu1-mu0)) {
      mu=(mu0+mu1)/2.;
    }
    if (fabs(phase[0]->x(T,mu0) - phase[0]->x(T,mu0-(mu-mu0))) < dx && fabs(phase[1]->x(T,mu1) - phase[1]->x(T,mu1+(mu1-mu))) < dx ) {
      PhaseBoundary *pb=new PhaseBoundary;
        pb->x[0]=phase[0]->x(T,mu0);
        pb->x[1]=phase[1]->x(T,mu1);
        pb->phase[0]=phase[0];
        pb->phase[1]=phase[1];
        (*pb_list) << pb;
    }
    else {
        recurse_phase_boundaries(pb_list, phases,T,mu0,mu,dx);
        recurse_phase_boundaries(pb_list, phases,T,mu,mu1,dx);
    }
  }
}

void calc_phase_boundaries(LinkedList<PhaseBoundary> *pb_list, const LinkedList<PhaseThermoData> &phases,
                           Real T, Real dx, const Array<Real> &mu_hints) {
  for (int i=1; i<mu_hints.get_size(); i++) {
    recurse_phase_boundaries(pb_list, phases,T,mu_hints[i-1],mu_hints[i],dx);
  }
}

int same_phases(const LinkedList<PhaseBoundary> &pb_list1, const LinkedList<PhaseBoundary> &pb_list2) {
  LinkedListIterator<PhaseBoundary> i1(pb_list1);
  LinkedListIterator<PhaseBoundary> i2(pb_list2);
  for (; i1 && i2; i1++, i2++) {
    if (i1->phase[0]!=i2->phase[0] || i1->phase[1]!=i2->phase[1]) break;
  }
  if (i1 || i2) {
    return 0;
  }
  else {
    return 1;
  }
}

void add_one_T(LinkedList<TwoPhaseRegion> *two_phase_list, Real T, const LinkedList<PhaseBoundary> &pb_list) {
  LinkedListIterator<PhaseBoundary> i(pb_list);
  for (; i; i++) {
    LinkedListIterator<TwoPhaseRegion> j(*two_phase_list);
    for (; j; j++) {
      if (i->phase[0]==j->phase[0] && i->phase[1]==j->phase[1]) {
        break;
      }
    }
    if (!j) {
      TwoPhaseRegion *tpr=new TwoPhaseRegion;
      tpr->phase[0]=i->phase[0];
      tpr->phase[1]=i->phase[1];
      two_phase_list->add(tpr,j);
    }
    Array<Real> *conc=new Array<Real>(2);
    (*conc)(0)=i->x[0];
    (*conc)(1)=i->x[1];
    j->x_list << conc;
    j->T_list << new Real(T);
  }
}

void calc_phase_boundaries(LinkedList<TwoPhaseRegion> *two_phase_list, const LinkedList<PhaseThermoData> &phases,
                           Real dT, Real dx, const Array<Real> &mu_hints) {
  LinkedList<PhaseBoundary> old_pb_list;
  LinkedList<PhaseBoundary> pb_list;
  calc_phase_boundaries(&old_pb_list, phases,dT,dx,mu_hints);
  Real T=2*dT;
  while (1) {
    calc_phase_boundaries(&pb_list, phases,T,dx,mu_hints);
    if (!same_phases(pb_list,old_pb_list)) {
      LinkedListIterator<PhaseBoundary> i(old_pb_list);
      LinkedListIterator<PhaseBoundary> j(old_pb_list);
      j++;
      while (j) {
        LinkedListIterator<PhaseBoundary> k(pb_list);
        for ( ; k; k++) {
          if (k->phase[0]==i->phase[1] || k->phase[1]==i->phase[1]) break;
        }
        if (!k) {
          Real mean=(i->x[1] + j->x[0])/2.;
          i->x[1]=mean;
          j->x[0]=mean;
        }
        i++;
        j++;
      }
      add_one_T(two_phase_list, T,old_pb_list);
    }
    add_one_T(two_phase_list, T,pb_list);
    if (pb_list.get_size()==0) break;
    T+=dT;
    old_pb_list.delete_all();
    transfer_list(&old_pb_list,&pb_list);
  }
}

void write_phase_diagram(const LinkedList<TwoPhaseRegion> &two_phase_list, ostream &file) {
  LinkedListIterator<TwoPhaseRegion> i(two_phase_list);
  for (; i; i++) {
    Array<ArrayReal> x;
    Array<Real> T;
    LinkedList_to_Array(&x,i->x_list);
    LinkedList_to_Array(&T,i->T_list);
    for (int j=0; j<T.get_size(); j++) {
      file << x(j)(0) << " " << T(j) << endl;
    }
    for (int j=T.get_size()-1; j>=0; j--) {
      file << x(j)(1) << " " << T(j) << endl;
    }
    file << endl;
  }
}

Real E_from_phi(Real T, Real mu, PhaseThermoData *p) {
  Real dT=MIN(T,1e-3);
  return p->phi(T,mu)-T*(p->phi(T+dT,mu)-p->phi(T-dT,mu))/(2*dT);
}
