#include "calcmf.h"

CalcMeanField::CalcMeanField(const Structure &_lat, const Array<Array<int> > &_site_type_list, const Structure &_str, 
		const SpaceGroup &space_group, const LinkedList<MultiCluster> &cluster_list, 
		const Array<Array<Array<Real> > > &_corrfunc): ref_prob(), clusters(), eci(), corrfunc(_corrfunc), which_is_empty(-1), empty_mult(1.), do_LTE(0) {
  LinkedListIterator<MultiCluster> ic(cluster_list);
  for ( int i=0; ic; ic++, i++ ) {
    if (ic->clus.get_size()==0) {which_is_empty=i;}
  }
  empty_mult=(Real)_str.atom_pos.get_size()/(Real)_lat.atom_pos.get_size();
  ref_prob.resize(_str.atom_pos.get_size());
  clusters.resize(_str.atom_pos.get_size());
  rMatrix3d istrcell=!_str.cell;
  rMatrix3d ilatcell=!_lat.cell;
  for (int i=0; i<ref_prob.get_size(); i++) {
    ref_prob(i).resize(_site_type_list(_lat.atom_type(which_atom(_lat.atom_pos, _str.atom_pos(i), ilatcell))).get_size());
    if (_str.atom_type(i)==-1) {
      for (int t=0; t<ref_prob(i).get_size(); t++) {
	ref_prob(i)(t)=1./(Real)(ref_prob(i).get_size());
      }
    } else {
      zero_array(&(ref_prob(i)));
      ref_prob(i)(_str.atom_type(i))=1.;
    }
  }
  for (int i=0; i<ref_prob.get_size(); i++) {
    Array<MultiCluster> o_cluster;
    Array<int> o_which;
    find_clusters_overlapping_site(&o_cluster, &o_which, _str.atom_pos(i), space_group, cluster_list);
    clusters(i).resize(o_cluster.get_size());
    for (int c=0; c<clusters(i).get_size(); c++) {
      int nbs=o_cluster(c).clus.get_size();
      clusters(i)(c).which_eci=o_which(c);
      clusters(i)(c).site.resize(nbs);
      clusters(i)(c).site_type.resize(nbs);
      clusters(i)(c).func.resize(nbs);
      for (int s=0; s<nbs; s++) {
        clusters(i)(c).site(s)=which_atom(_str.atom_pos, o_cluster(c).clus(s), istrcell);
        clusters(i)(c).site_type(s)=o_cluster(c).site_type(s);
        clusters(i)(c).func(s)=o_cluster(c).func(s);
      }
    }  
  }
}

Real CalcMeanField::calc_point_corr(const Array<Real> &siteprob, int site_type, int func) {
  Real sigma=0.;
  for (int t=0; t<siteprob.get_size(); t++) {
    sigma+=corrfunc(site_type)(func)(t)*siteprob(t);
  }
  return sigma;
}

void CalcMeanField::set_T_mu(Real _T, const Array<Real> &_mu) {
  Array<Real> ecimu(eci);
  for (int i=0; i<_mu.get_size(); i++) {ecimu(1+i)-=_mu(i);}

  Array<Array<Real> > prob;
  prob=ref_prob;

  Real F_old=MAXFLOAT;
  Real F=0.;
  Real F_tol=1e-6;
  int step=0;
  int maxstep=100;

  while (fabs(F-F_old)>F_tol && step<maxstep) {
    F_old=F;
    Real E0=(which_is_empty==-1 ? 0. : eci(which_is_empty))*empty_mult;
    for (int i=0; i<clusters.get_size(); i++) {
      for (int c=0; c<clusters(i).get_size(); c++) {
        Real rho=1.;
        for (int s=0; s<clusters(i)(c).site.get_size(); s++) {
          rho*=calc_point_corr(prob(clusters(i)(c).site(s)), clusters(i)(c).site_type(s), clusters(i)(c).func(s));
        }
        E0+=rho*ecimu(clusters(i)(c).which_eci)/((Real)clusters(i)(c).site.get_size());
      }
    }
    Array<Array<Real> > dE(prob.get_size());
    for (int i=0; i<prob.get_size(); i++) {
      dE(i).resize(prob(i).get_size());
      zero_array(&(dE(i)));
      for (int ex=0; ex<prob(i).get_size(); ex++) {
        Array<Real> dprobi(prob(i).get_size());
        for (int t=0; t<dprobi.get_size(); t++) {
          dprobi(t)=(t==ex ? 1. : 0.)-prob(i)(t);
        }
        for (int c=0; c<clusters(i).get_size(); c++) {
          Real rho=calc_point_corr(dprobi, clusters(i)(c).site_type(0), clusters(i)(c).func(0));
          for (int s=1; s<clusters(i)(c).site.get_size(); s++) {
            rho*=calc_point_corr(prob(clusters(i)(c).site(s)), clusters(i)(c).site_type(s), clusters(i)(c).func(s));
          }
          dE(i)(ex)+=rho*ecimu(clusters(i)(c).which_eci);
        }
      }
    }
    F=E0;
    if (_T==0) break;
    for (int i=0; i<prob.get_size(); i++) {
      Real pf=0.;
      for (int ex=0; ex<dE(i).get_size(); ex++) {
        Real p=exp(-(dE(i)(ex)/_T));
        prob(i)(ex)=p;
        pf+=p;
      }
      for (int ex=0; ex<dE(i).get_size(); ex++) {
        prob(i)(ex)/=pf;
      }
      F+=-_T*log(pf);
    }
    if (do_LTE) break;
    step++;
  }
  Real disord=0.;
  for (int i=0; i<ref_prob.get_size(); i++) {
    for (int t=0; t<ref_prob(i).get_size(); t++) {
      if (ref_prob(i)(t)==1.) {disord+=1.-prob(i)(t);}
    }
  }
  lro=1.-disord/(Real)(ref_prob.get_size());
  phi=F/(Real)(ref_prob.get_size());
}

void calc_conc_from_phi(Array<Real> *pconc, GCThermoData *ptd, Real T, const Array<Real> &mu, Real dmu) {
  Array<Real> nmu(mu.get_size());
  pconc->resize(mu.get_size());
  for (int i=0; i<mu.get_size(); i++) {
    Real phi_p,phi_m;
    nmu=mu;
    nmu(i)=mu(i)+dmu;
    phi_p=ptd->set_get_phi(T,nmu);
    nmu(i)=mu(i)-dmu;
    phi_m=ptd->set_get_phi(T,nmu);
    (*pconc)(i)=-(phi_p-phi_m)/(2.*dmu);
  }
}

Real calc_E_from_phi(GCThermoData *ptd, Real T, const Array<Real> &mu, Real dT) {
  Real phi,phi_p,phi_m;
  phi=ptd->set_get_phi(T,mu);
  phi_p=ptd->set_get_phi(T+dT,mu);
  phi_m=ptd->set_get_phi(T-dT,mu);
  return (phi-T*(phi_p-phi_m)/(2.*dT));
}

Real calc_E_from_phi(CalcMeanField *ptd, const PolyInterpolatorBig<Array<Real> > &teci, Real T, const Array<Real> &mu, Real dT) {
  Real phi,phi_p,phi_m;
  Array<Real> eci;
  teci.interpol(&eci,T);
  ptd->set_eci(eci);
  phi=ptd->set_get_phi(T,mu);
  teci.interpol(&eci,T+dT);
  ptd->set_eci(eci);
  phi_p=ptd->set_get_phi(T+dT,mu);
  teci.interpol(&eci,T-dT);
  ptd->set_eci(eci);
  phi_m=ptd->set_get_phi(T-dT,mu);
  return (phi-T*(phi_p-phi_m)/(2.*dT));
}
