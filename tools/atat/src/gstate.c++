#include "gstate.h"

void calc_convex_hull(Array<int> *p_gs_str,
                      const Array<Real> &concentration, const Array<Real> &energy) {
  Array<int> gs_str(concentration.get_size());
  int nb_in_hull=0;

  Real cur_concentration=MAXFLOAT;
  Real cur_energy=MAXFLOAT;
  int cur_gs;
  for (int i=0; i<energy.get_size(); i++) {
    if ((concentration(i)<cur_concentration) || (concentration(i)==cur_concentration && energy(i)<cur_energy)) {
      cur_concentration=concentration(i);
      cur_energy=energy(i);
      cur_gs=i;
    }
  }
  gs_str(nb_in_hull)=cur_gs;
  nb_in_hull++;

  while (1) {
    Real best_slope=MAXFLOAT;
    for (int i=0; i<energy.get_size(); i++) {
      if (concentration(i)>cur_concentration) {
        Real slope=(energy(i)-cur_energy)/(concentration(i)-cur_concentration);
        if (slope<best_slope) {
          best_slope=slope;
          cur_gs=i;
        }
      }
    }
    if (best_slope==MAXFLOAT) break;
    cur_concentration=concentration(cur_gs);
    cur_energy=energy(cur_gs);
    gs_str(nb_in_hull)=cur_gs;
    nb_in_hull++;
  }
  p_gs_str->resize(nb_in_hull);
  for (int i=0; i<nb_in_hull; i++) {
    (*p_gs_str)(i)=gs_str(i);
  }
}

void calc_convex_hull(Array<Real> *p_gs_concentration, Array<Real> *p_gs_energy,
			const Array<Real> &concentration, const Array<Real> &energy) {
  Array<int> gs_str;
  
  calc_convex_hull(&gs_str, concentration,energy);
  p_gs_concentration->resize(gs_str.get_size());
  p_gs_energy->resize(gs_str.get_size());
  for (int i=0; i<gs_str.get_size(); i++) {
    (*p_gs_concentration)(i)=concentration(gs_str(i));
    (*p_gs_energy)(i)=energy(gs_str(i));
  }
}

void convex_hull_from_gs(Array<Real> *p_gs_concentration, Array<Real> *p_gs_energy, const Array<int> &gs_str,
			const Array<Real> &concentration, const Array<Real> &energy) {
  p_gs_concentration->resize(gs_str.get_size());
  p_gs_energy->resize(gs_str.get_size());
  for (int i=0; i<gs_str.get_size(); i++) {
    (*p_gs_concentration)(i)=concentration(gs_str(i));
    (*p_gs_energy)(i)=energy(gs_str(i));
  }
}

int ground_state_problem(Array<int> *p_problem, const Array<int> &ground_state,
    const Array<Real> &concentration, const Array<Real> &fitted_energy,
    Real minc, Real maxc) {
  Real minminc=min(concentration);
  Real maxmaxc=max(concentration);
  if (minc<minminc) minc=minminc;
  if (maxc>maxmaxc) maxc=maxmaxc;
  int a_problem=0;
  p_problem->resize(fitted_energy.get_size());
  zero_array(p_problem);

  if (near_zero(maxc-minc)) {
    int index_min=-1;
    Real min_energy=MAXFLOAT;
    for (int str=0; str<concentration.get_size(); str++) {
      if (near_zero(concentration(str)-minc) && fitted_energy(str)<min_energy) {
	min_energy=fitted_energy(str);
	index_min=str;
      }
    }
    if (index_min!=-1) {
      int gs;
      for (gs=0; gs<ground_state.get_size(); gs++) {
	if (concentration(ground_state(gs)) > minc-zero_tolerance) break;
      }
      if (!near_zero(concentration(ground_state(gs))-minc)) {
	(*p_problem)(ground_state(gs))=1;
	(*p_problem)(ground_state(gs-1))=1;
	(*p_problem)(index_min)=1;
	a_problem=1;
      }
      else if (index_min!=ground_state(gs)) {
	(*p_problem)(ground_state(gs))=1;
	(*p_problem)(index_min)=1;
	a_problem=1;
      }
    }
  }
  else {
    int min_gs;
    for (min_gs=1; min_gs<ground_state.get_size()-1; min_gs++) {
      if (concentration(ground_state(min_gs))-zero_tolerance > minc) break;
    }
    int max_gs;
    for (max_gs=ground_state.get_size()-1; max_gs>1; max_gs--) {
      if (concentration(ground_state(max_gs-1))+zero_tolerance < maxc) break;
    }
    
    for (int gs=min_gs; gs<=max_gs; gs++) {
      Real c0=concentration(ground_state(gs-1));
      Real c1=concentration(ground_state(gs));
      Real e0=fitted_energy(ground_state(gs-1));
      Real e1=fitted_energy(ground_state(gs));
      Real slope=(e1-e0)/(c1-c0);
      if (gs==min_gs) {
	e0=e0+slope*(minminc-c0);
	c0=minminc;
      }
      if (gs==max_gs) {
	e1=e1+slope*(maxmaxc-c1);
	c1=maxmaxc;
      }
      if (gs>=min_gs+1) {
	Real c_1=concentration(ground_state(gs-2));
	Real e_1=fitted_energy(ground_state(gs-2));
	if ( (e1-e0)/(c1-c0) < (e0-e_1)/(c0-c_1) ) {
	  (*p_problem)(ground_state(gs-2))=1;
	  (*p_problem)(ground_state(gs-1))=1;
	  (*p_problem)(ground_state(gs))=1;
	  a_problem=1;
	}
      }
      Real a=(c1*e0-c0*e1)/(c1-c0);
      Real b=(e1-e0)/(c1-c0);
      for (int str=0; str<fitted_energy.get_size(); str++) {
	Real c=concentration(str);
	if (c>=c0 && c<=c1) {
	  Real hull_e=a+b*c;
	  if (str!=ground_state(gs-1) && str!=ground_state(gs) && fitted_energy(str)<=hull_e) {
	    (*p_problem)(str)=1;
	    (*p_problem)(ground_state(gs-1))=1;
	    (*p_problem)(ground_state(gs  ))=1;
	    a_problem=1;
	  }
	}
      }
    }
  }
  return a_problem;
}

int ground_state_problem(Array<int> *p_problem, const Array<Real> &gs_concentration, const Array<Real> &gs_energy,
       const Array<Real> &concentration, const Array<Real> &energy,
    Real minc, Real maxc) {
  if (concentration.get_size()==0) {
    p_problem->resize(0);
    return 0;
  }
  Real minminc=MIN(min(concentration),min(gs_concentration));
  Real maxmaxc=MAX(max(concentration),max(gs_concentration));
  if (minc<minminc) minc=minminc;
  if (maxc>maxmaxc) maxc=maxmaxc;
  int a_problem=0;
  p_problem->resize(energy.get_size());
  zero_array(p_problem);
  if (near_zero(maxc-minc)) {
    minc-=2*zero_tolerance;
    maxc+=2*zero_tolerance;
  }
  int min_gs;
  for (min_gs=1; min_gs<gs_concentration.get_size()-1; min_gs++) {
    if (gs_concentration(min_gs)-zero_tolerance > minc) break;
  }
  int max_gs;
  for (max_gs=gs_concentration.get_size()-1; max_gs>1; max_gs--) {
    if (gs_concentration(max_gs-1)+zero_tolerance < maxc) break;
  }
  for (int gs=min_gs; gs<=max_gs; gs++) {
    Real c0=gs_concentration(gs-1);
    Real c1=gs_concentration(gs);
    Real e0=gs_energy(gs-1);
    Real e1=gs_energy(gs);
    Real slope=(e1-e0)/(c1-c0);
    if (gs==min_gs) {
      e0=e0+slope*(minminc-c0);
      c0=minminc;
    }
    if (gs==max_gs) {
      e1=e1+slope*(maxmaxc-c1);
      c1=maxmaxc;
    }
    Real a=(c1*e0-c0*e1)/(c1-c0);
    Real b=(e1-e0)/(c1-c0);
    for (int str=0; str<energy.get_size(); str++) {
      Real c=concentration(str);
      Real hull_e=a+b*c;
      if (c>c0 && c<c1) {
	if (energy(str)<=hull_e) {
	  (*p_problem)(str)=1;
	  a_problem=1;
	}
      }
      else if (c==c0 || c==c1) {
	if (energy(str)<hull_e) {
	  (*p_problem)(str)=1;
	  a_problem=1;
	}
      }
    }
  }
  return a_problem;
}

void calc_formation(Array<Real> *p_form_value, const Array<Real> &concentration, const Array<Real> &value, Real *e_c0, Real *e_c1) {
  Real pure_A=MAXFLOAT;
  Real pure_B=MAXFLOAT;
  if (e_c0 && e_c1 && (*e_c0)!=MAXFLOAT && (*e_c1)!=MAXFLOAT) {
    pure_A=*e_c0;
    pure_B=*e_c1;
  }
  else {
    Real minc=min(concentration);
    Real maxc=max(concentration);
    Real mince,maxce;
    for (int i=0; i<value.get_size(); i++) {
      if (concentration(i)==minc) mince=value(i);
      if (concentration(i)==maxc) maxce=value(i);
    }
    if (near_zero(maxc-minc)) {
      pure_A=mince;
      pure_B=maxce;
    }
    else {
      pure_A=(maxce-mince)/(maxc-minc)*(-1.-minc)+mince;
      pure_B=(maxce-mince)/(maxc-minc)*( 1.-minc)+mince;
    }
    if (e_c0 && e_c1) {
      (*e_c0)=pure_A;
      (*e_c1)=pure_B;
    }
  }
  p_form_value->resize(value.get_size());
  Real a=(pure_B+pure_A)/2.;
  Real b=(pure_B-pure_A)/2.;
  for (int i=0; i<value.get_size(); i++) {
    (*p_form_value)(i)=value(i)-(a+b*concentration(i));
  }
}

void calc_hull_energy(Array<Real> *hull_energy, const Array<int> &ground_state, const Array<Real> &concentration, const Array<Real> &energy) {
  hull_energy->resize(energy.get_size());
  for (int gs=1; gs<ground_state.get_size(); gs++) {
    Real c0=concentration(ground_state(gs-1));
    Real c1=concentration(ground_state(gs));
    Real e0=energy(ground_state(gs-1));
    Real e1=energy(ground_state(gs));
    Real a=(c1*e0-c0*e1)/(c1-c0);
    Real b=(e1-e0)/(c1-c0);
    for (int str=0; str<energy.get_size(); str++) {
      Real c=concentration(str);
      if (c>=c0 && c<=c1) {
	(*hull_energy)(str)=a+b*c;
      }
    }
  }
}
