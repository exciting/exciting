#include <fstream.h>
#include <strstream.h>
#include "lstsqr.h"

#define MAXLINELEN 512

Real ground_state_badness(const Array<int> &ground_state, const Array<Real> &concentration, const Array<Real> &fitted_energy) {
  Real badness=0;
  for (int gs=1; gs<ground_state.get_size()-1; gs++) {
    Real c0=concentration(ground_state(gs-1));
    Real c1=concentration(ground_state(gs+1));
    Real e0=fitted_energy(ground_state(gs-1));
    Real e1=fitted_energy(ground_state(gs+1));
    Real a=(c1*e0-c0*e1)/(c1-c0);
    Real b=(e1-e0)/(c1-c0);
    int gs_str=ground_state(gs);
    Real de=a+b*concentration(gs_str)-fitted_energy(gs_str);
    if (de<0) badness+=-de;
  }
  for (int gs=1; gs<ground_state.get_size(); gs++) {
    Real c0=concentration(ground_state(gs-1));
    Real c1=concentration(ground_state(gs));
    Real e0=fitted_energy(ground_state(gs-1));
    Real e1=fitted_energy(ground_state(gs));
    Real a=(c1*e0-c0*e1)/(c1-c0);
    Real b=(e1-e0)/(c1-c0);
    for (int str=0; str<fitted_energy.get_size(); str++) {
      Real c=concentration(str);
      if (c>=c0 && c<=c1) {
	Real hull_e=a+b*c;
	if (fitted_energy(str)<hull_e) badness+=(hull_e-fitted_energy(str));
      }
    }
  }
  return badness;
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

Real ground_state_problem(Array<Real> *p_problem, const Array<int> &ground_state, const Array<Real> &concentration, const Array<Real> &fitted_energy) {
  p_problem->resize(fitted_energy.get_size());
  zero_array(p_problem);
  for (int gs=1; gs<ground_state.get_size(); gs++) {
    Real c0=concentration(ground_state(gs-1));
    Real c1=concentration(ground_state(gs));
    Real e0=fitted_energy(ground_state(gs-1));
    Real e1=fitted_energy(ground_state(gs));
    Real a=(c1*e0-c0*e1)/(c1-c0);
    Real b=(e1-e0)/(c1-c0);
    for (int str=0; str<fitted_energy.get_size(); str++) {
      Real c=concentration(str);
      if (c>=c0 && c<=c1) {
	Real hull_e=a+b*c;
	if (fitted_energy(str)<hull_e) {
	  (*p_problem)(str)=1.;
	  (*p_problem)(ground_state(gs-1))=1.;
	  (*p_problem)(ground_state(gs  ))=1.;
	}
      }
    }
  }
}

Real ground_state_problem(Array<Real> *p_problem, const Array<Real> &gs_concentration, const Array<Real> &gs_energy, const Array<Real> &concentration, const Array<Real> &energy) {
  p_problem->resize(energy.get_size());
  zero_array(p_problem);
  for (int gs=1; gs<gs_concentration.get_size(); gs++) {
    Real c0=gs_concentration(gs-1);
    Real c1=gs_concentration(gs);
    Real e0=gs_energy(gs-1);
    Real e1=gs_energy(gs);
    Real a=(c1*e0-c0*e1)/(c1-c0);
    Real b=(e1-e0)/(c1-c0);
    for (int str=0; str<energy.get_size(); str++) {
      Real c=concentration(str);
      if (c>=c0 && c<=c1) {
	Real hull_e=a+b*c;
	if (energy(str)<hull_e) {
	  (*p_problem)(str)=hull_e-energy(str);
	}
      }
    }
  }
}

Real gs_penalty=1000.;
Real gs_badness_mult=100.;

Real add_best_eci(Array<int> *selected_eci, int nb_eci,
                  const Array2d<Real> &all_rho, const Array<Real> &energy,
		  const Array<Real> &weight,
                  const Array<Real> &concentration,
		  const Array<int> &ground_state) {
  Array<Real> included(all_rho.get_size()(1));
  zero_array(&included);
  for (int eci=0; eci<nb_eci; eci++) {
    included((*selected_eci)(eci))=1;
  }
  Array2d<Real> rho;
  (*selected_eci)(nb_eci)=-1;
  extract_columns(&rho, all_rho, *selected_eci, nb_eci+1);
  int best_eci=-1;
  Real best_cv=MAXFLOAT;
  for (int new_eci=0; new_eci<all_rho.get_size()(1); new_eci++) {
    if (!included(new_eci)) {
      for (int s=0; s<all_rho.get_size()(0); s++) {
        rho(s,nb_eci)=all_rho(s,new_eci);
      }
      Real cv=calc_cv(rho,energy,weight);
      Array<Real> fitted_energy;
      predict_gls(&fitted_energy, rho,energy,weight);
      Real gs_badness=ground_state_badness(ground_state,concentration,fitted_energy);
      cv+=gs_badness_mult*gs_badness;
      if (gs_badness>0.) cv+=gs_penalty;
      if (cv<=best_cv) {
        best_eci=new_eci;
        best_cv=cv;
      }
    }
  }
  (*selected_eci)(nb_eci)=best_eci;
  return best_cv;
}

Real remove_worse_eci(Array<int> *selected_eci, int nb_eci,
                      const Array2d<Real> &all_rho, const Array<Real> &energy, const Array<Real> &weight,
                      const Array<Real> &concentration, const Array<int> &ground_state, int min_nb_eci=0) {
  int worse_eci=-1;
  Real best_cv=MAXFLOAT;
  Array<int> cur_eci(selected_eci->get_size());
  for (int bad_eci=min_nb_eci; bad_eci<nb_eci; bad_eci++) {
    for (int i=0; i<bad_eci; i++) {cur_eci(i)=(*selected_eci)(i);}
    for (int i=bad_eci+1; i<nb_eci; i++) {cur_eci(i-1)=(*selected_eci)(i);}
    Array2d<Real> rho;
    extract_columns(&rho, all_rho, cur_eci, nb_eci-1);
    Real cv=calc_cv(rho,energy,weight);
    Array<Real> fitted_energy;
    predict_gls(&fitted_energy, rho,energy,weight);
    Real gs_badness=ground_state_badness(ground_state,concentration,fitted_energy);
    cv+=gs_badness_mult*gs_badness;
    if (gs_badness>0.) cv+=gs_penalty;
    if (cv<=best_cv) {
      worse_eci=bad_eci;
      best_cv=cv;
    }
  }
  for (int i=worse_eci+1; i<nb_eci; i++) {(*selected_eci)(i-1)=(*selected_eci)(i);}
  (*selected_eci)(nb_eci-1)=0;
  return best_cv;
}

void calc_formation(Array<Real> *p_form_value, const Array<Real> &concentration, const Array<Real> &value, Real *e_c0=NULL, Real *e_c1=NULL) {
  Real pure_A=MAXFLOAT;
  Real pure_B=MAXFLOAT;
  if (e_c0 && e_c1 && (*e_c0)!=MAXFLOAT && (*e_c1)!=MAXFLOAT) {
    pure_A=*e_c0;
    pure_B=*e_c1;
  }
  else {
    for (int i=0; i<value.get_size(); i++) {
      if (concentration(i)==-1) pure_A=value(i);
      if (concentration(i)== 1) pure_B=value(i);
    }
    if (pure_A==MAXFLOAT || pure_B==MAXFLOAT)
      ERRORQUIT("No pure A or Pure B structure found in calc_formation.");
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

void calc_convex_hull(Array<int> *p_gs_str,
                      const Array<Real> &concentration, const Array<Real> &energy) {
  Array<int> gs_str(concentration.get_size());
  int nb_in_hull=0;

  Real cur_concentration=MAXFLOAT;
  Real cur_energy=MAXFLOAT;
  int cur_gs;
  for (int i=0; i<energy.get_size(); i++) {
    if (concentration(i)<=cur_concentration && energy(i)<=cur_energy) {
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

extern char *helpstring;

void main(int argc, char *argv[]) {
  if (argc>1) {
    cout << helpstring;
    exit(1);
  }
  int nb_str,nb_eci;
  system("grep -v '^[ 	]*$' corr.in | wc -l > size.out");
  system("grep -v '^[ 	]*$' corr.in | head -1 | wc -w >> size.out");
  {
    ifstream file("size.out");
    file >> nb_str >> nb_eci;
  }
  int add_empty=0;
  Array<Real> point_mult;
  Real gs_prec=0.;
  {
    ifstream file("info.in");
    if (!file) ERRORQUIT("Unable to open file info.in");
    file >> gs_penalty;
    file >> gs_badness_mult;
    file >> add_empty;
    file >> gs_prec;
    add_empty=(add_empty!=0);
    nb_eci+=add_empty;
    point_mult.resize(nb_eci);
    zero_array(&point_mult);
    while (!file.eof()) {
      int which_col;
      file >> which_col;
      if (file.eof()) break;
      file >> point_mult(which_col-1+add_empty);
    }      
  }

  Array<int> known(nb_str);
  zero_array(&known);
  int nb_known=0;
  {
    ifstream file("energy.in");
    if (!file) ERRORQUIT("Unable to open file energy.in");
    for (int str=0; str<nb_str; str++) {
      char s[MAXLINELEN];
      file >> s;
      istrstream num(s);
      Real tmp=MAXFLOAT;
      num >> tmp;
      if (tmp!=MAXFLOAT) {
        known(str)=1;
	nb_known++;
      }
      else {
        known(str)=0;
      }
    }
  }
  Array<Real> energy(nb_known);
  {
    ifstream file("energy.in");
    if (!file) ERRORQUIT("Unable to open file energy.in");
    int known_i=0;
    for (int str=0; str<nb_str; str++) {
      if (known(str)) {
        file >> energy(known_i);
	known_i++;
      }
      else {
	char tmp;
	file >> tmp;
      }
    }
  }
  {
    Real maxenergy=0.;
    for (int i=0; i<energy.get_size(); i++) {
      if (fabs(energy(i))>maxenergy) maxenergy=fabs(energy(i));
    }
    gs_penalty=(gs_penalty!=0)*maxenergy*10;
  }

  Array2d<Real> all_rho(nb_known,nb_eci);
  Array2d<Real> all_rho_unknown(nb_str-nb_known,nb_eci);
  {
    ifstream file("corr.in");
    if (!file) ERRORQUIT("Unable to open file corr.in");
    int known_i=0;
    int unknown_i=0;
    for (int str=0; str<nb_str; str++) {
      if (known(str)) {
        if (add_empty) all_rho(known_i,0)=1.;
        for (int eci=add_empty; eci<nb_eci; eci++) {
          file >> all_rho(known_i,eci);
        }
        known_i++;
      }
      else {
        if (add_empty) all_rho_unknown(unknown_i,0)=1.;
        for (int eci=add_empty; eci<nb_eci; eci++) {
          file >> all_rho_unknown(unknown_i,eci);
        }
        unknown_i++;
      }
    }
  }

  Array<Real> concentration;
  product(&concentration, all_rho,point_mult);
  Array<Real> concentration_unknown;
  product(&concentration_unknown, all_rho_unknown,point_mult);
    
  int max_nb_eci=MIN(nb_eci,nb_known-1);
  Array<int> selected_eci(max_nb_eci);
  int min_nb_eci=0;
  {
    if (add_empty) {
      selected_eci(0)=0;
      min_nb_eci++;
    }
    for (int i=0; i<point_mult.get_size(); i++) {
      if (point_mult(i)>0) {
	selected_eci(min_nb_eci)=i;
	min_nb_eci++;
      }
    }
  }
  {
    ifstream file("include.in");
    if (file) {
      while (!file.eof()) {
	int tmp=-1;
	file >> tmp;
	if (tmp==-1) break;
	selected_eci(min_nb_eci)=tmp-1+add_empty;
	min_nb_eci++;
      }
    }
  }
  Array<Real> in_weight(nb_known);
  {
    ifstream file("weight.in");
    if (file) {
      int known_i=0;
      for (int str=0; str<nb_str; str++) {
        if (known(str)) {
          file >> in_weight(known_i);
          known_i++;
        }
        else {
          Real tmp;
          file >> tmp;
        }
      }
    }
    else {
      for (int str=0; str<nb_known; str++) {
        in_weight(str)=1.;
      }
    }
  }

  Array<int> ground_state;
  calc_convex_hull(&ground_state, concentration,energy);

  Array<Real> hull_energy;
  calc_hull_energy(&hull_energy, ground_state,concentration,energy);

  Array<Real> weight(nb_known);
  {
    ofstream file("weight_used.out");
    for (int str=0; str<nb_known; str++) {
      weight(str)=in_weight(str)*gs_prec/(energy(str)-hull_energy(str)+gs_prec);
      file << weight(str) << endl;
    }
  }


  int nb_eci_incl=min_nb_eci;
  Real best_cv=MAXFLOAT;
  int done;

  do {
    for (int eci=nb_eci_incl; eci<max_nb_eci; eci++) {
      Real cv=add_best_eci(&selected_eci,eci,all_rho,energy,weight,concentration,ground_state);
        if (cv<best_cv) {
	  nb_eci_incl=eci+1;
	  best_cv=cv;
cerr << cv << endl;
        }
      }
      
      Array<int> cur_selected_eci(max_nb_eci);
      cur_selected_eci.copy(selected_eci);
      done=1;
      for (int eci=nb_eci_incl; eci>min_nb_eci; eci--) {
      Real cv=remove_worse_eci(&cur_selected_eci,eci,all_rho,energy,weight,concentration,ground_state,min_nb_eci);
      if (cv<best_cv && eci>=2) {
        selected_eci.copy(cur_selected_eci);
        nb_eci_incl=eci-1;
        best_cv=cv;
cerr << cv << endl;
        done=0;
      }
    }
  } while (!done);

  Array2d<Real> rho;
  extract_columns(&rho,all_rho,selected_eci,nb_eci_incl);
  Array<Real> eci;
  calc_gls(&eci,rho,energy,weight);
  Array2d<Real> var_eci;
  calc_gls_var(&var_eci,rho,energy,weight);
    
  Array<Real> fitted_energy;
  predict_gls(&fitted_energy,rho,energy,weight);

  Array<Real> all_eci(all_rho.get_size()(1));
  zero_array(&all_eci);
  for (int i=0; i<nb_eci_incl; i++) {
    all_eci(selected_eci(i))=eci(i);
  }
  Array<Real> outlier_test(energy.get_size());
  {
    Array2d<Real> tmp;
    product(&tmp, rho,var_eci);
    Array<Real> prediction_var;
    outer_product_diag(&prediction_var, rho,tmp);
    for (int i=0; i<energy.get_size(); i++) {
      outlier_test(i)=(fitted_energy(i)-energy(i))/sqrt(prediction_var(i));
    }
  }

  Array<Real> form_energy;
  calc_formation(&form_energy, concentration,energy);
  Array<Real> form_fitted_energy;
  Real e_c0=MAXFLOAT,e_c1=MAXFLOAT;
  calc_formation(&form_fitted_energy, concentration,fitted_energy,&e_c0,&e_c1);
  Array<Real> gs_concentration,gs_energy;
  calc_convex_hull(&gs_concentration,&gs_energy, concentration,form_energy);

  Array<Real> gs_fitted_concentration,gs_fitted_energy;
  calc_convex_hull(&gs_fitted_concentration,&gs_fitted_energy, concentration,form_fitted_energy);

  Array2d<Real> rho_unknown;
  extract_columns(&rho_unknown,all_rho_unknown,selected_eci,nb_eci_incl);
  Array<Real> fitted_energy_unknown;
  product(&fitted_energy_unknown, rho_unknown,eci);
  Array<Real> form_fitted_energy_unknown;
  calc_formation(&form_fitted_energy_unknown, concentration_unknown,fitted_energy_unknown,&e_c0,&e_c1);

  Array<Real> gs_problem_unknown;
  {
    ground_state_problem(&gs_problem_unknown, gs_fitted_concentration,gs_fitted_energy, concentration_unknown, form_fitted_energy_unknown);
    Array2d<Real> tmp;
    product(&tmp, rho_unknown,var_eci);
    Array<Real> prediction_var_unknown;
    outer_product_diag(&prediction_var_unknown, rho_unknown,tmp);
    for (int i=0; i<gs_problem_unknown.get_size(); i++) {
      gs_problem_unknown(i)/=sqrt(prediction_var_unknown(i));
    }
  }
  {
    ofstream file("energy.out");
    for (int i=0; i<energy.get_size(); i++) {
      file << (1.+concentration(i))/2. << " " << form_energy(i) << endl;
    }
  }
  {
    ofstream file("eci.out");
    int j=0;
    for (int i=0; i<all_eci.get_size(); i++) {
      if (all_eci(i)!=0) {
	file << all_eci(i) << endl;
	j++;
      }
      else {
	file << 0. << endl;
      }
    }
  }
  {
    ofstream file("fitted.out");
    for (int i=0; i<fitted_energy.get_size(); i++) {
      file << (1.+concentration(i))/2. << " " << form_fitted_energy(i) << endl;
    }
  }
  {
    ofstream file("fitted_u.out");
    for (int i=0; i<fitted_energy_unknown.get_size(); i++) {
      file << (1.+concentration_unknown(i))/2. << " " << form_fitted_energy_unknown(i) << " " << gs_problem_unknown(i) << endl;
    }
  }
  {
    ofstream file("gs.out");
    for (int i=0; i<gs_energy.get_size(); i++) {
      file << (1.+gs_concentration(i))/2. << " " << gs_energy(i) << endl;
    }
  }
  {
    ofstream file("gs_fitted.out");
    for (int i=0; i<gs_fitted_energy.get_size(); i++) {
      file << (1.+gs_fitted_concentration(i))/2. << " " << gs_fitted_energy(i) << endl;
    }
  }
  Array<Real> gs_problem;
  ground_state_problem(&gs_problem,ground_state,concentration,fitted_energy);
  {
    ofstream file("weight_sug.out");
    int j=0;
    for (int i=0; i<known.get_size(); i++) {
      if (known(i)) {
	file << in_weight(j)*(1.+gs_problem(j)) << endl;
	j=j+1;
      }
      else {
	file << 0 << endl;
      }
    }
  }
  {
    ofstream file("outlier.out");
    for (int i=0; i<outlier_test.get_size(); i++) {
      file << i << " " << outlier_test(i) << endl;
    }
  }
    
  {
    ofstream file("fit.gnu");
    file << "plot 'energy.out' w p, 'fitted.out' w p, 'fitted_u.out' u 1:2 w p, 'gs.out' w l, 'gs_fitted.out' w l" << endl;
    file << "pause -1" << endl;
  }
}

