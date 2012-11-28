#include <fstream.h>
#include "refine.h"
#include "ridge.h"
#include "getvalue.h"

// -----------------------------;


class CE_TLambda: public ClusterExpansion {
public:
  Array<Real> cur_regul_param;
  CE_TLambda (const Structure &_parent_lattice,
                   const Array<Arrayint> &_site_type_list,
                   const Array<AutoString> &_atom_label,
                   const SpaceGroup &_spacegroup): ClusterExpansion(_parent_lattice,_site_type_list,_atom_label,_spacegroup), cur_regul_param() {};

  // Main routine calling the cluster enumeration and crossvalidation routines;
  // Initialize *pfitinfo with results (see class CEFitInfo above);
  virtual void find_best_cluster_choice(CEFitInfo *pfitinfo);
  // Main routine calling the structure enumeration routines and the variance calculation codes.;
  // Return pointer to structure found (or NULL upon failure);
  virtual StructureInfo * find_best_structure(void);
  // Finds the first few structures needed to initiate the refinement process.;
  // Return pointer to structure found (or NULL upon failure);
  // virtual StructureInfo * find_initial_structures(void);
  // Returns the least computationally intensive structure still marked as unknown.
  // Return pointer to structure (never NULL);
  // virtual StructureInfo * find_first_unknown_structure(void);

  int find_next_multiplet_choice(void);
  void calc_regularization_vector(Array<Real> *regul, const Array<Real> &param);
};

SpecificPlugIn<ClusterExpansionCreator, CustomClusterExpansionCreator<CE_TLambda> > TLambdaClusterExpansionPlugIn("tl"); // register plug-in under name "tl";

void CE_TLambda::find_best_cluster_choice(CEFitInfo *pfitinfo) {
  int max_pair_len=0;
  Array<Real> regul_param(2);
  {
    ifstream infile("tlambda.in");
    if (!infile) {
      cerr << "You need to provide a tlambda.in file containing:" << endl
	   << "  the length of the longest pair" << endl
	   << "  the 't' parameter" << endl
	   << "  the 'lambda' parameter" << endl;
      ERRORQUIT("Aborting.");
    }
    infile >> max_pair_len;
    infile >> regul_param(0) >> regul_param(1);
  }
  cur_regul_param=regul_param; // save param for later;

  reset_cluster_choice(); // select minimal cluster expansion;
  while (pclusters(1)->get_current_length()<max_pair_len) {(*pclusters(1))++;} // include requested pairs;
  Array<int> pair_clusters;
  get_cluster_choice(&pair_clusters);
  set_cluster_choice(pair_clusters);

  // setup correlation matrix, and energy and concentration vectors;
  Array2d<Real> corr_matrix;
  Array<Real> regul;
  Array<Real> energy;
  Array<Real> predictor_energy;
  Array<Real> ce_energy;
  Array<Real> concentration;
  calc_regression_matrices(&corr_matrix,&energy);
  if (corr_matrix.get_size()(0)<=corr_matrix.get_size()(1)-pair_clusters(1)) { // enough structures?;
    if (pfitinfo) {pfitinfo->status=CEFitInfo::fit_impossible;}
    return;
  }
  if (pfitinfo) {pfitinfo->status=CEFitInfo::fit_ok;}
  calc_concentration(&concentration,corr_matrix);
  init_predictors();
  // find ground states;
  Array<int> true_gs;
  calc_convex_hull(&true_gs, concentration,energy);
  // init weights (to 1);
  weight.resize(concentration.get_size());
  one_array(&weight);
  Array<int> best_choice; // to remember best choice;
  Array<Real> best_eci_mult; // to remember best eci (times multiplicity);
  Array<int> best_choice_noweight; // if weighting fails, save the best noweight solution;
  Array<Real> best_eci_mult_noweight; // 
  {
    ifstream cluster_choice_file("nbclusters.in");
    if (cluster_choice_file) {
      LinkedList<int> cluster_choice;
      cluster_choice << new int(pclusters(0)->get_current_index());
      cerr << "nbclusters.in file detected, using user-specified cluster choice" << endl;
      while (skip_delim(cluster_choice_file)) {
	int nc;
	cluster_choice_file >> nc;
	cluster_choice << new int(nc);
      }
      LinkedList_to_Array(&best_choice,cluster_choice);
    }
  }
  if (best_choice.get_size()>0) { // if user-specified cluster choice;
    calc_predictor_energy(&predictor_energy);
  }
  else { // if automatically determined cluster choice;
    while (1) { // repeat until ground states ok;
      set_cluster_choice(pair_clusters); // include requested pairs;
      int best_gs_ok=0;
      Real best_cv=MAXFLOAT;
      do { // loop through cluster choices;
	//if you want to print out trace information;
	for (int i=0; i<pclusters.get_size(); i++) {
	  cerr << pclusters(i)->get_current_index() << " ";
	}

	// calc cross-validation score;
	calc_regression_matrices(&corr_matrix,&energy);
	calc_regularization_vector(&regul,regul_param);
	calc_predictor_energy(&predictor_energy);
	diff(&ce_energy,energy,predictor_energy);
	Real cv=calc_cv_regul(corr_matrix,ce_energy,weight,1,regul);
	// do regression and add back the predictor energy;
	Array<Real> cur_eci_mult;
	Array<Real> fitted_energy;
	calc_ols_regul(&cur_eci_mult,corr_matrix,ce_energy,weight,1,regul);
	product(&fitted_energy,corr_matrix,cur_eci_mult);
	sum(&fitted_energy,fitted_energy,predictor_energy);
	// check for qualitatively wrong ground states;
	Array<int> problems;
	int gs_ok=!ground_state_problem(&problems, true_gs,concentration,fitted_energy,minc_gs_ok,maxc_gs_ok);
	
	// save best choice so far;
	cerr << gs_ok << " " << cv << endl;
	if ((gs_ok>best_gs_ok && cv!=MAXFLOAT) || (gs_ok==best_gs_ok && cv<best_cv)) {
	  best_gs_ok=gs_ok;
	  best_eci_mult=cur_eci_mult;
	  best_cv=cv;
	  get_cluster_choice(&best_choice);
	}
      } while (find_next_multiplet_choice());
      if (best_choice_noweight.get_size()==0) {
	best_choice_noweight=best_choice;
	best_eci_mult_noweight=best_eci_mult;
      }
      if (best_choice.get_size()==0) {
	if (pfitinfo) {pfitinfo->status=CEFitInfo::fit_impossible;}
	return;
      }
      // do regression and predict energies with best choice found;
      set_cluster_choice(best_choice);
      calc_regression_matrices(&corr_matrix,NULL);
      Array<Real> fitted_energy;
      product(&fitted_energy,corr_matrix,best_eci_mult);
      sum(&fitted_energy,fitted_energy,predictor_energy);
      // check for qualitatively wrong ground states;
      Array<int> problems;
      if (!ground_state_problem(&problems, true_gs,concentration,fitted_energy,minc_gs_ok,maxc_gs_ok)) break; // no problems: exit loop;
      // double the weights of all structures that cause problem;
      Real max_weight=0.;
      for (int i=0; i<problems.get_size(); i++) {
	if (problems(i)) weight(i)=weight(i)+1.;
	if (weight(i)>max_weight) max_weight=weight(i);
      }
      if (max_weight>=16.) { // give up after too many retries;
	one_array(&weight); // reset weights to one;
	best_choice=best_choice_noweight;
	best_eci_mult=best_eci_mult_noweight;
	if (pfitinfo) {pfitinfo->status = CEFitInfo::gs_problem;}
	break;
      }
    }
  }
  set_cluster_choice(best_choice); // save best choice into ClusterExpansion object;
  eci_mult=best_eci_mult; // along with best eci;

  // save all the information in a CEFitInfo object; 
  if (pfitinfo) {
    Real volume=(Real)(lattice_only.atom_pos.get_size());
    calc_regression_matrices(&corr_matrix,&energy);
    calc_regularization_vector(&regul,regul_param);

    // create array of pointers to structures included in the fit;
    pfitinfo->pstr.resize(energy.get_size());
    LinkedListIterator<StructureInfo> s(structures.get_structure_list());
    for (int i=0; i<pfitinfo->pstr.get_size(); s++) {
      if (s->status==StructureInfo::calculated) {
        pfitinfo->pstr(i)=s;
        i++;
      }
    }
    // find (or read in) the energy of the end members;
    pfitinfo->pure_energy.resize(2);
    pfitinfo->pure_energy(0)=MAXFLOAT;
    pfitinfo->pure_energy(1)=MAXFLOAT;
    {
      ifstream pure_file("ref_energy.in");
      if (pure_file) {
	pure_file >> pfitinfo->pure_energy(0);
	pure_file >> pfitinfo->pure_energy(1);
	pfitinfo->pure_energy(0)*=volume;
	pfitinfo->pure_energy(1)*=volume;
      }
    }
    // consider formation energies;
    calc_formation(&energy,concentration,energy,&(pfitinfo->pure_energy(0)),&(pfitinfo->pure_energy(1)));
    diff(&ce_energy,energy,predictor_energy);
    // copy concentration, energy, weight, list of true ground states;
    pfitinfo->concentration.copy(concentration);
    pfitinfo->energy.copy(energy);
    pfitinfo->weight.copy(weight);
    pfitinfo->gs_list.copy(true_gs);
     // volume of the largest generated structure;
    pfitinfo->max_volume=structures.get_max_volume();

    // eci (will be divided by multiplicity below);
    calc_ols_regul(&(pfitinfo->eci), corr_matrix,ce_energy,weight,1,regul);
    // predicted energies;
    product(&(pfitinfo->fitted_energy), corr_matrix,pfitinfo->eci);
    sum(&(pfitinfo->fitted_energy), pfitinfo->fitted_energy,predictor_energy);
    // predicted ground state line (among structures whose energy is known) using fitted energy (not used later on);
    calc_convex_hull(&(pfitinfo->fitted_gs_list), concentration,pfitinfo->fitted_energy);
    // cv score;
    pfitinfo->cv=calc_cv_regul(corr_matrix,ce_energy,weight,1,regul)/volume;
    // concentration range where gs must be ok;
    pfitinfo->minc_gs_ok=minc_gs_ok;
    pfitinfo->maxc_gs_ok=maxc_gs_ok;

    // done with fitting, we can express energies per atom;
    for (int i=0; i<pfitinfo->energy.get_size(); i++) {
      pfitinfo->energy(i)/=volume;
      pfitinfo->fitted_energy(i)/=volume;
    }
    pfitinfo->pure_energy(0)/=volume;
    pfitinfo->pure_energy(1)/=volume;
    {
      // write out reference energies of pure end members
      ofstream pure_file("ref_energy.out");
      pure_file << pfitinfo->pure_energy(0) << endl;
      pure_file << pfitinfo->pure_energy(1) << endl;
    }

    // look for predicted ground states among all generated structures;
    // of unknown energies;
    {
      // construct ground state line by joining the fitted energies of the
      // true ground states;
      Array<Real> gs_fitted_energy;
      Array<Real> gs_concentration;
      convex_hull_from_gs(&gs_concentration,&gs_fitted_energy, pfitinfo->gs_list,concentration,pfitinfo->fitted_energy);
      
      // predict energies of all generated structures of unknown energy;
      Array2d<Real> all_corr;
      Array<Real> all_fitted_energy;
      Array<Real> all_predictor_energy;
      calc_regression_matrices(&all_corr,NULL, (StructureInfo::Status) (StructureInfo::unknown | StructureInfo::busy | StructureInfo::error));
      calc_predictor_energy(&all_predictor_energy, (StructureInfo::Status) (StructureInfo::unknown | StructureInfo::busy | StructureInfo::error));
      product(&all_fitted_energy, all_corr,pfitinfo->eci);
      sum(&all_fitted_energy, all_fitted_energy,all_predictor_energy);
      product(&all_fitted_energy, all_fitted_energy,1./volume);
      // find concentration of all generated structures;
      Array<Real> all_concentration;
      calc_concentration(&all_concentration,all_corr);
      // check for structures breaking the hull;
      Array<int> problem;
      if (ground_state_problem(&problem, gs_concentration,gs_fitted_energy,all_concentration,all_fitted_energy,minc_gs_ok,maxc_gs_ok)) {
        pfitinfo->status = (CEFitInfo::Status)(pfitinfo->status | CEFitInfo::new_gs);
      }
      // output all conc. and energy of all generated structures of unknown energy;
      {
	LinkedListIterator<StructureInfo> s(structures.get_structure_list());
        ofstream allstr("predstr.out");
	allstr.setf(ios::fixed);
	allstr.precision(6);
        for (int i=0; i<all_fitted_energy.get_size(); i++) {
	  while (!(s->status & (StructureInfo::Status) (StructureInfo::unknown | StructureInfo::busy | StructureInfo::error))) {s++;}
	  const char *slabel=s->label;
          allstr << (1.+all_concentration(i))/2. << " " << all_fitted_energy(i) << " " 
		 << (s->status & StructureInfo::unknown ? "?" : slabel) << " " 
		 << (s->status & StructureInfo::busy    ? "b" : "") 
		 << (s->status & StructureInfo::error   ? "e" : "") 
		 << (s->status & StructureInfo::unknown ? "u" : "") 
		 << (problem(i) ? "g" : "") 
		 << endl;
	  s++;
        }
      }
    }
    // save clusters and their multiplicity;
    pfitinfo->cluster.resize(nb_clusters);
    pfitinfo->multiplicity.resize(nb_clusters);
    pfitinfo->multiplicity(0)=1.;
    pfitinfo->cluster(0).resize(0);
    int c=1;
    for (int m=0; m<best_choice.get_size(); m++) {
      pclusters(m)->reset();
      for (int i=0; i<best_choice(m); i++, (*pclusters(m))++, c++) {
        pfitinfo->cluster(c).copy((*pclusters(m)));
        pfitinfo->multiplicity(c)=(Real)calc_multiplicity(pfitinfo->cluster(c),parent_lattice.cell,spacegroup.point_op,spacegroup.trans);
        pfitinfo->eci(c)/=pfitinfo->multiplicity(c);
      }
    }
  }
}

StructureInfo * CE_TLambda::find_best_structure(void) {
  // do the fit;
  Array2d<Real> raw_corr,corr;
  Array<Real> energy;
  Array<Real> ce_energy;
  Array<Real> predictor_energy;
  Array<Real> eci;
  Array<Real> fitted_energy;
  Array<Real> raw_regul,regul;
  calc_regularization_vector(&raw_regul,cur_regul_param); // new for tlambda;
  calc_regression_matrices(&raw_corr,&energy);
  calc_predictor_energy(&predictor_energy);
  diff(&ce_energy,energy,predictor_energy);
  Array<int> cols;
  list_nonredundant_columns(&cols,raw_corr);
  extract_columns(&corr,raw_corr,cols);
  extract_elements(&regul,raw_regul,cols);

  if (corr.get_size()(1)<3 || corr.get_size()(0)<corr.get_size()(1)) return NULL;
  if (weight.get_size()!=energy.get_size()) {
    weight.resize(energy.get_size());
    one_array(&weight);
  }
  
  calc_ols_regul(&eci, corr,ce_energy,weight,1,regul); // new for tlambda;
  product(&fitted_energy, corr,eci);
  sum(&fitted_energy, fitted_energy, predictor_energy);

  // find ground state line;
  Array<Real> concentration;
  calc_concentration(&concentration,corr);
  Array<int> true_gs;
  Array<Real> gs_fitted_energy;
  Array<Real> gs_concentration;
  calc_convex_hull(&true_gs, concentration,energy);
  convex_hull_from_gs(&gs_concentration,&gs_fitted_energy, true_gs,concentration,fitted_energy);

  // predict all unknown energies;
  Array2d<Real> raw_un_corr,un_corr;
  Array<Real> un_fitted_energy;
  Array<Real> un_predictor_energy;
  calc_regression_matrices(&raw_un_corr,NULL, StructureInfo::unknown);
  calc_predictor_energy(&un_predictor_energy, StructureInfo::unknown);
  extract_columns(&un_corr,raw_un_corr,cols);

  product(&un_fitted_energy, un_corr,eci);
  sum(&un_fitted_energy, un_fitted_energy,un_predictor_energy);
  // and concentrations;
  Array<Real> un_concentration;
  calc_concentration(&un_concentration,un_corr);

  // are there ground state problems?;
  Array<int> problem; // will contain 1's for structures which break the hull;
  if (ground_state_problem(&problem, gs_concentration,gs_fitted_energy,un_concentration,un_fitted_energy,minc_gs_ok,maxc_gs_ok)) {
    // if there are new predicted ground state, look among them for the structure that give the greatest var reduction;
    Array2d<Real> rhorho;
    // treat all structures as equally important;
    set_identity(&rhorho,corr.get_size()(1));

    Array<int> new_gs;
    calc_convex_hull(&new_gs,un_concentration,un_fitted_energy);
    for (int i=0; i<new_gs.get_size(); i++) {
      problem(new_gs(i))|=2; // structures that break the hull and are potential new ground states will be labelled 3;
    }
    
    // compute the matrices that enter the expression
    // of the covariance matrix of the ECI;
    Array2d<Real> w_x,ww_x;
    product_diag(&w_x, weight,corr);
    product_diag(&ww_x, weight,w_x);
    Array2d<Real> xw2x,xw4x;
    inner_product(&xw2x, w_x,w_x);
    for (int k=0; k<regul.get_size(); k++) {xw2x(k,k)+=regul(k);} // new for tlambda;
    inner_product(&xw4x, ww_x,ww_x);
    Array2d<Real> xw2xi;
    invert_matrix(&xw2xi, xw2x);
    // compute variance before adding a new structure;
    Real var_before; 
    {
      Array2d<Real> tmp,var_eci;
      product(&tmp, xw4x,xw2xi);
      product(&var_eci, xw2xi,tmp);
      product(&tmp, var_eci, rhorho);
      var_before=trace(tmp);
    }

    // save best structure;
    StructureInfo *p_best_str=NULL;
    Real max_benefit=0.; // benefit is variance reduction over computational cost;
    reset_structure();
    // loop over predicted ground states whose true energy is unknown;
    for (int istr=0; istr<problem.get_size(); istr++) {
      while (structures.get_current_structure().status != StructureInfo::unknown) {find_next_structure();}
      if (problem(istr)==3) { // if new structure breaks the hull and is a potential new ground state;
        // calculate variance after the current structure has been added to the fit;
        Real var;
        {
          Array<Real> raw_cur_rho,cur_rho;
          make_correlation_array(&raw_cur_rho, structures.get_current_structure().correlations);
          extract_elements(&cur_rho,raw_cur_rho,cols);
          Array2d<Real> cur_rhorho;
          outer_product(&cur_rhorho, cur_rho,cur_rho);
          Array2d<Real> cur_xw2x,cur_xw2xi,cur_xw4x;
          sum(&cur_xw2x, xw2x,cur_rhorho);
          sum(&cur_xw4x, xw4x,cur_rhorho);
          invert_matrix(&cur_xw2xi, cur_xw2x);
          Array2d<Real> tmp,var_eci;
          product(&tmp, cur_xw4x,cur_xw2xi);
          product(&var_eci, cur_xw2xi,tmp);
          product(&tmp, var_eci, rhorho);
          var=trace(tmp);
        }
        // calculate benefit of adding structure;
        Real benefit=(var_before-var)/structures.get_current_structure().cost;
        // save structure if it is the best so far;
        if (benefit>max_benefit) {
          max_benefit=benefit;
          p_best_str=&structures.get_current_structure();
        }
      }
      if (istr<problem.get_size()-1) find_next_structure();
    } // loop to next structure;
    return p_best_str;
  }
  else {
    // no ground state problems: look among all unknown structures for the most variance-reducing one;
    Array2d<Real> rhorho;
     // treat all structures as equally important;
    set_identity(&rhorho,corr.get_size()(1));

    // compute the matrices that enter the expression
    // of the covariance matrix of the ECI;
    Array2d<Real> w_x,ww_x;
    product_diag(&w_x, weight,corr);
    product_diag(&ww_x, weight,w_x);
    Array2d<Real> xw2x,xw4x;
    inner_product(&xw2x, w_x,w_x);
    for (int k=0; k<regul.get_size(); k++) {xw2x(k,k)+=regul(k);} // new for tlambda;
    inner_product(&xw4x, ww_x,ww_x);
    Array2d<Real> xw2xi;
    invert_matrix(&xw2xi, xw2x);
    
    // compute best possible correlation to add;
    Array<Real> best_dir(corr.get_size()(1));
    {
      // best_dir is the eigenvector associated with the smallest eigenvalue of xw2x;
      // find by iterating u=inverse(xw2x)*v; v=u/|u|;
      zero_array(&best_dir); 
      best_dir(best_dir.get_size()-1)=1.; // initial guess;
      Array<Real> old_dir;
      Array<Real> d_dir;
      do { // iterate;
        old_dir.copy(best_dir);
        Array<Real> tmp;
        product(&tmp, xw2xi,best_dir);
        product(&best_dir, tmp, 1./sqrt(inner_product(tmp,tmp)));
        diff(&d_dir ,best_dir,old_dir);
      } while (inner_product(d_dir,d_dir)>1e-4); // until tolerance met;
      // normalize so that largest correlation has modulus 1;
      Real m=MAX(fabs(max(best_dir)),fabs(min(best_dir)));
      product(&best_dir,best_dir,1./m);
    }
    
    // compute variance before adding a new structure;
    Real var_before; 
    {
      Array2d<Real> tmp,var_eci;
      product(&tmp, xw4x,xw2xi);
      product(&var_eci, xw2xi,tmp);
      product(&tmp, var_eci, rhorho);
      var_before=trace(tmp);
    }
    
    // compute variance after the best possible correlation has been added;
    Real var;
    {
      Array2d<Real> cur_rhorho;
      outer_product(&cur_rhorho, best_dir, best_dir);
      Array2d<Real> cur_xw2x,cur_xw2xi,cur_xw4x;
      sum(&cur_xw2x, xw2x,cur_rhorho);
      sum(&cur_xw4x, xw4x,cur_rhorho);
      invert_matrix(&cur_xw2xi, cur_xw2x);
      Array2d<Real> tmp,var_eci;
      product(&tmp, cur_xw4x,cur_xw2xi);
      product(&var_eci, cur_xw2xi,tmp);
      product(&tmp, var_eci, rhorho);
      var=trace(tmp);
    }
    
    // now we have the maximum possible variance reduction due to adding a structure;
    Real max_possible_dvar=var_before-var;
    
    // save best structure;
    StructureInfo *p_best_str=NULL;
    Real max_benefit=0.; // benefit is variance reduction over computational cost;
    reset_structure();
    // loop over unknown structures;
    while (1) {
      if (structures.get_current_structure().status==StructureInfo::unknown) {
        // stop search when current computational cost is so high that even if;
        // the maximum variance reduction is achieved the benefits are lower;
        // than what we already have;
        if (max_possible_dvar/structures.get_current_structure().cost < max_benefit) break;
        // calculate variance after the current structure has been added to the fit;
        Real var;
        Real cur_concentration;
        {
          Array<Real> raw_cur_rho,cur_rho;
          make_correlation_array(&raw_cur_rho, structures.get_current_structure().correlations);
          Array<Real> point_multiplicity;
          calc_point_to_concentration(&point_multiplicity, raw_cur_rho.get_size());
          cur_concentration=inner_product(raw_cur_rho,point_multiplicity);
          extract_elements(&cur_rho,raw_cur_rho,cols);
          Array2d<Real> cur_rhorho;
          outer_product(&cur_rhorho, cur_rho,cur_rho);
          Array2d<Real> cur_xw2x,cur_xw2xi,cur_xw4x;
          sum(&cur_xw2x, xw2x,cur_rhorho);
          sum(&cur_xw4x, xw4x,cur_rhorho);
          invert_matrix(&cur_xw2xi, cur_xw2x);
          Array2d<Real> tmp,var_eci;
          product(&tmp, cur_xw4x,cur_xw2xi);
          product(&var_eci, cur_xw2xi,tmp);
          product(&tmp, var_eci, rhorho);
          var=trace(tmp);
        }
        // calculate benefit of adding structure;
        Real benefit=(var_before-var)/structures.get_current_structure().cost;
        // save structure if it is the best so far;
        if (benefit>max_benefit &&
            cur_concentration>=minc_gs_ok-zero_tolerance &&
            cur_concentration<=maxc_gs_ok+zero_tolerance) {
          max_benefit=benefit;
          p_best_str=&structures.get_current_structure();
        }
      }
      find_next_structure();
    } // loop to next structure;
    return p_best_str;
  }
}

void CE_TLambda::calc_regularization_vector(Array<Real> *regul, const Array<Real> &param) {
  regul->resize(nb_clusters);
  (*regul)(0)=0.;
  int i=1;
  for (int m=0; m<pclusters.get_size(); m++) {
    LinkedListIterator<Cluster> ic(get_cluster_list(m));
    LinkedListIterator<Array<Cluster> > ieqc(*equiv_cluster(m));
    for (int index=0; index < pclusters(m)->get_current_index(); ic++, ieqc++, index++, i++) {
      if (m==1) {
	(*regul)(i)=param(0)*pow(get_length_quick(*ic),param(1))/(ieqc->get_size());
      }
      else {
	(*regul)(i)=0;
      }
    }
  }
}

int CE_TLambda::find_next_multiplet_choice(void) {
  // add new multiplets if needed;
  if (pclusters(pclusters.get_size()-1)->get_current_index()>0 && pclusters.get_size()<max_multiplet) {
    Array<ClusterBank *> tmp(pclusters);
    pclusters.resize(tmp.get_size()+1);
    for (int i=0; i<tmp.get_size(); i++) {
      pclusters(i)=tmp(i);
    }
    pclusters(tmp.get_size())=new ClusterBank(*(pclusters(1)),tmp.get_size()+1);
  }

  // main algorithm to find next cluster choice;
  //  A) go to the multiplet with the most points;
  //  B) add new larger clusters which have the number of point;
  //     (clusters of the same diameter are simultaneously added);
  //  C) if the largest cluster is a triplet, no more cluster choices exist;
  //     exit;
  //  D) if;
  //      -the diameter of the added clusters is smaller or equal to;
  //       the one of other clusters with less points;
  //      -and the total number of clusters is less than the nb of structures;
  //  E) then exit loop: we have found the next cluster choice;
  //  F) otherwise, get rid of all clusters with that number of points;
  //  G) repeat the process with a clusters having one less point;
  int c=1;
  // A;
  while (c<pclusters.get_size()-1 && pclusters(c)->get_current_index()>0) {c++;}
  // B;
  while (1) {
    Real cur_len=pclusters(c)->get_current_length();
    while (pclusters(c)->get_current_length()<cur_len+zero_tolerance) {
      (*pclusters(c))++;
      nb_clusters++;
    }
    // C;
    if (c==2) break;
    // D;
    if (nb_clusters < nb_calculated_str &&
        pclusters(c)->get_previous_length() < pclusters(c-1)->get_previous_length()+zero_tolerance) break; // E;
    // F;
    nb_clusters-=pclusters(c)->get_current_index();
    pclusters(c)->reset();
    // G;
    c--;
  }

  // update the lists of equivalent clusters;
  update_equivalent_clusters();
  // make sure that the corresponding corr. have been calculated;
  update_correlations(StructureInfo::calculated);
  // are there any other cluster choices to try? (check if number of non pair is less than number of structures);
  return (nb_clusters-(pclusters(1)->get_current_index()) < nb_calculated_str);
}
