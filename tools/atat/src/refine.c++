#include <fstream.h>
#include <strstream.h>
#include "refine.h"
#include "getvalue.h"

GenericPlugIn<EnergyPredictor> *GenericPlugIn<EnergyPredictor>::list=NULL;
GenericPlugIn<ClusterExpansionCreator> *GenericPlugIn<ClusterExpansionCreator>::list=NULL;
SpecificPlugIn<ClusterExpansionCreator, ClusterExpansionCreator> StdClusterExpansionPlugIn("std"); // register plug-in under name "std";


Real calc_structure_cost(const Structure &str, const rMatrix3d &cell, const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans, Real complexity_exp) {
  return pow((Real)(str.atom_pos.get_size()),complexity_exp);
}

ClusterExpansion::ClusterExpansion(const Structure &_parent_lattice,
                                   const Array<Arrayint> &_site_type_list,
                                   const Array<AutoString> &_atom_label,
                                   const SpaceGroup &_spacegroup):
   pclusters(), equiv_cluster(), spacegroup(_spacegroup), structures(_parent_lattice,_site_type_list,spacegroup), weight(), 
   parent_lattice(_parent_lattice), site_type_list(_site_type_list), atom_label(_atom_label), label_list(), i_label_list() {
     max_multiplet=4;
     complexity_exp=3.;
     minc_gs_ok=-MAXFLOAT;
     maxc_gs_ok=MAXFLOAT;
     predictor_labels="";

     // initialize lattice_only;
     int nbsite=0;
     for (nbsite=0; nbsite<parent_lattice.atom_pos.get_size(); nbsite++) {
       if (site_type_list(parent_lattice.atom_type(nbsite)).get_size()==1) break;
     }
     lattice_only.cell=parent_lattice.cell;
     lattice_only.atom_pos.resize(nbsite);
     lattice_only.atom_type.resize(nbsite);
     for (int i=0; i<lattice_only.atom_pos.get_size(); i++) {
       lattice_only.atom_pos(i) =parent_lattice.atom_pos(i);
       lattice_only.atom_type(i)=parent_lattice.atom_type(i);
     }

     pclusters.resize(2);
     // initialize point clusters;
     pclusters(0)=new ClusterBank(lattice_only.cell,lattice_only.atom_pos,1,spacegroup);
     int npoint=pclusters(0)->get_cluster_list().get_size();
     for (int i=0; i<npoint; i++) {
       (*pclusters(0))++;
     }

     // initialize 1nn pair clusters;
     radius_1nn=find_1nn_radius(lattice_only)+zero_tolerance;
     pclusters(1)=new ClusterBank(lattice_only.cell,lattice_only.atom_pos,2,spacegroup);
     while (pclusters(1)->get_current_length()<radius_1nn) {(*pclusters(1))++;}

     // initialize equivalent clusters lists;
     equiv_cluster.resize(2);
     for (int i=0; i<pclusters.get_size(); i++) {
       equiv_cluster(i)=new LinkedListArrayCluster;
       LinkedListIterator<Cluster> c(pclusters(i)->get_cluster_list());
       for ( ; c; c++) {
         Array<Cluster> *parraycluster=new Array<Cluster>;
         find_equivalent_clusters(parraycluster,*c,spacegroup.cell,spacegroup.point_op,spacegroup.trans);
         (*equiv_cluster(i)) << parraycluster;
       }
     }
     reset_cluster_choice();
     update_correlations(StructureInfo::anystatus);
     reset_structure();
}

ClusterExpansion::~ClusterExpansion(void) {
  for (int i=0; i<pclusters.get_size(); i++) {
    delete pclusters(i);
    delete equiv_cluster(i);
  }
}

void ClusterExpansion::update_correlations(StructureInfo *str) {
  // set structure label, if needed;
  if (strlen(str->label)==0) {
    while (i_label_list && *i_label_list == cur_label) {
      i_label_list++;
      cur_label++;
    }
    label_list.add(new int(cur_label),i_label_list);
    ostrstream tmp;
    tmp << cur_label << '\0';
    str->label.set(tmp.str());
  }
  // set structure cost, if needed;
  if (str->cost==0) str->cost=calc_structure_cost(*str,spacegroup.cell,spacegroup.point_op,spacegroup.trans,complexity_exp);

  // creates new lists for larger multiplet if needed;
  int multiplet_to_add=pclusters.get_size() - str->correlations.get_size();
  while (multiplet_to_add>0) {
    str->correlations << new LinkedList<Real>();
    multiplet_to_add--;
  }

  // loop through multiplets and add missing correlations;
  LinkedListIterator<LinkedListReal> corr_list(str->correlations);
  for (int clus_list=0; clus_list<pclusters.get_size(); clus_list++, corr_list++) {
    // for each multiplet, prepare to loop through clusters;
    LinkedListIterator<Real> corr(*corr_list);
    LinkedListIterator<Cluster> clus(pclusters(clus_list)->get_cluster_list());
    LinkedListIterator<ArrayCluster> equiv_clus(*equiv_cluster(clus_list));
    int idx=0;
    // go to the end of the corr. list;
    while (corr) {corr++; clus++; equiv_clus++; idx++;}\
    // and then add the missing corr.;
    for ( ; idx<pclusters(clus_list)->get_current_index(); clus++, equiv_clus++, idx++) {
      (*corr_list) << new Real(calc_correlation(*str,*equiv_clus,spacegroup.cell));
    }
  }

  // create predictor objects, if needed;
  if (str->predictor.get_size()==0) {
    make_plug_in_list(&(str->predictor),predictor_labels);
  }
}

void ClusterExpansion::update_correlations(StructureInfo::Status select) {
  reset_structure();
  // for each structure that maches the select mask, update correlations.;
  LinkedListIterator<StructureInfo> str(structures.get_structure_list());
  for ( ; str; str++) {
    if (str->status & select) update_correlations(str);
  }
}

void ClusterExpansion::init_predictors(void) {
  LinkedList<EnergyPredictor> l;
  make_plug_in_list(&l,predictor_labels);
  LinkedListIterator<EnergyPredictor> i(l);
  for (; i; i++) {
    i->static_init(*this);
  }
}

void ClusterExpansion::calc_point_to_concentration(Array<Real> *multiplicity, int nb_cols) {
  // figure out the weight to give to each point corr to get the concentration;
  LinkedListIterator<Cluster> i(pclusters(0)->get_cluster_list());
  // figure out the weight to give to each point corr to get the concentration;
  multiplicity->resize(nb_cols);
  zero_array(multiplicity);
  for (int j=1; i; i++,j++) {
    (*multiplicity)(j)=(Real)calc_multiplicity(*i,lattice_only.cell,spacegroup.point_op,spacegroup.trans)/
                    (Real)lattice_only.atom_pos.get_size();
  }
}

void ClusterExpansion::calc_concentration(Array<Real> *pconcentration, const Array2d<Real> &corr_matrix) {
  // figure out the weight to give to each point corr to get the concentration;
  Array<Real> multiplicity;
  calc_point_to_concentration(&multiplicity,corr_matrix.get_size()(1));
  product(pconcentration,corr_matrix,multiplicity);
}

void ClusterExpansion::calc_regression_matrices(Array2d<Real> *pcorr_matrix, Array<Real> *penergy,
        StructureInfo::Status select) {
  // update all selected structures;
  update_correlations(select);
  // first count the number of structures selected and resize relevant arrays;
  int nb_str=0;
  LinkedListIterator<StructureInfo> str(structures.get_structure_list());
  for (; str; str++) {
    if (str->status & select) nb_str++;
  }
  pcorr_matrix->resize(iVector2d(nb_str,nb_clusters));
  if (penergy) {penergy->resize(nb_str);}

  // loop through structures;
  int i_str=0;
  str.init(structures.get_structure_list());
  for (; str; str++) {
    if (str->status & select) { //if structure matches mask;
      if (penergy) {(*penergy)(i_str)=str->energy;} // copy energy;
      // copy correlations;
      (*pcorr_matrix)(i_str,0)=1.;
      int i_corr=1;
      LinkedListIterator<LinkedListReal> corr_list(str->correlations);
      for (int clus_list=0; corr_list; corr_list++, clus_list++) {
        LinkedListIterator<Real> corr(*corr_list);
        for (int i=0; i<pclusters(clus_list)->get_current_index(); i++) {
          (*pcorr_matrix)(i_str,i_corr)=*corr;
          corr++;
          i_corr++;
        }
      }
      i_str++;
    }
  }
}

void ClusterExpansion::calc_predictor_energy(Array<Real> *penergy, StructureInfo::Status select) {
  // update all selected structures;
  update_correlations(select);
  // first count the number of structures selected and resize relevant arrays;
  int nb_str=0;
  LinkedListIterator<StructureInfo> str(structures.get_structure_list());
  for (; str; str++) {
    if (str->status & select) nb_str++;
  }
  penergy->resize(nb_str);

  // loop through structures;
  int i_str=0;
  str.init(structures.get_structure_list());
  for (; str; str++) {
    if (str->status & select) { //if structure matches mask;
      (*penergy)(i_str)=0.;
      LinkedListIterator<EnergyPredictor> ip(str->predictor); // loop through user-defined predictor;
      for (; ip; ip++) {
	(*penergy)(i_str)+=ip->get_energy(*str);
      }
      // cout << str->label << " " << (*penergy)(i_str) << endl; // trace info, if needed;
      i_str++;
    }
  }
}

void ClusterExpansion::reset_cluster_choice(void) {
  // reset all cluster generators;
  for (int i=1; i<pclusters.get_size(); i++) {
    pclusters(i)->reset();
  }
  // include minimum number of pairs (1nn);
  while ( pclusters(1)->get_current_length()<radius_1nn ) {(*(pclusters(1)))++;}

  // initialize nb of structures in fit;
  nb_calculated_str=0;
  LinkedListIterator<StructureInfo> str(structures.get_structure_list());
  for (; str; str++) {
    if (str->status==StructureInfo::calculated) nb_calculated_str++;
  }
  // init nb of clusters in fit;
  nb_clusters=1+pclusters(0)->get_current_index()+pclusters(1)->get_current_index();
  // make sure we have enough corr;
  update_correlations(StructureInfo::calculated);
}

void ClusterExpansion::update_equivalent_clusters(void) {
  // update the lists of clusters equivalent by a space group symmetry operation (to speed up the calculation of the correlations);
  //  first add new multiplets if needed;
  if (pclusters.get_size()>equiv_cluster.get_size()) {
    Array<LinkedListArrayCluster *> tmp;
    tmp.copy(equiv_cluster);
    equiv_cluster.resize(pclusters.get_size());
    int m;
    for (m=0; m<tmp.get_size(); m++) {
      equiv_cluster(m)=tmp(m);
    }
    for (; m<pclusters.get_size(); m++) {
      equiv_cluster(m)=new LinkedListArrayCluster;
    }
  }
  // now add the missing equivalent clusters;
  for (int m=0; m<equiv_cluster.get_size(); m++) {
    LinkedListIterator<Cluster> c(pclusters(m)->get_cluster_list());
    LinkedListIterator<ArrayCluster> e(*equiv_cluster(m));
    while (e) {e++; c++;} // go to the end of the list;
    for (; c; c++) { // loop through clusters;
      Array<ArrayrVector3d> *pequiv=new Array<Cluster>;
      find_equivalent_clusters(pequiv,*c,spacegroup.cell,spacegroup.point_op,spacegroup.trans);
      (*equiv_cluster(m)) << pequiv;
    }
  }
}

int ClusterExpansion::find_next_cluster_choice(void) {
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
  //  C) if the largest cluster is a pair, no more cluster choices exist;
  //     exit;
  //  D) if;
  //      -the diameter of the added clusters is smaller or equal to;
  //       the one of other clusters with less points;
  //      -and the total number of clusters is less than the nb of structures;
  //  E) then exit loop: we have found the next cluster choice;
  //  F) otherwise, get rid of all clusters with that number of points;
  //  G) repeat the process with a clusters having one less point;
  int c=0;
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
    if (c==1) break;
    // D;
    if (nb_clusters < nb_calculated_str &&
        pclusters(c)->get_previous_length() < pclusters(c-1)->get_previous_length()+zero_tolerance) break; // E;
    // F;
    nb_clusters-=pclusters(c)->get_current_index();
    pclusters(c)->reset();
    // G;
    c--;
  }

  // update the lists of clusters equivalent under space group symmetry operation;
  update_equivalent_clusters();
  // make sure that the corresponding corr. have been calculated;
  update_correlations(StructureInfo::calculated);
  // are there any other cluster choices to try?;
  return (nb_clusters < nb_calculated_str);
}

void ClusterExpansion::set_cluster_choice(const Array<int> &choice) {
  // add new multiplets if needed;
  if (pclusters.get_size()<choice.get_size()) {
    Array<ClusterBank *> tmp(pclusters);
    pclusters.resize(choice.get_size());
    int m;
    for (m=0; m<tmp.get_size(); m++) {
      pclusters(m)=tmp(m);
    }
    for (; m<choice.get_size(); m++) {
      pclusters(m)=new ClusterBank(*(pclusters(1)),m+1);
    }
  }

  nb_clusters=1+pclusters(0)->get_current_index(); // count empty and point;
  // advance the index of cluster generator i to choice(i);
  int i;
  for (i=1; i<choice.get_size(); i++) {
    pclusters(i)->reset();
    while (pclusters(i)->get_current_index()<choice(i)) (*pclusters(i))++;
    nb_clusters+=pclusters(i)->get_current_index(); // keep count of clusters;
  }
  // reset any other multiplet not given in array choice;
  for (; i<pclusters.get_size(); i++) {
    pclusters(i)->reset();
  }
  update_equivalent_clusters();
}

void ClusterExpansion::get_cluster_choice(Array<int> *pchoice) {
  pchoice->resize(pclusters.get_size());
  for (int i=0; i<pchoice->get_size(); i++) {
    (*pchoice)(i)=pclusters(i)->get_current_index();
  }
}

void ClusterExpansion::find_best_cluster_choice(CEFitInfo *pfitinfo) {
  reset_cluster_choice(); // select minimal cluster expansion;

  // setup correlation matrix, and energy and concentration vectors;
  Array2d<Real> corr_matrix;
  Array<Real> energy;
  Array<Real> predictor_energy;
  Array<Real> ce_energy;
  Array<Real> concentration;
  calc_regression_matrices(&corr_matrix,&energy);
  if (corr_matrix.get_size()(0)<=corr_matrix.get_size()(1)) { // enough structures?;
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
  // if a user specified choice has been given;
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
      reset_cluster_choice(); // select minimal cluster expansion;
      int best_gs_ok=0;
      Real best_cv=MAXFLOAT;
      do { // loop through cluster choices;
	//if you want to print out trace information;
	for (int i=0; i<pclusters.get_size(); i++) {
	  cerr << pclusters(i)->get_current_index() << " ";
	}
	
	// calc cross-validation score;
	calc_regression_matrices(&corr_matrix,&energy);
	calc_predictor_energy(&predictor_energy);
	diff(&ce_energy,energy,predictor_energy);
	Real cv=calc_cv(corr_matrix,ce_energy,weight,1);
	// do regression and add back the predictor energy;
	Array<Real> cur_eci_mult;
	Array<Real> fitted_energy;
	calc_ols(&cur_eci_mult,corr_matrix,ce_energy,weight,1);
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
	if (cv==MAXFLOAT) {
	  Array<int> cur_choice;
	  get_cluster_choice(&cur_choice);
	  int i=2;
	  for ( ; i<cur_choice.get_size(); i++) {
	    if (cur_choice(i)!=0) break;
	  }
	  if (i==cur_choice.get_size()) break;
	}
      } while (find_next_cluster_choice());
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
      //cerr << "GS problems" << endl;
      // increase the weights of all structures that cause problem;
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
    calc_ols(&(pfitinfo->eci), corr_matrix,ce_energy,weight,1);
    // predicted energies;
    product(&(pfitinfo->fitted_energy), corr_matrix,pfitinfo->eci);
    sum(&(pfitinfo->fitted_energy), pfitinfo->fitted_energy,predictor_energy);
    // predicted ground state line (among structures whose energy is known) using fitted energy (not used later on);
    calc_convex_hull(&(pfitinfo->fitted_gs_list), concentration,pfitinfo->fitted_energy);
    Array<int> problems;
    if (ground_state_problem(&problems, pfitinfo->fitted_gs_list,concentration,pfitinfo->fitted_energy,minc_gs_ok,maxc_gs_ok)) {
      pfitinfo->status = CEFitInfo::gs_problem;
    }
      
    // cv score;
    pfitinfo->cv=calc_cv(corr_matrix,ce_energy,empty_rArray,1)/volume;
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
	  Real c=(1.+all_concentration(i))/2.;
	  Real form_e=s->energy/volume-(1-c)*pfitinfo->pure_energy(0)-c*pfitinfo->pure_energy(1);
          allstr << c << " " << (s->energy==0 ? 0 : form_e) << " " << all_fitted_energy(i) << " " 
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

void ClusterExpansion::make_correlation_array(Array<Real> *array_corr, const LinkedList<LinkedListReal> &list_corr) {
  array_corr->resize(nb_clusters);
  (*array_corr)(0)=1.;
  int dest=1;
  // loop through each multiplet
  LinkedListIterator<LinkedListReal> lc(list_corr);
  for (int i=0; i<pclusters.get_size(); i++, lc++) {
    // loop through each cluster;
    LinkedListIterator<Real> src(*lc);
    for (int c=0; c<pclusters(i)->get_current_index(); c++, dest++, src++) {
      (*array_corr)(dest)=*src;
    }
  }
}

void ClusterExpansion::reset_structure(void) {
  LinkedListIterator<int> i(label_list);
  while (i) delete label_list.detach(i);
  LinkedListIterator<StructureInfo> j(structures.get_structure_list());
  for ( ; j; j++) {
    int index=-1;
    istrstream tmp((const char *)(j->label));
    tmp >> index;
    if (index!=-1) {
      LinkedListIterator<int> i(label_list);
      for ( ; i; i++) {
	if (*i>index) break;
      }
      label_list.add(new int(index),i);
    }
  }
  i_label_list.init(label_list);
  cur_label=0;
  structures.reset();
  update_correlations(&(structures.get_current_structure()));
}

void ClusterExpansion::find_next_structure(void) {
  structures++;
  update_correlations(&(structures.get_current_structure()));
}

void ClusterExpansion::set_predictor_labels(const char *labels) {
  predictor_labels=labels;
  update_correlations(StructureInfo::anystatus);
}

StructureInfo * ClusterExpansion::find_best_structure(void) {
  // do the fit;
  Array2d<Real> raw_corr,corr;
  Array<Real> energy;
  Array<Real> ce_energy;
  Array<Real> predictor_energy;
  Array<Real> eci;
  Array<Real> fitted_energy;
  calc_regression_matrices(&raw_corr,&energy);
  calc_predictor_energy(&predictor_energy);
  diff(&ce_energy,energy,predictor_energy);
  Array<int> cols;
  list_nonredundant_columns(&cols,raw_corr);
  extract_columns(&corr,raw_corr,cols);

  if (corr.get_size()(1)<3 || corr.get_size()(0)<corr.get_size()(1)) return NULL;
  if (weight.get_size()!=energy.get_size()) {
    weight.resize(energy.get_size());
    one_array(&weight);
  }
  
  calc_ols(&eci, corr,ce_energy,weight);
  product(&fitted_energy, corr,eci);
  sum(&fitted_energy, fitted_energy, predictor_energy);

  // find ground state line;
  Array<Real> concentration;
  calc_concentration(&concentration,raw_corr);
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
  calc_concentration(&un_concentration,raw_un_corr);

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

StructureInfo * ClusterExpansion::find_initial_structures(void) {
  // find the space spanned by the correlations of existing structures;
  // basis contains the result;
  reset_cluster_choice();
  update_correlations(StructureInfo::anystatus);
  LinkedListIterator<StructureInfo> i(structures.get_structure_list());
  Array<ArrayReal> basis;
  int nb_in_basis=0;
  for (;i; i++) {
    if (i->status & (StructureInfo::calculated | StructureInfo::busy)) {
      Array<Real> cur_corr;
      make_correlation_array(&cur_corr,i->correlations);
      build_basis(&basis,&nb_in_basis,cur_corr);
    }
  }
  // if there are enough structures to fit the minimal cluster expansion, quit;
  if (nb_in_basis==nb_clusters) return NULL;

  // loop through unknown structures in search of a non-coplanar correlation;
  reset_structure();
  while (1) {
    if (structures.get_current_structure().status==StructureInfo::unknown) {
      Array<Real> cur_corr;
      make_correlation_array(&cur_corr,structures.get_current_structure().correlations);
      Array<Real> point_multiplicity;
      calc_point_to_concentration(&point_multiplicity, cur_corr.get_size());
      Real c=inner_product(cur_corr,point_multiplicity);
      if (c>=minc_gs_ok && c<=maxc_gs_ok) {
	// if correlation vector is not coplanar to the basis, we're done;
	if (build_basis(&basis,&nb_in_basis,cur_corr)) break;
      }
    }
    find_next_structure();
  }
  return &(structures.get_current_structure());
}

StructureInfo * ClusterExpansion::find_first_unknown_structure(void) {
  reset_structure();
  while (1) {
    while (structures.get_current_structure().status!=StructureInfo::unknown) {
      find_next_structure();
    }
    Array<Real> cur_corr;
    make_correlation_array(&cur_corr,structures.get_current_structure().correlations);
    Array<Real> point_multiplicity;
    calc_point_to_concentration(&point_multiplicity, cur_corr.get_size());
    Real c=inner_product(cur_corr,point_multiplicity);
    if (c>=minc_gs_ok && c<=maxc_gs_ok) break;
    find_next_structure();
  }
  return &(structures.get_current_structure());
}
