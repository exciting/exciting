#include <fstream.h>
#include <strstream.h>
#include "mmclib.h"
#include "parse.h"
#include "getvalue.h"
#include "lstsqr.h"
#include "version.h"
#include "mteci.h"
#include "equil.h"
#include "kmeci.h"
#include "plugin.h"
#include "calcmf.h"
#include "multipoly.h"
#include "chull.h"

extern char *helpstring;

GenericPlugIn<KSpaceECI> *GenericPlugIn<KSpaceECI>::list=NULL;

void print_vector(ostream &file, const Array<Real> &v, char *delim) {
  for (int i=0; i<v.get_size(); i++) {
    file << v(i) << delim;
  }
}

//ofstream debugfile("lrodebug.out");

Real linlro(Real lro) {
  return (log(1./lro-1.));
}

Real check_in_phase(const LinkedList<Array<Real> > &kept_control, const LinkedList<Real> &kept_lro, const Array<Real> &cur_control, Real cur_lro) {
  const int pow_limit=6;
  Array<Real> y;
  LinkedList_to_Array(&y,kept_lro);
  for (int i=0; i<y.get_size(); i++) {y(i)=linlro(y(i));}
  Real best_cv=MAXFLOAT;
  int best_pow=-1;
  Array2d<Real> best_x;
  LinearInequality ineq(cur_control.get_size());
  LinkedList<LinearInequality> ineqlist;
  ineqlist << &ineq;
  one_array(&ineq.v);
  Array<Real> max_pow;
  int nb_kept=kept_control.get_size();
  int cur_pow=0;
  while (1) {
    ineq.c=(Real)cur_pow;
    product(&max_pow,ineq.v,ineq.c);
    MultiDimPoly poly(max_pow,ineqlist);
    if (poly.get_dim_param()>(nb_kept-2) || cur_pow>pow_limit) break;
    Array<Real> param(poly.get_dim_param());
    zero_array(&param);
    Array2d<Real> x(nb_kept,poly.get_dim_param());
    LinkedListIterator<Array<Real> > ic(kept_control);
    for (int i=0; ic; i++, ic++) {
      for (int j=0; j<poly.get_dim_param(); j++) {
	param(j)=1.;
	x(i,j)=poly.eval(param,*ic);
	param(j)=0.;
      }
    }
    Real cv=calc_cv(x,y);
    //    debugfile << "pow=" << cur_pow << endl;
    if (cv<best_cv) {
      best_cv=cv;
      best_pow=cur_pow;
      best_x=x;
    }
    cur_pow++;
  }
  //  cerr << "pow=" << best_pow << endl;
  if (best_pow==-1) {
    ineqlist.detach_all();
    return 0.;
  }
  ineq.c=(Real)best_pow;
  product(&max_pow,ineq.v,ineq.c);
  MultiDimPoly poly(max_pow,ineqlist);

  Array<Real> beta,yhat,e,tmp;
  Array2d<Real> var;
  calc_ols(&beta,best_x,y);
  calc_ols_var(&var,best_x,y);
  Array<Real> xpow;
  poly.eval_powers(&xpow,cur_control);
  product(&tmp,var,xpow);
  Real varp=inner_product(xpow,tmp);

  product(&yhat,best_x,beta);
  diff(&e,y,yhat);
  //  cerr << "VAR=" << sqrt(varp) << " " << sqrt(inner_product(e,e)/(Real)(e.get_size()-beta.get_size())) << endl;
  Real sigma=sqrt(varp + inner_product(e,e)/(Real)(e.get_size()-beta.get_size()));

  Real cur_lro_hat=poly.eval(beta,cur_control);
  Real t=fabs(cur_lro_hat-linlro(cur_lro))/sigma;
  ineqlist.detach_all();
  return t;
}

int main(int argc, char *argv[]) {
  // set defaults for all parameters;
  char *delim="\t";
  int help=0;
  Real enclosed_radius=0.;
  char *controlfilename="control.in";
  int n_equil=-1;
  int n_step=-1;
  Real t_prec=0.;
  int which_col=1;
  Real mf_thresh=-1.;
  int mf_thresh_col=0;
  int init_gs=-2;
  Real phi0=MAXFLOAT;
  Real threshold=3.;
  int flip_span=-1;
  int rnd_walk=0;
  int sigdig=6;
  int quiet=0;
  char *outfile="mc.out";
  char *snapshotfile="mcsnapshot.out";
  char *snapshotnumfile="";
  Real kboltzman=1.;
  int keV=0;
  int seed=0;
  int droplast=0;
  int addmux=0;
  int half_shift=0;
  int include_last=0;
  char *conc_axes="";
  char *corrfunc_label="trigo";
  Real fdT=1e-2;
  Real fdmu=1e-2;
  char *my_init_str="";
  char *kspace_labels="";

  // parse command line;
  AskStruct options[]={
    {"","Multicomponent Eazy Monte Carlo Code " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-h","Help",BOOLVAL,&help},
    {"-er","Set the system  size so that a sphere of that radius must fit inside the simulation cell",REALVAL,&enclosed_radius},
    {"-cf","Control file specifying the ranges of temperature and chem. pot. scanned (default: control.in)",STRINGVAL,&controlfilename},
    {"-eq","Number of equilibration passes",INTVAL,&n_equil},
    {"-n","Number of averaging passes",INTVAL,&n_step},
    {"-tp","Target precision (optional, replaces -n and -eq)",REALVAL,&t_prec},
    {"-aq","Quantity that must meet the tolerance specified by -tp. "
           "1: energy (default), 2: long-range order, 3-: point correlations",INTVAL,&which_col},
    {"-gs","which ground state to use as initial config (-gs=-1 to use random state)",INTVAL,&init_gs},
    {"-k","Boltzman's constant (conversion factor from T to energy)",REALVAL,&kboltzman},
    {"-keV","Set Boltzman's constant to 8.617e-5 so that temperature is in K when energy is in eV",BOOLVAL,&keV},
    {"-g2c","Convert output to canonical rather than grand-canonical quantities",BOOLVAL,&addmux},
    {"-phi0","initial (grand) canonical potential (default: from mean field approx.)",REALVAL,&phi0},
    {"-tstat","Critical value of the test for discontinuity",REALVAL,&threshold},
    {"-mft","Mean field threshold (if |lte-mf|<mft, use mf values instead of mc)",REALVAL,&mf_thresh},
    {"-mftq","Quantity that must meet the tolerance specified by -mft."
           "0: phi (default) , 1: energy , 2: long-range order, 3-: concentrations",INTVAL,&mf_thresh_col},
    {"-sigdig","Number of significant digits printed",INTVAL,&sigdig},
    {"-q","Quiet (do not write to stdout)",BOOLVAL,&quiet},
    {"-o","Output file (default: mc.out)",STRINGVAL,&outfile},
    {"-oss","Output snapshot file (default: mcsnapshot.out)",STRINGVAL,&snapshotfile},
    {"-opss","Output periodic snapshot files (default: do not write)",STRINGVAL,&snapshotnumfile},
    {"-sd","Seed for random number generation (default: use clock)",INTVAL,&seed},
    {"-dl","Drop the last data point of each inner loop (after the phase transition occured)",BOOLVAL,&droplast},
    {"-is","File name containing a user-specified initial configuration (replaces -gs)",STRINGVAL,&my_init_str},
    {"-hf","Shift all coordinates by half a grid point in scan (e.g. 2 steps in [0,1] give 0.25,0.75 instead of 0,0.5)",BOOLVAL,&half_shift},
    {"-il","Include last coordinate in scan (e.g. 3 steps in [0,1] gives 0,0.5,1 instead of 0,0.333,0.666)",BOOLVAL,&include_last},
    {"-ts","Triangular scanning. Specify list of composition axes (e.g. -ts=1,2).",STRINGVAL,&conc_axes},
    {"-crf","Select correlation functions (default: trigo)",STRINGVAL,&corrfunc_label},
    {"-df","Max distance between flips (in unit cells) in grand-canonical mode",INTVAL,&flip_span},
    {"-rw","Use random walk algorithm (experimental)",BOOLVAL,&rnd_walk},
    {"-fdT","Temperature step for finite differences",REALVAL,&fdT},
    {"-fdmu","Chemical potential step for finite differences",REALVAL,&fdmu},
    {"-ks","Specify how k space ECI are calculated (e.g. -ks=cs).",STRINGVAL,&kspace_labels}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    return 1;
  }
  if (help) {
    cout << helpstring;
    exit(1);
  }

  // check and fix conflicts between options;
  if (strlen(my_init_str)>0) {
    init_gs=-1;
  }

  if (keV) {kboltzman=8.617e-5;}

  if (init_gs==-2)       ERRORQUIT("Specify -gs");
  if (! ((n_equil!=-1 && n_step!=-1) || t_prec>0.)) ERRORQUIT("Specify either (-eq and -n) or -tp");
  if (t_prec>0. && n_equil==-1) {n_equil=0;}
  if ( t_prec<0. ) ERRORQUIT("-tp value must be positive");

  rndseed(seed);
  int snapshotnum=0;

  // read the lattice file and do some error checking;
  Structure lattice;
  Array<Arrayint> labellookup;
  Array<AutoString> label;
  rMatrix3d axes;
  {
    ifstream latticefile("lat.in");
    if (!latticefile) ERRORQUIT("Unable to open lat.in");
    parse_lattice_file(&lattice.cell, &lattice.atom_pos, &lattice.atom_type, &labellookup, &label, latticefile, &axes);
    wrap_inside_cell(&lattice.atom_pos,lattice.atom_pos,lattice.cell);
  }

  // initialize a table of correlation functions;
  if (!check_plug_in(CorrFuncTable(),corrfunc_label)) {
    ERRORQUIT("Aborting");
  }
  CorrFuncTable *pcorrfunc=GenericPlugIn<CorrFuncTable>::create(corrfunc_label);
  pcorrfunc->init_from_site_type_list(labellookup);

  if (!quiet) {
    cerr << "Species:" << endl;
    for (int i=0; i<label.get_size(); i++) {
      cerr << label(i) << " ";
    }
    cerr << endl;
  }

  // Find spacegroup of the parent lattice;
  SpaceGroup spacegroup;
  spacegroup.cell=lattice.cell;
  find_spacegroup(&spacegroup.point_op,&spacegroup.trans,lattice.cell,lattice.atom_pos,lattice.atom_type);
  if (contains_pure_translations(spacegroup.point_op,spacegroup.trans)) {
    cerr << "Warning: unit cell is not primitive." << endl;
  }

  // setup matrices for conversion between correlations and concentrations;
  Array2d<Real> concmat;
  Array2d<Real> corr_to_conc;
  Array2d<Real> conc_to_fullconc;
  Array2d<Real> corr_to_fullconc;
  Array<Real> conc_to_fullconc_c;
  calc_corr_to_conc(&concmat,lattice,labellookup,spacegroup,*pcorrfunc,1);
  extract_nonredundant(&corr_to_conc,&conc_to_fullconc,&conc_to_fullconc_c, concmat);
  product(&corr_to_fullconc,conc_to_fullconc,corr_to_conc);

  // check for sites with multiple possible species;
  int nbsite=0;
  for (nbsite=0; nbsite<lattice.atom_pos.get_size(); nbsite++) {
    if (labellookup(lattice.atom_type(nbsite)).get_size()==1) break;
  }
  if (nbsite==0) ERRORQUIT("Need at least one site with multiple species in lattice file.");
  // extract only sites with multiple possible species since other sites useless for MC;
  Structure lattice_only;
  lattice_only.cell=lattice.cell;
  lattice_only.atom_pos.resize(nbsite);
  lattice_only.atom_type.resize(nbsite);
  for (int i=0; i<lattice_only.atom_pos.get_size(); i++) {
    lattice_only.atom_pos(i) =lattice.atom_pos(i);
    lattice_only.atom_type(i)=lattice.atom_type(i);
  }

  // read clusters;
  LinkedList<MultiCluster> clusterlist;
  {
    ifstream clusterfile("clusters.out");
    if (!clusterfile) ERRORQUIT("Unable to open clusters.out");
    Real maxlen=read_clusters_and_eci(&clusterlist,NULL,clusterfile,clusterfile,axes);
    if (enclosed_radius==0) enclosed_radius=maxlen;
  }

  // read the (possibly temperature-dependent) ECI;
  PolyInterpolatorBig<Array<Real> > teci;
  PolyInterpolatorBig<Array<Real> > teeci;
  read_t_dep_eci_file(&teci,&teeci,clusterlist.get_size(),kboltzman);

  // Write the file describing the columns of the output file;
  {
    int c=0;
    ofstream file("mcheader.out");
    c++; file << c << ":T" << endl;
    for (int i=0; i<label.get_size(); i++) {
      c++; file << c << ":mu(" << label(i) << ")" << endl;
    }
    char *methodname[]={"lte","mf","mc"};
    for (int m=0; m<3; m++) {
      c++; file << c << (addmux ? ":F_" : ":phi_") << methodname[m] << endl;
      c++; file << c << (addmux ? ":E_" : ":Egc_") << methodname[m] << endl;
      c++; file << c << ":lro_" << methodname[m] << endl;
      for (int i=0; i<label.get_size(); i++) {
	c++; file << c << ":x_" << methodname[m] << "(" << label(i) << ")" << endl;
      }
    }
    c++; file << c << ":use_mf" << endl;
    LinkedListIterator<MultiCluster> ic(clusterlist);
    for ( ; ic; ic++) {
      c++; file << c << ":corr_" << ic->clus.get_size() << "_" << get_cluster_length(ic->clus) << endl;
    }
  }

  // read the ground state structures;
  Array<Structure> gs_list;
  {
    ifstream strfile("gs_str.out");
    if (!strfile) ERRORQUIT("Unable to open gs_str.out");
    LinkedList<Structure> gs_llist;
    while (!strfile.eof()) {
      Structure str;
      parse_structure_file(&str.cell,&str.atom_pos,&str.atom_type,label,strfile);
      if (!fix_atom_type(&str, lattice, labellookup,1)) ERRORQUIT("Error reading gs_str.out");
      skip_to_next_structure(strfile);
      wrap_inside_cell(&str.atom_pos,str.atom_pos,str.cell);
      gs_llist << new Structure(str);
    }
    LinkedList_to_Array(&gs_list,gs_llist);
  }
  int nb_gs=gs_list.get_size();

  // setup the initial configuration;
  if (init_gs>=nb_gs || init_gs<-1) ERRORQUIT("-gs option is out of range.");
  Structure disordered=lattice_only;
  for (int i=0; i<disordered.atom_type.get_size(); i++) {
    disordered.atom_type(i)=-1;
  }
  Structure *p_init_str,*my_p_init_str;
  // either disordered alloy or ground state config;
  p_init_str=(init_gs==-1 ? &disordered : &gs_list(init_gs));
  Structure mystr;
  if (strlen(my_init_str)>0) {
    // or user-specified;
    ifstream strfile(my_init_str);
    if (!strfile) ERRORQUIT("Unable to open initial structure file");
    parse_structure_file(&mystr.cell,&mystr.atom_pos,&mystr.atom_type,label,strfile);
    if (!fix_atom_type(&mystr, lattice, labellookup,1)) ERRORQUIT("Error reading initial structure file");
    wrap_inside_cell(&mystr.atom_pos,mystr.atom_pos,mystr.cell);
    my_p_init_str=&mystr;
  }
  else {
    my_p_init_str=p_init_str;
  }

  // figure out the right supercell;
  iVector3d simple_supercell;
  find_common_simple_supercell(&simple_supercell,lattice_only.cell,my_p_init_str->cell);
  rMatrix3d mat_simple_supercell;
  mat_simple_supercell.diag(to_real(simple_supercell));
  iVector3d fit_sphere;
  fit_sphere=find_sphere_bounding_box(lattice_only.cell*mat_simple_supercell,2.*enclosed_radius);
  for (int i=0; i<3; i++) {
    simple_supercell(i)*=fit_sphere(i);
  }
  if (!quiet) cout << "Supercell size: " << simple_supercell << endl;

  // create the Monte Carlo object;
  MultiMonteCarlo *pmc;
  MultiKSpaceECI multi_kspace_eci;
  // if plug-in for k-space method specified, use it;
  if (strlen(kspace_labels)>0) {
    if (!check_plug_in(KSpaceECI(),kspace_labels)) {
      ERRORQUIT("Aborting");
    }
    make_plug_in_list(&multi_kspace_eci, kspace_labels);
    LinkedListIterator<KSpaceECI> i(multi_kspace_eci);
    for (; i; i++) {
      i->init(lattice_only,labellookup,label, simple_supercell, *pcorrfunc);
    }
    pmc=new KSpaceMultiMonteCarlo(lattice_only,labellookup,simple_supercell,spacegroup,clusterlist,*pcorrfunc,&multi_kspace_eci);
  }
  else {
    pmc=new MultiMonteCarlo(lattice_only,labellookup,simple_supercell,spacegroup,clusterlist,*pcorrfunc);
  }

  // find all combinations of spin flips allowed by  constrains;
  int flipmode=1;
  if (file_exists("conccons.in")) {
    ifstream file("conccons.in");
    LinkedList<LinearInequality> eq_list;
    read_equalities(&eq_list,label,file);
    Array2d<Real> conc_cons;
    ineq_list_to_Array(&conc_cons,NULL,eq_list);
    Array2d<Real> corr_cons;
    product(&corr_cons,conc_cons,corr_to_fullconc);
    pmc->find_all_allowed_flips(corr_cons,iVector3d(flip_span,flip_span,flip_span));
    flipmode=0;
  }


  // create object providing thermo properties using mean-field approx;
  CalcMeanField *pmf;
  if (strlen(my_init_str)==0 && init_gs!=-1) {
    pmf=new CalcMeanField(lattice_only, labellookup, *my_p_init_str, spacegroup, clusterlist, *pcorrfunc);
  }
  else {
    pmf=NULL;
  }

  // find reference energy scale (to make some numerical parameters less system-dependent);
  Real Escale=0.;
  {
    Array<Real> mult;
    pmc->get_cluster_mult(&mult);
    Array<Real> eci;
    teci.interpol(&eci,0.);
    for (int i=0; i<eci.get_size(); i++) {
      Escale+=mult(i)*fabs(eci(i));
    }
  }

  ofstream mcfile(outfile);

  if (!rnd_walk) { // do standard grid scan;
    // first setup the grid;
    Array<Real> control_org;
    Array<Real> maxdiv;
    Array2d<Real> dcontrol;
    {
      LinkedList<Array<Real> > dcontrol_list;
      LinkedList<Real> maxdiv_list;
      ifstream controlfile(controlfilename);
      if (!controlfile) {
	cerr << "Error reading control file: " << controlfilename << endl;
	ERRORQUIT("Either create file or specify other file with -cf option.");
      }
      while (skip_delim(controlfile)) {
	Array<Real> cur_cont(label.get_size()+1);
	for (int i=0; i<cur_cont.get_size(); i++) {
	  controlfile >> cur_cont(i);
	}
	cur_cont(0)*=kboltzman;
	if (control_org.get_size()==0) {
	  control_org=cur_cont;
	}
	else {
	  Real curmaxdiv;
	  controlfile >> curmaxdiv;
	  maxdiv_list << new Real(fabs(curmaxdiv));
	  diff(&cur_cont,cur_cont,control_org);
	  product(&cur_cont,cur_cont,sgn(curmaxdiv)/MAX(fabs(curmaxdiv)-(include_last ? 1. : 0.),1.));
	  dcontrol_list << new Array<Real>(cur_cont);
	}
      }
      LinkedList_to_Array(&maxdiv,maxdiv_list);
      dcontrol.resize(control_org.get_size(),dcontrol_list.get_size());
      LinkedListIterator<Array<Real> > idc(dcontrol_list);
      for (int i=0; idc; i++,idc++) {
	for (int j=0; j<control_org.get_size(); j++) {
	  dcontrol(j,i)=(*idc)(j);
	}
      }
      if (half_shift) {
        Array<Real> half(dcontrol_list.get_size());
	for (int i=0; i<half.get_size(); i++) {half(i)=0.5;}
	Array<Real> control_org_shift;
	product(&control_org_shift,dcontrol,half);
	sum(&control_org,control_org,control_org_shift);
      }
    }
    
    Array<int> conc_axes_a;
    {
      LinkedList<int> conc_axes_list;
      istrstream line(conc_axes);
      while (!line.eof()) {
	int i=-1;
	line >> i;
	if (i==-1) break;
	conc_axes_list << new int(i-1);
	skip_delim(line,",; ");
      }
      LinkedList_to_Array(&conc_axes_a,conc_axes_list);
    }
    
    Array<LinkedList<Real> > lro_list(maxdiv.get_size());
    Array<Real> phi(maxdiv.get_size());
    
    MultiDimIterator<Array<Real> > curpt(maxdiv);
    while (curpt) {
      Real sumc=0.;
      for (int a=0; a<conc_axes_a.get_size(); a++) {
	sumc+=((Array<Real> &)curpt)(conc_axes_a(a))/maxdiv(conc_axes_a(a));
      }
      if (sumc > 1.+zero_tolerance) {
	curpt++;
	continue;
      }
      
      if (((Array<Real> &)curpt)(0)==0) {
	pmc->init_structure(*my_p_init_str);
      }
      
      Array<Real> cur_control;
      product(&cur_control,dcontrol,curpt);
      sum(&cur_control,cur_control,control_org);
      Real T=cur_control(0);
      Array<Real> full_mu;
      extract_elements(&full_mu,cur_control,1,cur_control.get_size());
      Array<Real> mu;
      product(&mu,full_mu,corr_to_fullconc);
      Real muxc=inner_product(full_mu,conc_to_fullconc_c);
      
      Array<Real> eci;
      teci.interpol(&eci,T);
      
      const int phi_offset=0;
      const int E_offset=1;
      const int lro_offset=2;
      const int x_offset=3;
      
      Array<Real> lte_val(3+conc_to_fullconc_c.get_size());
      zero_array(&lte_val);
      if (pmf) {
	pmf->set_eci(eci);
	pmf->set_LTE(1);
	lte_val(phi_offset)=pmf->set_get_phi(T,mu)-muxc;
	lte_val(lro_offset)=pmf->get_lro();
	lte_val(E_offset)=calc_E_from_phi(pmf,teci,T,mu,fdT*Escale)-muxc;
	Array<Real> x_lte,ptcorr_lte;
	pmf->set_eci(eci);
	calc_conc_from_phi(&ptcorr_lte, pmf,T,mu,fdmu*Escale);
	product(&x_lte,corr_to_fullconc,ptcorr_lte);
	sum(&x_lte,x_lte,conc_to_fullconc_c);
	extract_elements(&lte_val, x_offset,x_lte,0,x_lte.get_size());
	if (addmux) {
	  Real mux=inner_product(full_mu,x_lte);
	  lte_val(phi_offset)+=mux;
	  lte_val(E_offset)+=mux;
	}
      }
      
      Array<Real> mf_val(3+conc_to_fullconc_c.get_size());
      zero_array(&mf_val);
      if (pmf) {
	pmf->set_eci(eci);
	pmf->set_LTE(0);
	mf_val(phi_offset)=pmf->set_get_phi(T,mu)-muxc;
	mf_val(lro_offset)=pmf->get_lro();
	mf_val(E_offset)=calc_E_from_phi(pmf,teci,T,mu,fdT*Escale)-muxc;
	Array<Real> x_mf,ptcorr_mf;
	pmf->set_eci(eci);
	calc_conc_from_phi(&ptcorr_mf, pmf,T,mu,fdmu*Escale);
	product(&x_mf,corr_to_fullconc,ptcorr_mf);
	sum(&x_mf,x_mf,conc_to_fullconc_c);
	extract_elements(&mf_val, x_offset,x_mf,0,x_mf.get_size());
	if (addmux) {
	  Real mux=inner_product(full_mu,x_mf);
	  mf_val(phi_offset)+=mux;
	  mf_val(E_offset)+=mux;
	}
      }
      
      int do_quick=(fabs(lte_val(mf_thresh_col)-mf_val(mf_thresh_col))<mf_thresh);
      
      Array<Real> mc_val(3+conc_to_fullconc_c.get_size());
      pmc->set_eci(eci);
      pmc->set_T_mu(T,mu);
      
      GenericAccumulator *paccum;
      if (do_quick) {
	paccum=create_accum(0,0.,0);
      } else {
	if (((Array<Real> &)curpt)(0)==0) {
	  pmc->run(n_equil,flipmode);
	}
	paccum=create_accum(n_step,t_prec,which_col-1);
      }
      run_mc(paccum,pmc,flipmode);
      if (strlen(snapshotnumfile)>0) {
          char *pdot=strchr(snapshotnumfile,'.');
          if (pdot) {
            ostrstream num;
            num << snapshotnum << '\0';
            const char *pnum=num.str();
            strcpy(pdot-strlen(pnum),pnum);
            *pdot='.';
            ofstream file(snapshotnumfile);
            file.setf(ios::fixed);
            file.precision(sigdig);
            pmc->view(labellookup,label,file,axes);
            snapshotnum++;
         }
      }
      const Array<Real> &mcdata=paccum->get_mean();
      const Array<Real> &mcdata2=paccum->get_var();
      
      mc_val(E_offset)=mcdata(0)-muxc;
      mc_val(lro_offset)=1.-mcdata(1);
      Array<Real> corr,mult;
      extract_elements(&corr,mcdata,2+mu.get_size(),2+mu.get_size()+eci.get_size());
      pmc->get_cluster_mult(&mult);
      fix_energy(&(mc_val(E_offset)), mult,corr,teci,T);
      
      Array<Real> x_mc,ptcorr_mc;
      extract_elements(&ptcorr_mc,mcdata,2,2+mu.get_size());
      product(&x_mc,corr_to_fullconc,ptcorr_mc);
      sum(&x_mc,x_mc,conc_to_fullconc_c);
      extract_elements(&mc_val, x_offset,x_mc,0,x_mc.get_size());
      if (addmux) {
	Real mux=inner_product(full_mu,x_mc);
	mc_val(phi_offset)+=mux;
	mc_val(E_offset)+=mux;
      }
      
      delete paccum;
      
      int looplevel=0;
      while (looplevel<maxdiv.get_size()) {
	if (((Array<Real> &)curpt)(looplevel)!=0) break;
	looplevel++;
      }
      
      if (init_gs==-1) {
	lte_val(lro_offset)=0.;
	mf_val(lro_offset)=0.;
	mc_val(lro_offset)=0.;
      }
      
      if (looplevel==maxdiv.get_size()) {
	for (int ll=0; ll<phi.get_size(); ll++) {
	  phi(ll)=(phi0==MAXFLOAT ? mf_val(phi_offset) : phi0);
	}
      } else {
	Array<Real> dphi(x_mc.get_size()+1);
	dphi(0)=-(mc_val(E_offset)-phi(looplevel))/T;
	for (int i=0; i<x_mc.get_size(); i++) {
	  dphi(1+i)=-x_mc(i);
	}
	Array<Real> dcol;
	extract_column(&dcol,dcontrol,looplevel);
	Real new_phi=phi(looplevel)+inner_product(dphi,dcol);
	if (do_quick) {new_phi=mf_val(phi_offset);}
	for (int ll=looplevel; ll>=0; ll--) {phi(ll)=new_phi;}
      }
      mc_val(phi_offset)=phi(0);
      
      if (do_quick) {mc_val=mf_val;}
      
      for (int ll=looplevel-1; ll>=0; ll--) {
	lro_list(ll).delete_all();
      }
      int quitloop=0;
      if (looplevel<maxdiv.get_size()) {
	Array<Real> lro_array;
	int ll=looplevel;
	do {
	  LinkedList_to_Array(&lro_array,lro_list(ll));
	  if (lro_array.get_size()>=3) break;
	  ll++;
	} while (ll<lro_list.get_size());
	Real tt;
	quitloop=(threshold!=0 && (tt=detect_discontinuity(lro_array,mc_val(lro_offset)))>threshold);
	// cout << lro_array.get_size() << " " << ll << " " << tt << endl;
      }
      
      if (!(droplast && quitloop)) {
	ostrstream line;
	line.setf(ios::fixed);
	line.precision(sigdig);
	if (((Array<Real> &)curpt)(0)==0) { line << endl; }
	line << T/kboltzman << delim;
	print_vector(line,full_mu,delim);
	print_vector(line,lte_val,delim);
	print_vector(line,mf_val,delim);
	print_vector(line,mc_val,delim);
	line << do_quick << delim;
	print_vector(line,corr,delim);
	line << endl << '\0';
	mcfile << line.str() << flush;
	if (!quiet) {cout << line.str() << flush;}
      }
      if (quitloop) {
	curpt.bump_up(looplevel);
      }
      else {
	if (looplevel==maxdiv.get_size()) {looplevel--;}
	for (int ll=looplevel; ll>=0; ll--) {
	  lro_list(ll) << new Real(mc_val(lro_offset));
	}
	curpt++;
      }
    }
    {
      ofstream file(snapshotfile);
      file.setf(ios::fixed);
      file.precision(sigdig);
      pmc->view(labellookup,label,file,axes);
    }
  }
  else { //randomwalk algorithm;
    Array<Real> org_control(label.get_size()+1);
    Array<Real> dcontrol(label.get_size()+1);
    Array<Real> lim_control(label.get_size()+1);
    {
      ifstream controlfile(controlfilename);
      if (!controlfile) {
	cerr << "Error reading control file: " << controlfilename << endl;
	ERRORQUIT("Either create file or specify other file with -cf option.");
      }
      for (int i=0; i<org_control.get_size(); i++) {
	controlfile >> org_control(i);
      }
      org_control(0)*=kboltzman;
      for (int i=0; i<dcontrol.get_size(); i++) {
	controlfile >> dcontrol(i);
      }
      dcontrol(0)*=kboltzman;
      for (int i=0; i<lim_control.get_size(); i++) {
	controlfile >> lim_control(i);
      }
      lim_control(0)*=kboltzman;
    }

    if (!pmf) {
      ERRORQUIT("Cannot use -gs=-1 or -is options with random walk algorithm");
    }
    pmc->init_structure(*my_p_init_str);

    Array<Real> cur_control,old_control;
    old_control=org_control;
    cur_control=org_control;

    LinkedList<Array<Real> > kept_control;
    LinkedList<Real> kept_lro;

    while (1) {
      Real T=cur_control(0);
      Array<Real> full_mu;
      extract_elements(&full_mu,cur_control,1,cur_control.get_size());
      Array<Real> mu;
      product(&mu,full_mu,corr_to_fullconc);
      Real muxc=inner_product(full_mu,conc_to_fullconc_c);
      
      Array<Real> eci;
      teci.interpol(&eci,T);
      
      const int phi_offset=0;
      const int E_offset=1;
      const int lro_offset=2;
      const int x_offset=3;
      
      Array<Real> lte_val(3+conc_to_fullconc_c.get_size());
      zero_array(&lte_val);
      pmf->set_eci(eci);
      pmf->set_LTE(1);
      lte_val(phi_offset)=pmf->set_get_phi(T,mu)-muxc;
      lte_val(lro_offset)=pmf->get_lro();
      lte_val(E_offset)=calc_E_from_phi(pmf,teci,T,mu,fdT*Escale)-muxc;
      Array<Real> x_lte,ptcorr_lte;
      pmf->set_eci(eci);
      calc_conc_from_phi(&ptcorr_lte, pmf,T,mu,fdmu*Escale);
      product(&x_lte,corr_to_fullconc,ptcorr_lte);
      sum(&x_lte,x_lte,conc_to_fullconc_c);
      extract_elements(&lte_val, x_offset,x_lte,0,x_lte.get_size());
      if (addmux) {
	Real mux=inner_product(full_mu,x_lte);
	lte_val(phi_offset)+=mux;
	lte_val(E_offset)+=mux;
      }
      
      Array<Real> mf_val(3+conc_to_fullconc_c.get_size());
      zero_array(&mf_val);
      pmf->set_eci(eci);
      pmf->set_LTE(0);
      mf_val(phi_offset)=pmf->set_get_phi(T,mu)-muxc;
      mf_val(lro_offset)=pmf->get_lro();
      mf_val(E_offset)=calc_E_from_phi(pmf,teci,T,mu,fdT*Escale)-muxc;
      Array<Real> x_mf,ptcorr_mf;
      pmf->set_eci(eci);
      calc_conc_from_phi(&ptcorr_mf, pmf,T,mu,fdmu*Escale);
      product(&x_mf,corr_to_fullconc,ptcorr_mf);
      sum(&x_mf,x_mf,conc_to_fullconc_c);
      extract_elements(&mf_val, x_offset,x_mf,0,x_mf.get_size());
      if (addmux) {
	Real mux=inner_product(full_mu,x_mf);
	mf_val(phi_offset)+=mux;
	mf_val(E_offset)+=mux;
      }
      
      int do_quick=(fabs(lte_val(mf_thresh_col)-mf_val(mf_thresh_col))<mf_thresh);
      
      Array<Real> mc_val(3+conc_to_fullconc_c.get_size());
      pmc->set_eci(eci);
      pmc->set_T_mu(T,mu);
      
      GenericAccumulator *paccum;
      if (do_quick) {
	paccum=create_accum(0,0.,0);
      } else {
	paccum=create_accum(n_step,t_prec,which_col-1);
      }
      run_mc(paccum,pmc,flipmode);
      const Array<Real> &mcdata=paccum->get_mean();
      const Array<Real> &mcdata2=paccum->get_var();
      
      mc_val(E_offset)=mcdata(0)-muxc;
      mc_val(lro_offset)=1.-mcdata(1);
      Array<Real> corr,mult;
      extract_elements(&corr,mcdata,2+mu.get_size(),2+mu.get_size()+eci.get_size());
      pmc->get_cluster_mult(&mult);
      fix_energy(&(mc_val(E_offset)), mult,corr,teci,T);
      
      Array<Real> x_mc,ptcorr_mc;
      extract_elements(&ptcorr_mc,mcdata,2,2+mu.get_size());
      product(&x_mc,corr_to_fullconc,ptcorr_mc);
      sum(&x_mc,x_mc,conc_to_fullconc_c);
      extract_elements(&mc_val, x_offset,x_mc,0,x_mc.get_size());
      if (addmux) {
	Real mux=inner_product(full_mu,x_mc);
	mc_val(phi_offset)+=mux;
	mc_val(E_offset)+=mux;
      }
      
      delete paccum;
      
      mc_val(phi_offset)=0.;
      if (do_quick) {mc_val=lte_val;}

      int reject=0;
      for (int i=0; i<lim_control.get_size()-1; i++) {
	if (mc_val(x_offset+i)<lim_control(1+i)) {reject=1;}
      }
      Array<Real> Tmu(1+mu.get_size());
      Tmu(0)=T;
      for (int i=0; i<mu.get_size(); i++) {Tmu(i+1)=mu(i);}
      if (!reject) {
	reject=(check_in_phase(kept_control,kept_lro,Tmu,mc_val(lro_offset))>threshold);
	//	debugfile << reject << endl;
      }
      {
	ostrstream line;
	line.setf(ios::fixed);
	line.precision(sigdig);
	line << T/kboltzman << delim;
	print_vector(line,full_mu,delim);
	print_vector(line,lte_val,delim);
	print_vector(line,mf_val,delim);
	print_vector(line,mc_val,delim);
	line << do_quick << delim;
	print_vector(line,corr,delim);
	line << endl << '\0';
	if (!reject) {mcfile << line.str() << flush;}
	if (!quiet) {cout << (reject ? 'R' : 'A') << delim << line.str() << flush;}
      }
      if (reject) {
	if (!do_quick) {pmc->init_structure(*my_p_init_str);}
      }
      else {
	kept_control << new Array<Real>(Tmu);
	kept_lro << new Real(mc_val(lro_offset));
	old_control=cur_control;
      }
      do {
	cur_control(0)=old_control(0)+dcontrol(0)*(2.*uniform01()-1.);
      } while (cur_control(0)>lim_control(0) || cur_control(0)<0.);
      for (int i=1; i<old_control.get_size(); i++) {
	cur_control(i)=old_control(i)+dcontrol(i)*(2.*uniform01()-1.);
      }
    }
  }
  delete pcorrfunc;
  delete pmc;
  delete pmf;
  }
