#include <fstream.h>
#include <strstream.h>
#include "mclib.h"
#include "drawpd.h"
#include "parse.h"
#include "getvalue.h"
#include "lstsqr.h"
#include "version.h"
#include "teci.h"
#include "keci.h"
#include "plugin.h"

extern char *helpstring;

GenericPlugIn<KSpaceECI> *GenericPlugIn<KSpaceECI>::list=NULL;

int main(int argc, char *argv[]) {
  char *delim="\t";
  int help=0;
  Real mu_l[]={MAXFLOAT,MAXFLOAT};
  Real T_l[]={MAXFLOAT,MAXFLOAT};
  Real dmu=MAXFLOAT;
  Real dT=MAXFLOAT;
  Real db=MAXFLOAT;
  int do_abs=0;
  Real enclosed_radius=0.;
  Real phi0=MAXFLOAT;
  int n_equil=-1;
  int n_step=-1;
  Real x_prec=0.;
  int which_col=1;
  int can_mode=0;
  Real init_conc=MAXFLOAT;
  int init_gs=-2;
  int innerT=0;
  Real threshold=3.;
  int sigdig=6;
  int quiet=0;
  int dodb=0;
  char *outfile="mc.out";
  char *snapshotfile="mcsnapshot.out";
  Real kboltzman=1.;
  int keV=0;
  int seed=0;
  int droplast=0;
  int addmux=0;
  char *my_init_str="";
  char *kspace_labels="";
  AskStruct options[]={
    {"","Eazy Monte Carlo Code " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-h","Help",BOOLVAL,&help},
    {"-mu0","initial chemical potential",REALVAL,&mu_l[0]},
    {"-T0","initial temperature",REALVAL,&T_l[0]},
    {"-mu1","final chemical potential",REALVAL,&mu_l[1]},
    {"-T1","final temperature",REALVAL,&T_l[1]},
    {"-dmu","chemical potential step",REALVAL,&dmu},
    {"-dT","temperature step",REALVAL,&dT},
    {"-db","inverse temperature step",REALVAL,&db},
    {"-cm","Set Canonical mode",BOOLVAL,&can_mode},
    {"-x","Set concentration (implies -cm)",REALVAL,&init_conc},
    {"-abs","take chemical potentials as absolute quantities (as in mc.out)",BOOLVAL,&do_abs},
    {"-phi0","initial (grand) canonical potential",REALVAL,&phi0},
    {"-er","set the system  size so that a sphere of that radius must fit inside the simulation cell",REALVAL,&enclosed_radius},
    {"-eq","number of equilibration passes",INTVAL,&n_equil},
    {"-n","number of averaging passes",INTVAL,&n_step},
    {"-dx","Target precision for the average concentration (optional, replaces -n and -eq)",REALVAL,&x_prec},
    {"-aq","Alternative quantity that must meet the tolerance specified by -dx. "
           "0: energy, 1: concentration (default), 2: long-range order, 3- correlations",INTVAL,&which_col},
    {"-gs","which ground state to use as initial config (-gs=-1 to use random state, c=1/2)",INTVAL,&init_gs},
    {"-innerT","inner loop over T",BOOLVAL,&innerT},
    {"-tstat","Critical value of the test for discontinuity",REALVAL,&threshold},
    {"-sigdig","Number of significant digits printed",INTVAL,&sigdig},
    {"-q","Quiet (do not write to stdout)",BOOLVAL,&quiet},
    {"-o","Output file (default: mc.out)",STRINGVAL,&outfile},
    {"-oss","Output snapshot file (default: mcsnapshot.out)",STRINGVAL,&snapshotfile},
    {"-k","Boltzman's constant (conversion factor from T to energy)",REALVAL,&kboltzman},
    {"-keV","Set Boltzman's constant to 8.617e-5 so that temperature is in K when energy is in eV",BOOLVAL,&keV},
    {"-sd","Seed for random number generation (default: use clock)",INTVAL,&seed},
    {"-dl","Drop the last data point of each inner loop (after the phase transition occured)",BOOLVAL,&droplast},
    {"-g2c","Convert output to canonical rather than grand-canonical quantities",BOOLVAL,&addmux},
    {"-is","File name containing a user-specified initial configuration (replaces -gs)",STRINGVAL,&my_init_str},
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
  if (init_conc!=MAXFLOAT) {
    if (!is_between(init_conc,-1.,1.)) ERRORQUIT("Initial composition out of ]-1,1[ range.");
    can_mode=1;
  }
  if (can_mode) {
    mu_l[0]=0.;
    mu_l[1]=0.;
    do_abs=1;
    cerr << "Running in canonical mode: limited features available" << endl;
  }
  if (strlen(my_init_str)>0) {
    init_gs=-1;
  }

  if (keV) {kboltzman=8.617e-5;}
  if (T_l[0]!=MAXFLOAT && T_l[1]==MAXFLOAT && mu_l[0]!=MAXFLOAT && mu_l[1]==MAXFLOAT) {mu_l[1]=mu_l[0]; T_l[1]=T_l[0];}
  if (T_l[0]!=MAXFLOAT && T_l[1]!=MAXFLOAT && mu_l[0]!=MAXFLOAT && mu_l[1]==MAXFLOAT) mu_l[1]=mu_l[0];
  if (mu_l[0]!=MAXFLOAT && mu_l[1]!=MAXFLOAT && T_l[0]!=MAXFLOAT && T_l[1]==MAXFLOAT) T_l[1]=T_l[0];
  if (mu_l[0]==mu_l[1]) {dmu=1.; innerT=1;}
  if (T_l[0]==T_l[1]) {dT=1.; db=MAXFLOAT;}
  if (mu_l[0]==MAXFLOAT) ERRORQUIT("Specify -mu0");
  if (mu_l[1]==MAXFLOAT) ERRORQUIT("Specify -mu1");
  if ( T_l[0]==MAXFLOAT) ERRORQUIT("Specify -T0");
  if ( T_l[1]==MAXFLOAT) ERRORQUIT("Specify -T1");
  if (    dmu==MAXFLOAT) ERRORQUIT("Specify -dmu");
  if (! ((n_equil!=-1 && n_step!=-1) || (x_prec>0.))) ERRORQUIT("Specify either (-eq and -n) or -dx");
  if ( ((n_equil!=-1 || n_step!=-1) && (x_prec>0.))) ERRORQUIT("Specify either (-eq and -n) or -dx");
  if (init_gs==-2)       ERRORQUIT("Specify -gs (use -gs=-1 for disordered state)");
  if ( ! ((dT==MAXFLOAT) ^ (db==MAXFLOAT)) ) ERRORQUIT("Specify either -dT or -db");
  if ( x_prec<0 ) ERRORQUIT("-dx value must be positive");

  T_l[0]*=kboltzman;
  T_l[1]*=kboltzman;
  if (dT!=MAXFLOAT) {
    dT*=kboltzman;
  }
  else {
    db/=kboltzman;
    dodb=1;
  }

  rndseed(seed);

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
  for (int i=0; i<labellookup.get_size(); i++) {
    if (labellookup(i).get_size()>2)
      ERRORQUIT ("For multicomponent systems, use memc2.");
  }
  SpaceGroup spacegroup;
  spacegroup.cell=lattice.cell;
  find_spacegroup(&spacegroup.point_op,&spacegroup.trans,lattice.cell,lattice.atom_pos,lattice.atom_type);
  if (contains_pure_translations(spacegroup.point_op,spacegroup.trans)) {
    cerr << "Warning: unit cell is not primitive." << endl;
  }

  int nbsite=0;
  for (nbsite=0; nbsite<lattice.atom_pos.get_size(); nbsite++) {
    if (labellookup(lattice.atom_type(nbsite)).get_size()==1) break;
  }
  if (nbsite==0) ERRORQUIT("Need at least one site with multiple species in lattice file.");
  Structure lattice_only;
  lattice_only.cell=lattice.cell;
  lattice_only.atom_pos.resize(nbsite);
  lattice_only.atom_type.resize(nbsite);
  for (int i=0; i<lattice_only.atom_pos.get_size(); i++) {
    lattice_only.atom_pos(i) =lattice.atom_pos(i);
    lattice_only.atom_type(i)=lattice.atom_type(i);
  }

  LinkedList<Cluster> clusterlist;
  {
    ifstream clusterfile("clusters.out");
    if (!clusterfile) ERRORQUIT("Unable to open clusters.out");
    Real maxlen=read_clusters_and_eci(&clusterlist,NULL,clusterfile,clusterfile,axes);
    if (enclosed_radius==0) enclosed_radius=maxlen;
  }

  PolyInterpolatorBig<Array<Real> > teci;
  PolyInterpolatorBig<Array<Real> > teeci;
  read_t_dep_eci_file(&teci,&teeci,clusterlist.get_size(),kboltzman);

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
      for (int s=0; s<str.atom_type.get_size(); s++) {
        str.atom_type(s)=-1+2*str.atom_type(s);
      }
      gs_llist << new Structure(str);
    }
    LinkedList_to_Array(&gs_list,gs_llist);
  }
  int nb_gs=gs_list.get_size();
  if (init_gs>=nb_gs || init_gs<-1) ERRORQUIT("-gs option is out of range.");

  Array<Real> gs_concentration(nb_gs),gs_energy(nb_gs);

  for (int gs=0; gs<nb_gs; gs++) {
    LinkedList<Real> ecilist;
    make_eci_list(&ecilist,teci,0.);
    LTE_ThermoData gs_lte(lattice_only,spacegroup,gs_list(gs),clusterlist,ecilist);
    gs_concentration(gs)=gs_lte.x(0,0);
    gs_energy(gs)=gs_lte.phi(0,0);
  }
  for (int gs=0; gs<nb_gs-1; gs++) {
    if (gs_concentration(gs)>gs_concentration(gs+1)) ERRORQUIT("Ground states must be in increasing order of concentration in gs_str.out.");
  }
  Array<Real> mu_flat(nb_gs-1);
  for (int gs=0; gs<nb_gs-1; gs++) {
    mu_flat(gs)=(gs_energy(gs+1)-gs_energy(gs))/(gs_concentration(gs+1)-gs_concentration(gs));
  }
  Structure disordered=lattice_only;
  for (int i=0; i<disordered.atom_type.get_size(); i++) {
    disordered.atom_type(i)=0;
  }
  if (!do_abs) {
    if (nb_gs==2) {
      mu_l[0]+=mu_flat(0);
      mu_l[1]+=mu_flat(0);
    }
    else if (nb_gs>2) {
      mu_l[0]-=1.;
      mu_l[1]-=1.;
      for (int s=0; s<2; s++) {
	if (init_gs!=-1 && fabs(mu_l[s]-((Real)init_gs-0.5))>0.5 && mu_l[s]>-0.5 && mu_l[s]<(Real)nb_gs-1.5) {
	  cerr << "Warning: Chemical potential " << mu_l[s]+1. << " may not stabilize ground state " << init_gs << endl;
	}
	int gs=(int)floor(mu_l[s]);
	if (gs<0) gs=0;
	if (gs>nb_gs-3) gs=nb_gs-3;
	Real frac_gs=mu_l[s]-(Real)gs;
	mu_l[s]=(1-frac_gs)*mu_flat(gs)+frac_gs*mu_flat(gs+1);
      }
      dmu=dmu*(mu_flat(nb_gs-2)-mu_flat(0))/(nb_gs-2);
    }
  }
  if (dodb) {
    if (!is_between(1./(1./T_l[0]+db),T_l[0],T_l[1])) {db=-db;}
  }
  else {
    if (!is_between(T_l[0]+dT,T_l[0],T_l[1])) {dT=-dT;}
  }
  if (!is_between(mu_l[0]+dmu,mu_l[0],mu_l[1])) {dmu=-dmu;}

  Structure *p_init_str,*my_p_init_str;
  p_init_str=(init_gs==-1 ? &disordered : &gs_list(init_gs));
  Structure mystr;
  if (strlen(my_init_str)>0) {
    ifstream strfile(my_init_str);
    if (!strfile) ERRORQUIT("Unable to open initial structure file");
    parse_structure_file(&mystr.cell,&mystr.atom_pos,&mystr.atom_type,label,strfile);
    if (!fix_atom_type(&mystr, lattice, labellookup,1)) ERRORQUIT("Error reading initial structure file");
    wrap_inside_cell(&mystr.atom_pos,mystr.atom_pos,mystr.cell);
    for (int s=0; s<mystr.atom_type.get_size(); s++) {
      mystr.atom_type(s)=-1+2*mystr.atom_type(s);
    }
    my_p_init_str=&mystr;
  }
  else {
    my_p_init_str=p_init_str;
  }

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
  MonteCarlo *pmc;
  MultiKSpaceECI multi_kspace_eci;
  if (strlen(kspace_labels)>0) {
    if (!check_plug_in(KSpaceECI(),kspace_labels)) {
      ERRORQUIT("Aborting");
    }
    make_plug_in_list(&multi_kspace_eci, kspace_labels);
    multi_kspace_eci.check_only_one();
    LinkedListIterator<KSpaceECI> i(multi_kspace_eci);
    for (; i; i++) {
      i->static_init(lattice_only);
    }
    pmc=new KSpaceMonteCarlo(lattice_only,simple_supercell,spacegroup,clusterlist,&multi_kspace_eci);
  }
  else {
    pmc=new MonteCarlo(lattice_only,simple_supercell,spacegroup,clusterlist);
  }

  ofstream mcfile(outfile);
  mcfile.setf(ios::fixed);
  mcfile.precision(sigdig);
  cout.setf(ios::fixed);
  cout.precision(sigdig);

  if (can_mode) {
    Real T=T_l[0];
    Real phi=(phi0==MAXFLOAT ? (init_gs==-1 ? 0 : gs_energy(init_gs)) : phi0);
    Real old_E;
    int new_scan=1;
    LinkedList<Real> order_list;

    pmc->init_structure(*my_p_init_str);
    if (init_conc!=MAXFLOAT) pmc->set_concentration(init_conc);

    while (1) {
      Array<Real> eci;
      teci.interpol(&eci,T);
      pmc->set_eci(eci);
      pmc->init_run(T,0.);
      if (new_scan) {
	if (x_prec==0.) pmc->run(n_equil,2);
      }
      
      MCOutputData mcdata;
      run_mc(&mcdata,pmc,2,n_step,x_prec,which_col);
      if (init_gs==-1) mcdata.lro=0.;
      
      //fix_energy(&mcdata,teci,T);
      fix_energy(&(mcdata.E),mcdata.mult,mcdata.corr,teci,T);

      Array<Real> order_array;
      LinkedList_to_Array(&order_array,order_list);
      int discont=(threshold==0 ? 0 : detect_discontinuity(order_array,mcdata.lro)>threshold);
      order_list << new Real(mcdata.lro);

      if (!new_scan) {
	if (dodb) {
	  Real b0=(1./T)-db;
	  Real b1=(1./T);
	  phi=(phi*b0+(mcdata.E+old_E)*db/2)/b1;
	}
	else {
	  Real b0=1./(T-dT);
	  Real b1=1./T;
	  phi=(phi*b0+(mcdata.E+old_E)*(b1-b0)/2)/b1;
	}
      }
      
      ostrstream line;
      line.setf(ios::fixed);
      line.precision(sigdig);

      if (!discont || !droplast) {
	line << T/kboltzman << delim
	     << 0 << delim
	     << mcdata.E << delim
	     << mcdata.x << delim
	     << phi << delim
	     << mcdata.heatcap << delim
	     << 0 << delim
	     << 0 << delim
	     << 0 << delim
	     << 0 << delim
	     << 0 << delim
	     << 0 << delim
	     << 0 << delim
	     << 0 << delim
	     << 0 << delim
	     << 0 << delim
	     << mcdata.lro << delim;
	for (int i=0; i<mcdata.corr.get_size(); i++) {line << mcdata.corr(i) << delim;}
	line << endl << '\0';
	mcfile << line.str() << flush;
	if (!quiet) {cout << line.str() << flush;}
      }

      old_E=mcdata.E;
      if (dodb) {
	T=1./((1./T)+db);
      }
      else {
	T+=dT;
      }
      new_scan=0;
      if (discont || !is_between(T,T_l[0],T_l[1])) break;
    }
  }
  else {
    TDependentECIPhase lte("LTE",lattice_only,spacegroup,*p_init_str,clusterlist,teci);
    TDependentECIPhase  mf("MF" ,lattice_only,spacegroup,*p_init_str,clusterlist,teci);
    TDependentECIPhase hte("HTE",lattice_only,spacegroup,*p_init_str,clusterlist,teci);
    
    Real T=T_l[0];
    Real mu=mu_l[0];
    Real phi=(phi0==MAXFLOAT ? (init_gs==-1 ? hte.phi(T,mu): lte.phi(T,mu)): phi0);
    int very_new_scan=1;
    int new_scan=1;
    Real old_E,old_x,very_old_E,very_old_x,very_old_phi;
    LinkedList<Real> order_list;
    LinkedList<Real> outer_order_list;
    while (1) {
      if (new_scan) {
	pmc->init_structure(*my_p_init_str);
	if (init_conc!=MAXFLOAT) pmc->set_concentration(init_conc);
	pmc->init_run(T,mu);
	order_list.delete_all();
      }
      
      Array<Real> eci;
      teci.interpol(&eci,T);
      pmc->set_eci(eci);
      pmc->init_run(T,mu);

      if (new_scan) {
	if (x_prec==0.) pmc->run(n_equil,1);
      }
      
      MCOutputData mcdata;
      run_mc(&mcdata,pmc,1,n_step,x_prec,which_col);
      if (init_gs==-1) mcdata.lro=0.;
      
      //fix_energy(&mcdata,teci,T);
      fix_energy(&(mcdata.E),mcdata.mult,mcdata.corr,teci,T);
      
      Array<Real> order_array;
      LinkedList_to_Array(&order_array,order_list);
      int discont=(threshold==0 ? 0 : detect_discontinuity(order_array,mcdata.lro)>threshold);
      order_list << new Real(mcdata.lro);
      
      if (new_scan) {
	if (!very_new_scan) {
	  if (innerT) {
	    phi=very_old_phi-(very_old_x+mcdata.x)/2*dmu;
	  }
	  else {
	    if (dodb) {
	      Real b0=(1./T)-db;
	      Real b1=(1./T);
	      phi=(very_old_phi*b0+(mcdata.E+very_old_E)*db/2)/b1;
	    }
	    else {
	      Real b0=1./(T-dT);
	      Real b1=1./T;
	      phi=(very_old_phi*b0+(mcdata.E+very_old_E)*(b1-b0)/2)/b1;
	    }
	  }
	}
	very_new_scan=0;
	very_old_E=mcdata.E;
	very_old_x=mcdata.x;
	very_old_phi=phi;
	LinkedList_to_Array(&order_array,outer_order_list);
	if (threshold!=0 && detect_discontinuity(order_array,mcdata.lro)>threshold) break;
	outer_order_list << new Real(mcdata.lro);
      }
      else { /* !new_scan */
	if (innerT) {
	  if (dodb) {
	    Real b0=(1./T)-db;
	    Real b1=(1./T);
	    phi=(phi*b0+(mcdata.E+old_E)*db/2)/b1;
	  }
	  else {
	    Real b0=1./(T-dT);
	    Real b1=1./T;
	    phi=(phi*b0+(mcdata.E+old_E)*(b1-b0)/2)/b1;
	  }
	}
	else {
	  phi+=-(old_x+mcdata.x)/2*dmu;
	}
      }
      
      Real E_lte=E_from_phi(T,mu,&lte);
      Real  E_mf=E_from_phi(T,mu,&mf);
      Real E_hte=E_from_phi(T,mu,&hte);
      
      ostrstream line;
      line.setf(ios::fixed);
      line.precision(sigdig);
      if (!discont || !droplast) {
	Real lte_x=lte.x(T,mu);
	Real mf_x=mf.x(T,mu);
	Real hte_x=hte.x(T,mu);
	{
	  ofstream ltefile("ltedat.out");
	  ofstream mffile("mfdat.out");
	  ofstream htefile("htedat.out");
	  lte.write(ltefile);
	  mf.write(mffile);
	  hte.write(htefile);
	}
	line << T/kboltzman << delim
	     << mu << delim
	     << mcdata.E+(addmux ? mcdata.x*mu : 0.) << delim
	     << mcdata.x << delim
	     << phi+(addmux ? mcdata.x*mu : 0.) << delim
	     << mcdata.heatcap << delim
	     << mcdata.suscept << delim
	     << E_lte+(addmux ? lte_x*mu : 0.) << delim
	     << lte_x << delim
	     << lte.phi(T,mu)+(addmux ? lte_x*mu : 0.) << delim
	     << E_mf+(addmux ? mf_x*mu : 0.) << delim
	     << mf_x << delim
	     << mf.phi(T,mu)+(addmux ? mf_x*mu : 0.) << delim
	     << E_hte+(addmux ? hte_x*mu : 0.) << delim
	     << hte_x << delim
	     << hte.phi(T,mu)+(addmux ? hte_x*mu : 0.) << delim
	     << mcdata.lro << delim;
	for (int i=0; i<mcdata.corr.get_size(); i++) {line << mcdata.corr(i) << delim;}
	line << endl << '\0';
	mcfile << line.str() << flush;
	if (!quiet) {cout << line.str() << flush;}
      }
      
      if (innerT) {
	old_E=mcdata.E;
	if (dodb) {
	  T=1./((1./T)+db);
	}
	else {
	  T+=dT;
	}
	new_scan=0;
	if (discont || !is_between(T,T_l[0],T_l[1])) {
	  mcfile << endl << flush;
	  if (!quiet) {cout << endl << flush;}
	  T=T_l[0];
	  mu+=dmu;
	  if (!is_between(mu,mu_l[0],mu_l[1])) break;
	  new_scan=1;
	}
      }
      else {
	old_x=mcdata.x;
	mu+=dmu;
	new_scan=0;
	if (discont || !is_between(mu,mu_l[0],mu_l[1])) {
	  mcfile << endl << flush;
	  if (!quiet) {cout << endl << flush;}
	  mu=mu_l[0];
	  if (dodb) {
	    T=1./((1./T)+db);
	  }
	  else {
	    T+=dT;
	  }
	  if (!is_between(T,T_l[0],T_l[1])) break;
	  new_scan=1;
	}
      }
    }
  }
  {
    ofstream file(snapshotfile);
    file.setf(ios::fixed);
    file.precision(sigdig);
    pmc->view(labellookup,label,file,axes);
  }
  delete pmc;
}

