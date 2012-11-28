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
  Real mu0=MAXFLOAT;
  Real T0=MAXFLOAT;
  Real maxdmu=MAXFLOAT;
  Real mu_gap=0;
  Real maxdT=MAXFLOAT;
  Real enclosed_radius=0.;
  int init_gs[2]={-2,-2};
  char *lat_dir[]={".","."};
  Real threshold=3.;
  int sigdig=6;
  int quiet=0;
  char *outfile="mc.out";
  Real kboltzman=1.;
  int keV=0;
  int seed=0;
  Real lte_prec=0;
  int go_down=0;
  Real large_step=1e50;
  Real x_prec=0;
  char *kspace_labels="";
  AskStruct options[]={
    {"","PHase Boundary " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-h","Help",BOOLVAL,&help},
    {"-mu","initial chemical potential",REALVAL,&mu0},
    {"-T","initial temperature",REALVAL,&T0},
    {"-dmu","chemical potential adjustment step",REALVAL,&maxdmu},
    {"-dT","temperature step",REALVAL,&maxdT},
    {"-mug","Gap between the mu in phase 1 and mu in phase 2 (default: 0)",REALVAL,&mu_gap},
    {"-ltep","threshold free energy precision to use MC instead of LTE (in units of T) (default: always MC)",REALVAL,&lte_prec},
    {"-er","enclosed radius",REALVAL,&enclosed_radius},
    {"-gs1","ground state for phase #1",INTVAL,&(init_gs[0])},
    {"-gs2","ground state for phase #2",INTVAL,&(init_gs[1])},
    {"-d1","directory for phase #1 (default: current dir)",STRINGVAL,&(lat_dir[0])},
    {"-d2","directory for phase #2 (default: current dir)",STRINGVAL,&(lat_dir[1])},
    {"-tstat","Critical value of the test for discontinuity",REALVAL,&threshold},
    {"-smax","Maximum step (experimental feature)",REALVAL,&large_step},
    {"-sigdig","Number of significant digits printed",INTVAL,&sigdig},
    {"-q","Quiet (do not write to stdout)",BOOLVAL,&quiet},
    {"-o","Output file (default: mc.out)",STRINGVAL,&outfile},
    {"-k","Boltzman's constant (conversion factor from T to energy)",REALVAL,&kboltzman},
    {"-keV","Set Boltzman's constant to 8.617e-5 so that temperature is in K when energy is in eV",BOOLVAL,&keV},
    {"-sd","Seed for random number generation (default: use clock)",INTVAL,&seed},
    {"-dn","Go down in temperature",BOOLVAL,&go_down},
    {"-dx","Concentration Precision",REALVAL,&x_prec},
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
  if (keV) {kboltzman=8.617e-5;}
  if (mu0==MAXFLOAT ^ T0==MAXFLOAT)  ERRORQUIT("Specify both -mu and -T   or   neither of them");
  if ((T0==MAXFLOAT || T0==0) && lte_prec==0) ERRORQUIT("Specify the -lte option: e.g. -ltep=5e-3");
  if (maxdT==MAXFLOAT)   ERRORQUIT("Specify -dT");
  if (init_gs[0]==-2) ERRORQUIT("Specify -gs1");
  if (init_gs[1]==-2) ERRORQUIT("Specify -gs2");
  if (enclosed_radius==0) ERRORQUIT("Specify -er");
  if (! (x_prec>0.)) ERRORQUIT("Specify -dx");
  T0*=kboltzman;
  maxdT*=kboltzman;
  lte_prec*=kboltzman;

  rndseed(seed);

  Structure init_str[2];
  PhaseThermoData *lte[2];
  PhaseThermoData *mf[2];
  MonteCarlo *mc[2];
  MultiKSpaceECI multi_kspace_eci[2];
  PolyInterpolatorBig<Array<Real> > teci[2];
  PolyInterpolatorBig<Array<Real> > teeci[2];
  Real e_scale[2];

  for (int phase=0; phase<2; phase++) {
    char cur_dir[MAX_LINE_LEN];
    getcwd(cur_dir,MAX_LINE_LEN);
    chdir(lat_dir[phase]);

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
	ERRORQUIT ("Sorry, only binary alloys implemented");
    }
    SpaceGroup spacegroup;
    spacegroup.cell=lattice.cell;
    find_spacegroup(&spacegroup.point_op,&spacegroup.trans,lattice.cell,lattice.atom_pos,lattice.atom_type);

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
      read_clusters_and_eci(&clusterlist,NULL,clusterfile,clusterfile,axes);
      read_t_dep_eci_file(&(teci[phase]),&(teeci[phase]),clusterlist.get_size(),kboltzman);
    }
    
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
    if (init_gs[phase]>=nb_gs || init_gs[phase]<-1) ERRORQUIT("-gs option is out of range.");

    Structure disordered=lattice_only;
    for (int i=0; i<disordered.atom_type.get_size(); i++) {disordered.atom_type(i)=0;}

    init_str[phase]=(init_gs[phase]==-1 ? disordered : gs_list(init_gs[phase]));
    
    lte[phase]=new TDependentECIPhase("LTE",lattice_only,spacegroup,init_str[phase],clusterlist,teci[phase]);
    mf[phase] =new TDependentECIPhase("MF",lattice_only,spacegroup,init_str[phase],clusterlist,teci[phase]);

    {
      Array<Real> eci;
      teci[phase].interpol(&eci,0.);
      LinkedList<Real> eci_l;
      eci_l << eci;
      HTE_ThermoData hte(lattice_only,spacegroup,init_str[phase],clusterlist,eci_l);
      e_scale[phase]=hte.F_of_x(0.,0.)-(hte.F_of_x(0.,-1.)+hte.F_of_x(0.,1.))/2.;
    }
    iVector3d simple_supercell;
    find_common_simple_supercell(&simple_supercell,lattice_only.cell,init_str[phase].cell);
    rMatrix3d mat_simple_supercell;
    mat_simple_supercell.diag(to_real(simple_supercell));
    iVector3d fit_sphere;
    fit_sphere=find_sphere_bounding_box(lattice_only.cell*mat_simple_supercell,2.*enclosed_radius);
    for (int i=0; i<3; i++) {simple_supercell(i)*=fit_sphere(i);}
    cout << "Phase " << phase+1 << " size: " << simple_supercell << endl;
    if (strlen(kspace_labels)>0) {
      if (lte_prec!=0) ERRORQUIT("LTE not implemented with k-space cluster expansion.");
      if (!check_plug_in(KSpaceECI(),kspace_labels)) {
	ERRORQUIT("Aborting");
      }
      make_plug_in_list(&(multi_kspace_eci[phase]), kspace_labels);
      multi_kspace_eci[phase].check_only_one();
      LinkedListIterator<KSpaceECI> i(multi_kspace_eci[phase]);
      for (; i; i++) {
	i->static_init(lattice_only);
      }
      mc[phase]=new KSpaceMonteCarlo(lattice_only,simple_supercell,spacegroup,clusterlist,&(multi_kspace_eci[phase]));
    }
    else {
      mc[phase]=new MonteCarlo(lattice_only,simple_supercell,spacegroup,clusterlist);
    }
    chdir(cur_dir);
  }
  if (maxdmu==MAXFLOAT) {maxdmu=MAX(fabs(e_scale[0]),fabs(e_scale[1]))/100.;}

  ofstream mcfile(outfile);
  mcfile.setf(ios::fixed);
  mcfile.precision(sigdig);
  cout.setf(ios::fixed);
  cout.precision(sigdig);

  Real old_dT=MAXFLOAT;
  Real old_dmu=MAXFLOAT;

  int left0_right1=(lte[1]->x(0,0) > lte[0]->x(0,0));

  if (mu0==MAXFLOAT) {
    mu0=(lte[1]->phi(0,0)-lte[0]->phi(0,0))/(lte[1]->x(0,0)-lte[0]->x(0,0));
    T0=maxdT;

    ostrstream line;
    line.setf(ios::fixed);
    line.precision(sigdig);
    line << 0 << delim << mu0 << delim
	 << lte[0]->x(0,mu0) << delim << lte[1]->x(0,mu0) << delim
	 << lte[0]->phi(0,mu0) << delim << lte[1]->phi(0,mu0) << delim
	 << endl << '\0';
    mcfile << line.str() << flush;
    if (!quiet) {cout << line.str() << flush;}
    old_dT=maxdT;
    old_dmu=0;
  }

  Real T=T0;
  Real mu=mu0;
  LinkedList<Real> order_list[2];

  int mc_started[2]={0,0};

  while (1) {
    Real E[2],x[2];
    int discont[2];
    Array<Real> order_array[2];
    for (int phase=0; phase<2; phase++) {
      LinkedList_to_Array(&order_array[phase],order_list[phase]);
    }
    int secondpass=0;
    while (1) {
      for (int phase=0; phase<2; phase++) {
	if (fabs(lte[phase]->phi(T,mu)-mf[phase]->phi(T,mu)) < lte_prec && lte[phase]->phi(T,mu)!=MAXFLOAT) {
	  E[phase]=E_from_phi(T,mu,lte[phase]);
	  x[phase]=lte[phase]->x(T,mu);
	  discont[phase]=0;
	}
	else {
	  if (!mc_started[phase]) {
	    mc[phase]->init_structure(init_str[phase]);
	    mc_started[phase]=1;
	  }
	  Array<Real> eci;
	  teci[phase].interpol(&eci,T);
	  mc[phase]->set_eci(eci);
	  mc[phase]->init_run(T,mu+mu_gap*(phase==0 ? -1.:1. ));
	  MCOutputData mcdata;
	  run_mc(&mcdata,mc[phase],1,0,x_prec);
	  //	  fix_energy(&mcdata,teci[phase],T);
	  fix_energy(&(mcdata.E),mcdata.mult,mcdata.corr,teci[phase],T);
	  E[phase]=mcdata.E;
	  x[phase]=mcdata.x;
	  cout << "Phase " << phase+1 << " n_equil= " << mcdata.n_equil << " n_avg= " << mcdata.n_step << endl;
	  discont[phase]=(detect_discontinuity(order_array[phase],x[phase])
                       > detect_discontinuity(order_array[1-phase],x[phase]));
	}
      }
      if (!discont[0] && !discont[1]) {break;}
      cout << discont[0] << delim << discont[1] << endl;
      if (discont[0] && discont[1]) {
	mc_started[0]=0;
	mc_started[1]=0;
	continue;
      }
      cout << "Looking for phase transition..." << endl;
      cout << "mu           x" << endl;
      if (secondpass) maxdmu/=2.;
      secondpass=1;
      Real mu_bounds[2];
      Real new_mu=mu;
      Real dmu=maxdmu*(discont[1] ^ left0_right1 ? -1 : 1);
      int phase=discont[1];
      Real xx;
      do {
        new_mu+=dmu;
	mc[phase]->init_run(T,new_mu);
	MCOutputData mcdata;
	run_mc(&mcdata,mc[phase],1,0,x_prec);
	xx=mcdata.x;
	cout << new_mu << delim << xx << endl;
      } while (detect_discontinuity(order_array[phase],xx)>detect_discontinuity(order_array[1-phase],xx));
      mu_bounds[0]=new_mu;
      dmu=-dmu;
      do {
        new_mu+=dmu;
	mc[phase]->init_run(T,new_mu);
	MCOutputData mcdata;
	run_mc(&mcdata,mc[phase],1,0,x_prec);
	xx=mcdata.x;
	cout << new_mu << delim << xx << endl;
      } while (detect_discontinuity(order_array[1-phase],xx)>detect_discontinuity(order_array[phase],xx));
      mu_bounds[1]=new_mu;
      mu=(mu_bounds[0]+mu_bounds[1])/2.;
      mc_started[phase]=0;
      
      old_dT=MAXFLOAT;
      old_dmu=MAXFLOAT;
    }
    if (discont[0] && discont[1]) {break;}
    for (int phase=0; phase<2; phase++) {
      order_list[phase] << new Real(x[phase]);
    }
    
    ostrstream line;
    line << T/kboltzman << delim << mu << delim
         << x[0] << delim << x[1] << delim
         << E[0] << delim << E[1] << delim
         << endl << '\0';
    mcfile << line.str() << flush;
    if (!quiet) {cout << line.str() << flush;}

    Real dT,dmu;
    if (fabs(x[1]-x[0])<threshold*x_prec) {break;}
    if (fabs(E[1]-E[0])/sqr(T)<1./(large_step*maxdT)) {
      dT=maxdT;
      dmu=0.;
    }
    else {
      dT=sqr(T)/(E[1]-E[0]);
      dmu=-T/(x[1]-x[0]);
      Real scale=dT/maxdT;
      dT/=scale;
      dmu/=scale;
    }
    if (dT<0 ^ go_down) {dT=-dT; dmu=-dmu;}

    if (old_dT==MAXFLOAT) {
      T+=dT;
      mu+=dmu;
    }
    else {
      T +=1.5*dT -0.5*old_dT;
      mu+=1.5*dmu-0.5*old_dmu;
    }
    if (T<0) {break;}
    old_dT=dT;
    old_dmu=dmu;
  }
  for (int phase=0; phase<2; phase++) {
    delete lte[phase];
    delete mf[phase];
    delete mc[phase];
  }
}
