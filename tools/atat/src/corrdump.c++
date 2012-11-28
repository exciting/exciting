#include <fstream.h>
#include "parse.h"
#include "clus_str.h"
#include "getvalue.h"
#include "ctype.h"
#include "version.h"
#include "plugin.h"

#define MAXMULTIPLET 6

extern char *helpstring;

int main(int argc, char *argv[]) {
  char *delim="\t";
  int dohelp=0;
  Array<Real> maxd(MAXMULTIPLET+1);
  zero_array(&maxd);
  char *latticefilename="lat.in";
  char *strfilename="str.out";
  char *dilutefilename="";
  int doconc=0;
  int doconcmat=0;
  int dosym=0;
  int doclus=0;
  int readclusters=0;
  int sigdig=5;
  int noempty=0;
  int nopoint=0;
  int rndcorr=0;
  char *writeunrel="";
  char *ecifile="";
  int multincl=0;
  char *corrfunc_label="trigo";
  int printnum=0;
  AskStruct options[]={
    {"","CORRelation DUMPer " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-h","Display more help",BOOLVAL,&dohelp},
    {"-2","Maximum distance between two points within a pair",REALVAL,&maxd(2)},
    {"-3","Maximum distance between two points within a triplet",REALVAL,&maxd(3)},
    {"-4","Maximum distance between two points within a quadruplet",REALVAL,&maxd(4)},
    {"-5","Maximum distance between two points within a quintuplet",REALVAL,&maxd(5)},
    {"-6","Maximum distance between two points within a sextuplet",REALVAL,&maxd(6)},
    {"-l","Input file defining the lattice   (default: lat.in)",STRINGVAL,&latticefilename},
    {"-s","Input file defining the structure (default: str.out)",STRINGVAL,&strfilename},
    {"-pc","Print composition only",BOOLVAL,&doconc},
    {"-pcm","Print composition matrix only",BOOLVAL,&doconcmat},
    {"-sym","Just find space group",BOOLVAL,&dosym},
    {"-clus","Just find clusters",BOOLVAL,&doclus},
    {"-c","Read clusters.out file instead of writing it",BOOLVAL,&readclusters},
    {"-z","Tolerance for finding symmetry operations (default: 1e-3)",REALVAL,&zero_tolerance},
    {"-sig","Number of significant digits printed (default: 5)",INTVAL,&sigdig},
    {"-noe","Do not include empty cluster",BOOLVAL,&noempty},
    {"-nop","Do not include point cluster(s)",BOOLVAL,&nopoint},
    {"-wu","Write Unrelaxed structure into specified file",STRINGVAL,&writeunrel},
    {"-rnd","Print correlation of the random state of the same composition as the input structure",BOOLVAL,&rndcorr},
    {"-eci","Predict quantity using ECI in specified file",STRINGVAL,&ecifile},
    {"-mi","Multiplicities are already included in ECI file",BOOLVAL,&multincl},
//    {"-dil","Input file defining the dilute sites",STRINGVAL,&dilutefilename},
    {"-crf","Select correlation functions (default: trigo)",STRINGVAL,&corrfunc_label},
    {"-nb","Print structure number",BOOLVAL,&printnum}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    return 1;
  }
  if (dohelp) {
    cout << helpstring;
    return 1;
  }
  Structure lattice;
  Array<Arrayint> labellookup;
  Array<AutoString> label;
  rMatrix3d axes;
  {
    ifstream latticefile(latticefilename);
    if (!latticefile) ERRORQUIT("Unable to open lattice file");
    parse_lattice_file(&lattice.cell, &lattice.atom_pos, &lattice.atom_type, &labellookup, &label, latticefile, &axes);
    wrap_inside_cell(&lattice.atom_pos,lattice.atom_pos,lattice.cell);
  }
  SpaceGroup spacegroup;
  spacegroup.cell=lattice.cell;
  find_spacegroup(&spacegroup.point_op,&spacegroup.trans,lattice.cell,lattice.atom_pos,lattice.atom_type);
  if (contains_pure_translations(spacegroup.point_op,spacegroup.trans)) {
    cerr << "Warning: unit cell is not primitive." << endl;
  }
  {
    ofstream symfile("sym.out");
    symfile.setf(ios::fixed);
    symfile.precision(sigdig);
    symfile << spacegroup.point_op.get_size() << endl;
    for (int i=0; i<spacegroup.point_op.get_size(); i++) {
      symfile << ((!axes)*(spacegroup.point_op(i))*axes) << endl;
      symfile << ((!axes)*(spacegroup.trans(i))) << endl;
      symfile << endl;
    }
  }
  if (dosym) return 0;
  int nbsite=0;
  for (nbsite=0; nbsite<lattice.atom_pos.get_size(); nbsite++) {
    if (labellookup(lattice.atom_type(nbsite)).get_size()==1) break;
  }
  if (nbsite==0 && strlen(dilutefilename)==0) ERRORQUIT("Need at least one site with multiple species in lattice file.");
  Structure lattice_only;
  lattice_only.cell=lattice.cell;
  lattice_only.atom_pos.resize(nbsite);
  lattice_only.atom_type.resize(nbsite);
  for (int i=0; i<lattice_only.atom_pos.get_size(); i++) {
    lattice_only.atom_pos(i) =lattice.atom_pos(i);
    lattice_only.atom_type(i)=lattice.atom_type(i);
  }
  Structure lat_func=lattice_only;
  for (int at=0; at<lat_func.atom_type.get_size(); at++) {
    lat_func.atom_type(at)=labellookup(lat_func.atom_type(at)).get_size();
  }
  if (strlen(dilutefilename)>0) {
    LinkedList<Array<rVector3d> > dil_list;
    rMatrix3d icell=!lat_func.cell;
    Array<int> done(lat_func.atom_pos.get_size());
    zero_array(&done);
    ifstream dilfile(dilutefilename);
    if (!dilfile) ERRORQUIT("Unable to open dilute sites file.");
    while (!dilfile.eof()) {
      AutoString tmp;
      get_string(&tmp,dilfile,"\n");
      if (dilfile.eof()) break;
      ifstream line(tmp);
      rVector3d pos;
      line >> pos;
      for (int op=0; op<spacegroup.point_op.get_size(); op++) {
	rVector3d tpos=spacegroup.point_op(op)*pos+spacegroup.trans(op);
	int at=which_atom(lat_func.atom_pos,tpos,icell);
	if (!done(at)) {
	  lat_func.atom_type(at)++;
	  done(at)=1;
	}
      }
    }
  }

  LinkedList<MultiCluster> clusterlist;
  if (doconc || doconcmat) {
    ofstream file("atoms.out");
    for (int i=0; i<label.get_size(); i++) {
      file << label(i) << endl;
    }
  }
  if ((!doconc && !doconcmat) || strlen(ecifile)>0) {
    if (!readclusters) {
      int maxp=MAXMULTIPLET;
      while (maxd(maxp)==0 && maxp>0) maxp--;
      if (maxp==0) ERRORQUIT("Specify -2=range -3=range etc.");
      Array<MultiClusterBank *> pclusterbanks(maxp+1);
      Array<int> nbcluster(maxp+1);
      pclusterbanks(1)=new MultiClusterBank(lat_func,1,spacegroup);
      nbcluster(1)=pclusterbanks(1)->get_cluster_list().get_size();
      for (int m=2; m<=maxp; m++) {
	pclusterbanks(m)=new MultiClusterBank(lat_func,m,spacegroup);
	if (maxd(m)>0) {
	  while (get_length_quick((*pclusterbanks(m))->clus)<maxd(m)+zero_tolerance) (*pclusterbanks(m))++;
	}
	nbcluster(m)=pclusterbanks(m)->get_current_index();
      }
      if (! noempty) clusterlist << new MultiCluster(0);
      for (int m=(nopoint ? 2 : 1); m<=maxp; m++) {
	pclusterbanks(m)->reset();
	for (int c=0; c<nbcluster(m); c++) {
	  clusterlist << new MultiCluster(*pclusterbanks(m));
	  (*pclusterbanks(m))++;
	}
	delete pclusterbanks(m);
      }
      {
	ofstream clusterfile("clusters.out");
	clusterfile.setf(ios::fixed);
	clusterfile.precision(sigdig);
	LinkedListIterator<MultiCluster> icluster(clusterlist);
	for ( ; icluster; icluster++) {
	  int mult=calc_multiplicity(*icluster, spacegroup.cell, spacegroup.point_op, spacegroup.trans);
	  clusterfile << mult << endl;
	  clusterfile << get_length_quick(icluster->clus) << endl;
	  clusterfile << icluster->clus.get_size() << endl;
	  for (int i=0; i<icluster->clus.get_size(); i++) {
	    clusterfile << ((!axes)*(icluster->clus(i))) << " " << icluster->site_type(i) << " " << icluster->func(i) <<  endl;
	  }
	  clusterfile << endl;
	}
      }
    }
    else {
      ifstream clusterfile("clusters.out");
      if (!clusterfile) ERRORQUIT("Unable to open clusters.out file");
      read_clusters_and_eci(&clusterlist,NULL,clusterfile,clusterfile,axes);
    }
    if (doclus) return 0;
  }
  if (strlen(dilutefilename)>0) return 0;

  LinkedList<Array<MultiCluster> > eq_clusterlist;
  LinkedListIterator<MultiCluster> icluster(clusterlist);
  for ( ; icluster; icluster++) {
    Array<MultiCluster> *pmulticlus=new Array<MultiCluster>;
    find_equivalent_clusters(pmulticlus, *icluster, spacegroup.cell, spacegroup.point_op, spacegroup.trans);
    eq_clusterlist << pmulticlus;
  }
  
  // initialize a table of correlation functions;
  if (!check_plug_in(CorrFuncTable(),corrfunc_label)) {
    ERRORQUIT("Aborting");
  }
  CorrFuncTable *pcorrfunc=GenericPlugIn<CorrFuncTable>::create(corrfunc_label);
  // pcorrfunc->init(max(lat_func.atom_type));
  pcorrfunc->init_from_site_type_list(labellookup);

  if (doconcmat) {
    Array2d<Real> corr_to_fullconc; // converts correlation to redundant concentration;
    Array2d<Real> corr_to_conc; // converts correlation to non redundant concentration;
    Array2d<Real> conc_to_fullconc; // converts nonredundant concentration to full redundant concentration;
    Array<Real> conc_to_fullconc_c; // constant terms of previous matrix;
    calc_corr_to_conc(&corr_to_fullconc,lattice,labellookup,spacegroup,*pcorrfunc);
    extract_nonredundant(&corr_to_conc,&conc_to_fullconc,&conc_to_fullconc_c, corr_to_fullconc);
    cout.setf(ios::fixed);
    cout.precision(sigdig);
    cout << "corr_to_fullconc" << endl;
    cout << corr_to_fullconc << endl;
    cout << "corr_to_conc" << endl;
    cout << corr_to_conc << endl;
    cout << "conc_to_fullconc" << endl;
    cout << conc_to_fullconc << endl;
    cout << "conc_to_fullconc_c" << endl;
    cout << conc_to_fullconc_c << endl;
    return 0;
  }

  Array<Real> eci;
  if (strlen(ecifile)>0) {
    ifstream file(ecifile);
    if (!file) ERRORQUIT("Unable to open ECI file.");
    LinkedList<Real> ecilist;
    while (skip_delim(file)) {
      Real e;
      file >> e;
      ecilist << new Real(e);
    }
    LinkedList_to_Array(&eci,ecilist);
    if (clusterlist.get_size()!=eci.get_size()) ERRORQUIT("Number of ECI does not match number of clusters."); 
  }

  Structure str;
  ifstream strfile(strfilename);
  if (!strfile) ERRORQUIT("Unable to open structure file");
  int strnum=1;
  while (!strfile.eof()) {
    parse_structure_file(&str.cell,&str.atom_pos,&str.atom_type,label,strfile);
    strfile.clear();
    char c;
    while (!isdigit(c=strfile.get()) && !strfile.eof());
    strfile.putback(c);
    wrap_inside_cell(&str.atom_pos,str.atom_pos,str.cell);

    ofstream logfile("corrdump.log");
    logfile.setf(ios::fixed);
    logfile.precision(sigdig);
    rMatrix3d supercell=(!lattice.cell)*str.cell;
    rMatrix3d rounded_supercell=to_real(to_int(supercell));
    rMatrix3d transfo=lattice.cell*rounded_supercell*(!str.cell);
    logfile << "strain=" << endl << (FixedMatrix<Real,3> &)transfo << endl;
    str.cell=transfo*str.cell;
    for (int i=0; i<str.atom_pos.get_size(); i++) {
      str.atom_pos(i)=transfo*str.atom_pos(i);
    }

    Structure ideal_str;
    ideal_str.cell=str.cell;
    find_all_atom_in_supercell(&ideal_str.atom_pos,&ideal_str.atom_type, lattice_only.atom_pos,lattice_only.atom_type,lattice_only.cell, ideal_str.cell);

    Array<int> spin(ideal_str.atom_type.get_size());
    zero_array(&spin);

    for (int i=0; i<str.atom_pos.get_size(); i++) {
      Real min_d=MAXFLOAT;
      int best_j,best_k;
      for (int j=0; j<str.atom_pos.get_size(); j++) {
	for (int k=0; k<ideal_str.atom_pos.get_size(); k++) {
	  if (str.atom_type(j)!=-1 && ideal_str.atom_type(k)!=-1) {
	    if (is_in_array(labellookup(ideal_str.atom_type(k)),str.atom_type(j))) {
	      rVector3d v=ideal_str.atom_pos(k)-str.atom_pos(j);
	      v=ideal_str.cell*cylinder((!ideal_str.cell)*v);
	      Real d=norm(v);
	      if (d<min_d) {
		best_j=j;
		best_k=k;
		min_d=d;
	      }
	    }
	  }
	}
      }
      if (min_d!=MAXFLOAT) {
	logfile << label(str.atom_type(best_j)) << " atom at " << ((!axes)*(str.atom_pos(best_j))) << " mapped onto lattice site at " << ((!axes)*(ideal_str.atom_pos(best_k))) << " distance= " << min_d << endl;
	spin(best_k)=index_in_array(labellookup(ideal_str.atom_type(best_k)),str.atom_type(best_j));
	str.atom_type(best_j)=-1;
	ideal_str.atom_type(best_k)=-1;
      }
    }
    int vac;
    for (vac=0; vac<label.get_size(); vac++) {
      if (strcmp(label(vac),"Vac")==0) break;
    }
    for (int k=0; k<ideal_str.atom_pos.get_size(); k++) {
      if (ideal_str.atom_type(k)!=-1) {
	if (!is_in_array(labellookup(ideal_str.atom_type(k)),vac)) ERRORQUIT("Unexpected vacancy");
	spin(k)=index_in_array(labellookup(ideal_str.atom_type(k)),vac);
	logfile << "Vacancy at " << ((!axes)*(ideal_str.atom_pos(k))) << endl;
      }
    }
    ideal_str.atom_type=spin;
    for (int i=0; i<str.atom_type.get_size(); i++) {
      if (str.atom_type(i)!=-1) {
	logfile << label(str.atom_type(i)) << " atom ignored at " << ((!axes)*(str.atom_pos(i))) << endl;
      }
    }

    if (strlen(writeunrel)>0) {
      ofstream file(writeunrel);
      if (!file) {ERRORQUIT("Unable to write unrelaxed structure");}
      Structure full_ideal_str;
      full_ideal_str.cell=ideal_str.cell;
      find_all_atom_in_supercell(&full_ideal_str.atom_pos,&full_ideal_str.atom_type, lattice.atom_pos,lattice.atom_type,lattice.cell, full_ideal_str.cell);
      int at0=0;
      for (int at=0; at<full_ideal_str.atom_type.get_size(); at++) {
	if (labellookup(full_ideal_str.atom_type(at)).get_size()>1) {
	  if (norm(full_ideal_str.atom_pos(at)-ideal_str.atom_pos(at0))>zero_tolerance) {ERRORQUIT("Position mismatch in write unrelaxed feature");}
	  full_ideal_str.atom_type(at)=ideal_str.atom_type(at0);
	  at0++;
	}
	else {
	  full_ideal_str.atom_type(at)=0;
	}
      }
      file.setf(ios::fixed);
      file.precision(sigdig);
      write_structure(full_ideal_str, lattice, labellookup, label, axes, file);
    }
    
    cout.setf(ios::fixed);
    //    cout.setf(ios::showpos);
    cout.precision(sigdig);
    if (printnum) cout << strnum << delim;
    strnum++;
    Real pred=0.;
    if (doconc) {
      Array<Real> conc;
      Structure cstr(ideal_str);
      calc_concentration(&conc,lattice,labellookup,cstr);
      for (int i=0; i<conc.get_size(); i++) {
	cout << conc(i) << " ";
      }
    }
    if (rndcorr) {
      LinkedList<Real> pointcorr;
      LinkedListIterator<Array<MultiCluster> > ipcluster(eq_clusterlist);
      for ( ; ipcluster; ipcluster++) {
	if ((*ipcluster)(0).clus.get_size() > 1) break;
	if ((*ipcluster)(0).clus.get_size()==1) {
	  pointcorr << new Real(calc_correlation(ideal_str, *ipcluster, spacegroup.cell, *pcorrfunc));
	}
      }
      int ieci=0;
      LinkedListIterator<MultiCluster> icluster(clusterlist);
      LinkedListIterator<Array<MultiCluster> > ieqcluster(eq_clusterlist);
      for ( ; icluster; icluster++, ieqcluster++, ieci++) {
	Real rho=1.;
	for (int s=0; s<icluster->clus.get_size(); s++) {
	  LinkedListIterator<MultiCluster> ipcluster(clusterlist);
	  while (ipcluster->clus.get_size()==0) {ipcluster++;}
	  LinkedListIterator<Real> ipcorr(pointcorr);
	  while (1) {
	    if (ipcluster->site_type(0)==icluster->site_type(s) && ipcluster->func(0)==icluster->func(s)) {
	      if (equivalent_by_symmetry(ipcluster->clus(0),icluster->clus(s),spacegroup.cell,spacegroup.point_op,spacegroup.trans)) {
		break;
	      }
	    }
	    ipcluster++;
	    ipcorr++;
	  }
	  rho*=(*ipcorr);
	}
	if (strlen(ecifile)>0) {
	  pred+=eci(ieci)*(multincl ? 1 : ieqcluster->get_size())*rho;
	}
	else {
	  if (!doconc) {cout << rho << delim;}
	}
      }
    }
    else {
      int ieci=0;
      LinkedListIterator<Array<MultiCluster> > icluster(eq_clusterlist);
      for ( ; icluster; icluster++, ieci++) {
	Real rho=calc_correlation(ideal_str, *icluster, spacegroup.cell, *pcorrfunc);
        if (strlen(ecifile)>0) {
	  pred+=eci(ieci)*(multincl ? 1 : icluster->get_size())*rho;
	}
	else {
	  if (!doconc) {cout << rho << delim;}
	}
      }
    }
    if (strlen(ecifile)>0) {
      cout << pred;
    }
    cout << endl;
  }
}

