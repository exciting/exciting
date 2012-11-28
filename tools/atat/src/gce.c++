#include <fstream.h>
#include "parse.h"
#include "clus_str.h"
#include "getvalue.h"
#include "ctype.h"
#include "version.h"
#include "plugin.h"
#include "gceutil.h"

#define MAXMULTIPLET 6

extern char *helpstring;

int main(int argc, char *argv[]) {
  char *delim="\t";
  int dohelp=0;
  Array<Real> maxd(MAXMULTIPLET+1);
  zero_array(&maxd);
  char *latticefilename="lat.in";
  char *strfilename="str.out";
  int doconc=0;
  int dosym=0;
  int doclus=0;
  int readclusters=0;
  int sigdig=5;
  int noempty=0;
  int nopoint=0;
  int rndcorr=0;
  char *ecifile="";
  int multincl=0;
  char *corrfunc_label="trigo";
  char *gce_label="tensor";
  int printnum=0;
  AskStruct options[]={
    {"","Generalized Cluster Expansion " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-h","Display more help",BOOLVAL,&dohelp},
    {"-2","Maximum distance between two points within a pair",REALVAL,&maxd(2)},
    {"-3","Maximum distance between two points within a triplet",REALVAL,&maxd(3)},
    {"-4","Maximum distance between two points within a quadruplet",REALVAL,&maxd(4)},
    {"-5","Maximum distance between two points within a quintuplet",REALVAL,&maxd(5)},
    {"-6","Maximum distance between two points within a sextuplet",REALVAL,&maxd(6)},
    {"-l","Input file defining the lattice   (default: lat.in)",STRINGVAL,&latticefilename},
    {"-s","Input file defining the structure (default: str.out)",STRINGVAL,&strfilename},
    {"-pc","Print composition only",BOOLVAL,&doconc},
    {"-sym","Just find space group",BOOLVAL,&dosym},
    {"-clus","Just find clusters",BOOLVAL,&doclus},
    {"-c","Read clusters.out file instead of writing it",BOOLVAL,&readclusters},
    {"-z","Tolerance for finding symmetry operations (default: 1e-3)",REALVAL,&zero_tolerance},
    {"-sig","Number of significant digits printed (default: 5)",INTVAL,&sigdig},
    {"-noe","Do not include empty cluster",BOOLVAL,&noempty},
    {"-nop","Do not include point cluster(s)",BOOLVAL,&nopoint},
    {"-rnd","Print correlation of the random state of the same composition as the input structure",BOOLVAL,&rndcorr},
    {"-eci","Predict quantity using ECI in specified file",STRINGVAL,&ecifile},
    {"-mi","Multiplicities are already included in ECI file",BOOLVAL,&multincl},
    {"-crf","Select correlation functions (default: trigo)",STRINGVAL,&corrfunc_label},
    {"-gce","Select type of Generalized Cluster Expansion (default: tensor)",STRINGVAL,&gce_label},
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
  if (nbsite==0) ERRORQUIT("Need at least one site with multiple species in lattice file.");
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

  LinkedList<GeneralizedCluster> clusterlist;
  if (doconc) {
    ofstream file("atoms.out");
    for (int i=0; i<label.get_size(); i++) {
      file << label(i) << endl;
    }
  }
  else {
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
      if (!check_plug_in(SymmetryObeyingObject(),gce_label)) {
	ERRORQUIT("Wrong GCE label. Aborting");
      }
      SymmetryObeyingObject *pdummysymobj=GenericPlugIn<SymmetryObeyingObject>::create(gce_label);

      if (! noempty) {
	LinkedList<SymmetryObeyingObject> listobj;
	pdummysymobj->find_basis(&listobj, spacegroup);
	LinkedListIterator<SymmetryObeyingObject> i(listobj);
	for( ; i; i++) {
	  clusterlist << new GeneralizedCluster(MultiCluster(0),i);
	}
	listobj.detach_all();
      }
      for (int m=(nopoint ? 2 : 1); m<=maxp; m++) {
	pclusterbanks(m)->reset();
	for (int c=0; c<nbcluster(m); c++) {
	  SpaceGroup subgroup;
	  find_clusters_symmetry(&subgroup, *pclusterbanks(m),spacegroup);
	  LinkedList<SymmetryObeyingObject> listobj;
	  pdummysymobj->find_basis(&listobj, subgroup);
	  LinkedListIterator<SymmetryObeyingObject> i(listobj);
	  for( ; i; i++) {
	    clusterlist << new GeneralizedCluster(*pclusterbanks(m),i);
	  }
	  listobj.detach_all();
	  (*pclusterbanks(m))++;
	}
	delete pclusterbanks(m);
      }
      delete pdummysymobj;
      {
	ofstream clusterfile("clusters.out");
	clusterfile.setf(ios::fixed);
	clusterfile.precision(sigdig);
	LinkedListIterator<GeneralizedCluster> icluster(clusterlist);
	for ( ; icluster; icluster++) {
	  int mult=calc_multiplicity(icluster->clus, spacegroup.cell, spacegroup.point_op, spacegroup.trans);
	  clusterfile << mult << endl;
	  clusterfile << get_length_quick(icluster->clus.clus) << endl;
	  clusterfile << icluster->clus.clus.get_size() << endl;
	  for (int i=0; i<icluster->clus.clus.get_size(); i++) {
	    clusterfile << ((!axes)*(icluster->clus.clus(i))) << " " << icluster->clus.site_type(i) << " " << icluster->clus.func(i) <<  endl;
	  }
	  clusterfile << gce_label << endl;
	  icluster->psymobj->write(clusterfile);
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

  // calccorr;

  LinkedList<Array<GeneralizedCluster> > eq_clusterlist;
  LinkedListIterator<GeneralizedCluster> icluster(clusterlist);
  for ( ; icluster; icluster++) {
    Array<GeneralizedCluster> *pmulticlus=new Array<GeneralizedCluster>;
    find_equivalent_clusters(pmulticlus, *icluster, spacegroup.cell, spacegroup.point_op, spacegroup.trans);
    eq_clusterlist << pmulticlus;
  }
  
  // initialize a table of correlation functions;
  if (!check_plug_in(CorrFuncTable(),corrfunc_label)) {
    ERRORQUIT("Aborting: invalid correlation function label");
  }
  CorrFuncTable *pcorrfunc=GenericPlugIn<CorrFuncTable>::create(corrfunc_label);
  pcorrfunc->init_from_site_type_list(labellookup);

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
    
    cout.setf(ios::fixed);
    //    cout.setf(ios::showpos);
    cout.precision(sigdig);
    if (printnum) cout << strnum << delim;
    strnum++;
    Real pred=0.;
    if (doconc) {
      Array<Real> conc;
      Structure cstr(ideal_str);
      fix_atom_type(&cstr,lattice,labellookup,0);
      calc_concentration(&conc,lattice,labellookup,cstr);
      for (int i=0; i<conc.get_size(); i++) {
	cout << conc(i) << " ";
      }
    }
    else if (rndcorr) {
      ERRORQUIT("Not implemented!");
      /*
      LinkedList<Real> pointcorr;
      LinkedListIterator<Array<GeneralizedCluster> > ipcluster(eq_clusterlist);
      for ( ; ipcluster; ipcluster++) {
	if ((*ipcluster)(0).clus.get_size() > 1) break;
	if ((*ipcluster)(0).clus.get_size()==1) {
	  pointcorr << new Real(calc_correlation(ideal_str, *ipcluster, spacegroup.cell, *pcorrfunc));
	}
      }
      int ieci=0;
      LinkedListIterator<GeneralizedCluster> icluster(clusterlist);
      LinkedListIterator<Array<GeneralizedCluster> > ieqcluster(eq_clusterlist);
      for ( ; icluster; icluster++, ieqcluster++, ieci++) {
	Real rho=1.;
	for (int s=0; s<icluster->clus.get_size(); s++) {
	  LinkedListIterator<GeneralizedCluster> ipcluster(clusterlist);
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
	  cout << rho << delim;
	}
      }
      */
    }
    else {
      SymmetryObeyingObject *pdummysymobj=GenericPlugIn<SymmetryObeyingObject>::create(gce_label);
      LinkedList<SymmetryObeyingObject> rholist;
      SymmetryObeyingObject *ppred=pdummysymobj->clone();
      ppred->zero();
      int ieci=0;
      LinkedListIterator<Array<GeneralizedCluster> > icluster(eq_clusterlist);
      for ( ; icluster; icluster++, ieci++) {
	SymmetryObeyingObject *prho=pdummysymobj->clone();
	calc_correlation(prho, ideal_str, *icluster, spacegroup.cell, *pcorrfunc);
	rholist << prho;
        if (strlen(ecifile)>0) {
	  ppred->add(prho,eci(ieci)*(multincl ? 1 : icluster->get_size()));
	}
      }
      Array<Real> dummydata;
      ppred->get_data(&dummydata);
      int dimsymobj=dummydata.get_size();

      if (strlen(ecifile)>0) {
	ppred->write(cout);
      }
      else {
	for (int i=0; i<dimsymobj; i++) {
	  LinkedListIterator<Array<GeneralizedCluster> > icluster(eq_clusterlist);
	  LinkedListIterator<SymmetryObeyingObject> irholist(rholist);
	  for (; icluster; icluster++, irholist++) {
	    Array<Real> data;
	    irholist->get_data(&data);
	    cout << data(i) << delim;
	  }
	  if (i<dimsymobj-1) {cout << endl;}
	}
      }
      delete ppred;
      delete pdummysymobj;
    }
    cout << endl;
  }

}

