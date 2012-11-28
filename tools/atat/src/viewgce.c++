#include <fstream.h>
#include "tensor.h"
#include "getvalue.h"
#include "version.h"
#include "parse.h"

int equal_multiclus(const MultiCluster &a, const MultiCluster &b) {
  if (a.clus.get_size()!=b.clus.get_size()) return 0;
  for (int i=0; i<a.clus.get_size(); i++) {
    if ( (a.clus(i)!=b.clus(i)) || (a.site_type(i)!=b.site_type(i)) || (a.func(i)!=b.func(i)) ) return 0;
  }
  return 1;
}

void view_tensor_cluster(const rTensor &sumtensor, const MultiCluster &oldclus, int c, char *plotcommand, int normalize=0) {
  {
    ofstream pfile("clusp.tmp");
    ofstream lfile("clusl.tmp");
    lfile << 0 << " " << 0 << " " << 0 << endl << endl;
    for (int i=0; i<oldclus.clus.get_size(); i++) {
      pfile << oldclus.clus(i) << endl;
      lfile << 0 << " " << 0 << " " << 0 << endl;
      lfile << oldclus.clus(i)(0) << " " << 0 << " " << 0 << endl;
      lfile << oldclus.clus(i)(0) << " " << oldclus.clus(i)(1) << " " << 0 << endl;
      lfile << oldclus.clus(i)(0) << " " << oldclus.clus(i)(1) << " " << oldclus.clus(i)(2) << endl;
      lfile << endl << endl;
    }
  }
  LinkedList<rVector3d> fixx_list;
  if (file_exists("fixdir.in")) {
    ifstream file("fixdir.in");
    while (!file.eof()) {
      rVector3d tmp(MAXFLOAT,MAXFLOAT,MAXFLOAT);
      file >> tmp;
      if (tmp(0)==MAXFLOAT) break;
      fixx_list << new rVector3d(tmp);
    }
  }
  {
    ofstream shfile("sh.tmp");
    Real maxr=0.;
    for (int pass=0; pass<2; pass++) {
      Real dt=M_PI/48;
      Real dp=M_PI/48;
      for (Real t=0; t<M_PI+dt/2; t+=dt) {
	for (Real p=0; p<2*M_PI+dp/2; p+=dp) {
	  Array<Real> x(3);
	  x(0)=sin(t)*cos(p);
	  x(1)=sin(t)*sin(p);
	  x(2)=cos(t);
	  Real r=0.;
	  const Array<int> &size=sumtensor.get_size();
	  MultiDimIterator<Array<int> > i(size);
	  for ( ; i; i++) {
	    Real a=sumtensor(i);
	    LinkedListIterator<rVector3d> ifixx(fixx_list);
	    for (int j=0; j<size.get_size(); j++) {
	      if (ifixx) {
		a*=(*ifixx)(((Array<int> &)i)(j));
		ifixx++;
	      }
	      else {
		a*=x(((Array<int> &)i)(j));
	      }
	    }
	    r+=a;
	  }
	  maxr=MAX(maxr,fabs(r));
	  if (pass==1) {
	    shfile << x(0) << " " << x(1) << " " << x(2) << " " << (normalize ? r/maxr : r) << endl;
	  }
	}
	if (pass==1) {
	  shfile << endl;
	}
      }
    }
  }
  ostrstream line;
  line << plotcommand << " ";
  line.width(3);
  line.fill('0');
  line << c << '\0';
  system(line.str());
}

int main(int argc, char *argv[]) {
  char *latticefilename="lat.in";
  char *ecifilename="";
  char *rmcommand="rm -f shclus???.pnm";
  char *plotcommand="viewtensorclus";
  char *catcommand="pnmcat -topbottom shclus???.pnm | pnmtopng > all.png";
  int normalize=0;
  int dummy=0;
  AskStruct options[]={
    {"","Viewer for a Generalized Cluster Expansion " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-l","Input file defining the lattice   (default: lat.in)",STRINGVAL,&latticefilename},
    {"-eci","Combine terms using ECI in specified file",STRINGVAL,&ecifilename},
    {"-rc","Remove command to run beforehand",STRINGVAL,&rmcommand},
    {"-cc","Command to catenate all output",STRINGVAL,&catcommand},
    {"-pc","Command to plot tensor and cluster",STRINGVAL,&plotcommand},
    {"-n","Normalize plots by max",BOOLVAL,&normalize},
    {"-d","All default values",BOOLVAL,&dummy}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
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

  system(rmcommand);

  LinkedList<GeneralizedCluster> clusterlist;
  LinkedList<Real> ecilist;

  {
    ifstream clusterfile("clusters.out");
    if (!clusterfile) ERRORQUIT("Unable to open clusters.out file");

    if (strlen(ecifilename)>0) {
      ifstream ecifile(ecifilename);
      if (!ecifile) ERRORQUIT("Unable to open eci file");
      read_clusters_and_eci(&clusterlist,&ecilist,clusterfile,ecifile,axes);
    }
    else {
      read_clusters_and_eci(&clusterlist,NULL,clusterfile,clusterfile,axes);
    }
  }

  LinkedListIterator<GeneralizedCluster> ic(clusterlist);
  LinkedListIterator<Real> ie(ecilist);
  if (ie) {
    for (int c=0; ic; c++) {
      rTensor sumtensor(((TensorSymmetryObeyingObject *)(ic->psymobj))->tensor.get_size());
      sumtensor.zero();
      LinkedListIterator<GeneralizedCluster> jc=ic;
      while (jc) {
	if (!equal_multiclus(ic->clus,jc->clus)) break;
	rTensor tensor=((TensorSymmetryObeyingObject *)(jc->psymobj))->tensor;
	tensor*=(*ie);
	sumtensor+=tensor;
	jc++;
	ie++;
      }
      view_tensor_cluster(sumtensor,ic->clus,c,plotcommand,normalize);
      ic=jc;
    }
  }
  else {
    for (int c=0; ic; ic++,c++) {
      view_tensor_cluster(((TensorSymmetryObeyingObject *)(ic->psymobj))->tensor,ic->clus,c,plotcommand,normalize);
    }
  }
  system(catcommand);
}
