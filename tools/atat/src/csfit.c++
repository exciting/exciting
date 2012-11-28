#include <strstream.h>
#include <sys/stat.h>
#include "arraylist.h"
#include "vectmac.h"
#include "getvalue.h"
#include "lstsqr.h"
#include "version.h"
#include "stringo.h"
#include "parse.h"

#define NO_OBJECTS_IN_KSPACE_CS
#include "kspacecs.cc"

void distord_cell(rMatrix3d *pcell, const rVector3d &dir, Real e_paral, Real e_perp) {
  rMatrix3d strain,basis;
  strain.zero();
  strain(0,0)=e_paral;
  strain(1,1)=e_perp;
  strain(2,2)=e_perp;
  basis.set_column(0,dir/norm(dir));
  basis.set_column(1,find_perpendicular(dir));
  basis.set_column(2,basis.get_column(0)^basis.get_column(1));
  (*pcell)=basis*strain*(!basis)*(*pcell);
}

extern char *helpstring;

int main(int argc, char *argv[]) {
  // parse command line;
  int dohelp=0;
  int conc_mesh=50;
  int search_mesh=100;
  int perp_stretch_mesh=5;
  int paral_stretch_mesh=3;
  Real max_paral_stretch=0.05;
  char *puredir[2]={"0","1"};
  char *dirfilename="dir.in";
  char *latticefilename="lat.in";
  int dummy=0;
  int sigdig=6;
  int polltime=10;
  AskStruct options[]={
    {"","Constituent Strain FITter " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-h","Display more help",BOOLVAL,&dohelp},
    {"-nc","Number of points in concentration mesh (default 50)",INTVAL,&conc_mesh},
    {"-ns","Number of points in mesh used to look for energy minimum (default 100)",INTVAL,&search_mesh},
    {"-np","Number of points in stretching mesh in the direction perpendicular to the k-vector(default 5)",INTVAL,&perp_stretch_mesh},
    {"-nl","Number of points in stretching mesh in the direction parallel to the k-vector (default 3)",INTVAL,&paral_stretch_mesh},
    {"-ml","Maximum parallel stretching (default 0.05)",REALVAL,&max_paral_stretch},
    {"-l","Input file defining the lattice (default: lat.in)",STRINGVAL,&latticefilename},
    {"-pa","Directory containing the pure A calculations (default 0/)",STRINGVAL,&puredir[0]},
    {"-pb","Directory containing the pure B calculations (default 1/)",STRINGVAL,&puredir[1]},
    {"-ds","File containing a list of stretching directions (default: dir.in)",STRINGVAL,&dirfilename},
    {"-t","Time between disk reads in sec (default: 10 sec)",INTVAL,&polltime},
    {"-sig","Number of significant digits to print in output files",INTVAL,&sigdig},
    {"-d","Use all default values",BOOLVAL,&dummy}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    return 1;
  }
  if (dohelp) {
    cout << helpstring;
    return 1;
  }
  if ( perp_stretch_mesh<3) ERRORQUIT("-np must be greater than 2");
  if (paral_stretch_mesh<3) ERRORQUIT("-nl must be greater than 2");

  // read in the lattice file (so that we know the conventional cell and associated reciprocal cell);
  Structure lat;
  Array<Arrayint> site_type_list;
  Array<AutoString> atom_label;
  rMatrix3d axes;
  ifstream file(latticefilename);
  if (!file) ERRORQUIT("Unable to open lattice file.");
  parse_lattice_file(&lat.cell, &lat.atom_pos, &lat.atom_type, &site_type_list, &atom_label, file, &axes);
  rMatrix3d rec_axes=~(!axes);

  // read in specified direction in k space;
  LinkedList<rVector3d> dirlist;
  {
    ifstream dirfile(dirfilename);
    if (!dirfile) ERRORQUIT("Unable to open direction file");
    while (1) {
      rVector3d dir(MAXFLOAT,MAXFLOAT,MAXFLOAT);
      dirfile >> dir;
      if (dir(0)==MAXFLOAT || !dirfile || dirfile.eof()) break;
      dirlist << new rVector3d(dir);
    }
  }

  // read in the relaxed geometry of the pure end members;
  Real scale[2];
  Structure str[2];
  for (int elem=0; elem<2; elem++) {
    if (chdir(puredir[elem])!=0) {cerr << "Directory " << puredir[elem]; ERRORQUIT(" does not exist.");}
    {
      ifstream file("str_relax.out");
      if (!file)  {cerr << "Unable to open file" << puredir[elem] << "/str-relax.out" ; ERRORQUIT("");}
      parse_structure_file(&str[elem].cell,&str[elem].atom_pos,&str[elem].atom_type, atom_label,file);
    }
    scale[elem]=pow(det(str[elem].cell)/det(lat.cell),1./3.); // for cubic only!;
    chdir("..");
  }

  // generate distorted structures;
  {
    LinkedListIterator<rVector3d> i_dir(dirlist);
    for (; i_dir; i_dir++) {
      for (int elem=0; elem<2; elem++) {
	chdir(puredir[elem]);
	for (int perp=0; perp<perp_stretch_mesh; perp++) {
	  Real t=(Real)perp/(Real)(perp_stretch_mesh-1);
	  Real stretch_perp=(1.-t)*(scale[0]/scale[elem])+t*(scale[1]/scale[elem]); // interpolate between lattice parameter of each end member;
	  for (int paral=0; paral<paral_stretch_mesh; paral++) {
	    // calc the stretching along k vector: take volume preserving value and add or substract a fraction of maximum parallel strecthing;
	    Real stretch_paral=pow(stretch_perp,-2.) + ( 2.*(Real)paral/(Real)(paral_stretch_mesh-1) - 1. ) * max_paral_stretch;
	    Structure stretched_str(str[elem]);
	    rVector3d dir=rec_axes*(*i_dir);
	    distord_cell(&(stretched_str.cell),dir,stretch_paral,stretch_perp);
	    rMatrix3d transfo=stretched_str.cell*(!str[elem].cell);
	    for (int i=0; i<str[elem].atom_pos.get_size(); i++) {stretched_str.atom_pos(i)=transfo*str[elem].atom_pos(i);}
	    // create directory and write the structure geometry file;
	    ostrstream cur_dir;
	    cur_dir << "stretch_" << (*i_dir)(0) << "_" << (*i_dir)(1) << "_" << (*i_dir)(2) << "_" << perp << "_" << paral << '\0';
	    mkdir(cur_dir.str(),S_IRWXU | S_IRWXG | S_IRWXO);
	    if (chdir(cur_dir.str())!=0) {cerr << "Unable to create or cd into " << puredir[elem] << "/" << cur_dir.str(); ERRORQUIT("");}
	    {
	      ofstream file("str.out");
	      file.setf(ios::fixed);
	      file.precision(sigdig);
	      write_structure(stretched_str,atom_label,transfo*scale[elem]*axes,file);
	    }
	    if ( ! file_exists("energy")) {
	      ofstream file("wait"); ///flag for pollmach;
	    }
	    chdir("..");
	  }
	}
	chdir("..");
      }
    }
  }
  cerr << "Structure generation completed." << endl;

  while ( file_exists("ready") ) {unlink("ready");} // flag for pollmach;

  // read in energies, retrying until all of them are available;
  while (1) {
    int nomissing=1;
    LinkedListIterator<rVector3d> i_dir(dirlist);
    for (; nomissing && i_dir; i_dir++) {
      for (int elem=0; nomissing && elem<2; elem++) {
	chdir(puredir[elem]);
	for (int perp=0; nomissing && perp<perp_stretch_mesh; perp++) {
	  for (int paral=0; nomissing && paral<paral_stretch_mesh; paral++) {
	    ostrstream cur_dir;
	    cur_dir << "stretch_" << (*i_dir)(0) << "_" << (*i_dir)(1) << "_" << (*i_dir)(2) << "_" << perp << "_" << paral << '\0';
	    if (chdir(cur_dir.str())!=0) {cerr << "Unable to create or cd into " << puredir[elem] << "/" << cur_dir.str(); ERRORQUIT("");}
	    if (!file_exists("energy")) {nomissing=0;}
	    else {
	      ifstream file("energy");
	      Real e=MAXFLOAT;
	      file >> e;
	      if (e==MAXFLOAT || file.fail()) {
		cerr << "Error reading energy file in " << cur_dir.str() << endl;
		nomissing=0;
	      }
	    }
	    chdir("..");
	  }
	}
	chdir("..");
      }
    }
    if (nomissing) break;
    sleep(polltime);
  }

  cerr << "Calculations done, fitting the data" << endl;
  LinkedList<Array<Array<Real> > > perp_e_list;
  {
    ofstream debugfile("csdebug.out");
    LinkedListIterator<rVector3d> i_dir(dirlist);
    for (; i_dir; i_dir++) {
      Array<Array<Real> > perp_e(2);
      for (int elem=0; elem<2; elem++) {
	perp_e(elem).resize(perp_stretch_mesh);
	chdir(puredir[elem]);
	for (int perp=0; perp<perp_stretch_mesh; perp++) {
	  Array<Real> paral_e(paral_stretch_mesh);
	  for (int paral=0; paral<paral_stretch_mesh; paral++) {
	    ostrstream cur_dir;
	    cur_dir << "stretch_" << (*i_dir)(0) << "_" << (*i_dir)(1) << "_" << (*i_dir)(2) << "_" << perp << "_" << paral << '\0';
	    chdir(cur_dir.str());
	    {
	      ifstream file("energy");
	      file >> paral_e(paral);
	      debugfile << paral_e(paral) << " ";
	    }
	    chdir("..");
	  }
	  perp_e(elem)(perp)=find_minimum(paral_e,search_mesh);
	  debugfile << perp_e(elem)(perp) << endl;
	}
	chdir("..");
      }
      perp_e_list << new Array<Array<Real> >(perp_e);
    }
  }

  Array<Real> pure_e[2];
  for (int elem=0; elem<2; elem++) {
    pure_e[elem].resize(perp_e_list.get_size());
    LinkedListIterator<Array<Array<Real> > > i_perp_e(perp_e_list);
    for (int index_dir=0; i_perp_e; i_perp_e++,index_dir++) {
      pure_e[elem](index_dir)=find_minimum((*i_perp_e)(elem),search_mesh);
    }
  }

  Array<Array<Real> > best_kubic_coef;
  Real best_cv=MAXFLOAT;
  for (int n_kubic=2; n_kubic<=MIN(4,dirlist.get_size()-1); n_kubic++) {
    ofstream file("cs.log");
    Array<Array<Real> > kubic_coef(conc_mesh);
    Real cv=0.;
    for (int ic=0; ic<conc_mesh; ic++) {
      Real c=(Real)ic/(Real)(conc_mesh-1);
      Array<Real> cs_e(dirlist.get_size());
      Array2d<Real> kubic_val(dirlist.get_size(),n_kubic);
      LinkedListIterator<rVector3d> i_dir(dirlist);
      LinkedListIterator<Array<Array<Real> > > i_perp_e(perp_e_list);
      file << c << " ";
      for (int index_dir=0; i_dir; i_dir++,i_perp_e++,index_dir++) {
	Array<Real> perp_e((*i_perp_e)(0).get_size());
	for (int i=0; i<perp_e.get_size(); i++) {
	  perp_e(i)=(*i_perp_e)(0)(i)*(1.-c)+(*i_perp_e)(1)(i)*c;
	}
	cs_e(index_dir)=find_minimum(perp_e,search_mesh)-(pure_e[0](index_dir)*(1.-c)+pure_e[1](index_dir)*c);
	file << cs_e(index_dir) << " ";
	if (near_zero(c) || near_zero(1.-c)) {
	  cs_e(index_dir)=0.;
	}
	else {
	  cs_e(index_dir)/=4.*c*(1.-c);
	}
	rVector3d dir=rec_axes*(*i_dir);
	for (int k=0; k<n_kubic; k++) {
	  kubic_val(index_dir,k)=kubic_harm(dir,k);
	}
      }
      cv+=calc_cv(kubic_val,cs_e);
      calc_ols(&kubic_coef(ic),kubic_val,cs_e);
      //      for (int k=0;  k<n_kubic; k++) {
      //	file << kubic_coef(ic)(k) << " ";
      //      }
      file << endl;
    }
    if (cv<best_cv) {
      best_cv=cv;
      best_kubic_coef=kubic_coef;
    }
  }

  {
    ofstream file("cs.in");
    file.setf(ios::fixed);
    file.precision(sigdig);
    file << 0. << endl;
    file << best_kubic_coef(0).get_size() << endl;
    file << best_kubic_coef.get_size() << endl;
    for (int i=0; i<best_kubic_coef(0).get_size(); i++) {
      for (int j=0; j<best_kubic_coef.get_size(); j++) {
	file << best_kubic_coef(j)(i) << endl;
      }
      file << endl;
    }
  }
}
