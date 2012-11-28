#include <fstream.h>
#include "parse.h"
#include "getvalue.h"
#include "version.h"
#include "lstsqr.h"
#include "xtalutil.h"

int unroll_tri(int s,int t) {
    int i=MAX(s,t);
    int j=MIN(s,t);
    return (i*(i+1)/2+j);
}

int main(int argc, char *argv[]) {
  char *delim="\t";
  int dohelp=0;
  char *strfilename="str.out";
  int sigdig=5;
  Real kr=0.5;
  int nosym=0;
  AskStruct options[]={
    {"","Electron and Hole Masses " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-h","Display more help",BOOLVAL,&dohelp},
    {"-s","Input file defining the structure (default: str.out)",STRINGVAL,&strfilename},
    {"-ns","Turn off symmetry",BOOLVAL,&nosym},
    {"-z","Tolerance for finding symmetry operations (default: 1e-3)",REALVAL,&zero_tolerance},
    {"-kr","K-space radius",REALVAL,&kr},
    {"-sig","Number of significant digits printed (default: 5)",INTVAL,&sigdig}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    return 1;
  }
  ofstream logfile("ehmass.log");
  logfile.setf(ios::fixed);
  logfile.precision(sigdig);
  cout.setf(ios::fixed);
  cout.precision(sigdig);

  Structure str;
  Array<Arrayint> labellookup;
  Array<AutoString> label;
  rMatrix3d axes;
  {
    ifstream strfile(strfilename);
    if (!strfile) ERRORQUIT("Unable to open structure file");
    parse_lattice_file(&str.cell, &str.atom_pos, &str.atom_type, &labellookup, &label, strfile, &axes);
  }
  SpaceGroup spacegroup;
  spacegroup.cell=str.cell;
  find_spacegroup(&spacegroup.point_op,&spacegroup.trans,str.cell,str.atom_pos,str.atom_type);
//cerr << spacegroup.point_op.get_size() << endl;
  rMatrix3d rec_lat=!(~str.cell);
  rMatrix3d inv_rec_lat=!rec_lat;
  LinkedList<rVector3d> kpoint_list;
  LinkedList<rVector2d> E_list;
  while (1) {
      rVector3d kp;
      rVector2d E;
      cin >> kp;
//cerr << kp << endl;
      if (cin.eof()) break;
      cin >> E;
      kp=flip_into_brillouin_1(rec_lat*kp,rec_lat);
      if (norm(kp)<kr) {
	  if (nosym) {
	      kpoint_list << new rVector3d(kp);
	      E_list << new rVector2d(E);
	  }
	  else {
//cerr << "in" << endl;
	      for (int op=0; op<spacegroup.point_op.get_size(); op++) {
		  for (Real s=-1.; s<=1.+zero_tolerance; s+=2.) {
		      rVector3d tkp=s*spacegroup.point_op(op)*kp;
		      LinkedListIterator<rVector3d> i(kpoint_list);
		      for (; i; i++) {
			  if (cylinder_norm(inv_rec_lat*(tkp-(*i)))<zero_tolerance) break;
		      }
		      if (!i) {
			  kpoint_list << new rVector3d(tkp);
			  E_list << new rVector2d(E);
//		      cerr << tkp << " " << E(0) << " " << E(1) << endl;
		      }
		  }
	      }
	  }
      }
  }
//cerr << "kp read" << endl;
  Array2d<Real> mat(kpoint_list.get_size(),7);
  logfile << "nb of k points: " << kpoint_list.get_size() << endl << endl;;
  LinkedListIterator<rVector3d> ik(kpoint_list);
  for (int j=0; ik; ik++,j++) {
      int i=0;
      for (int t=0; t<3; t++) {
	  for (int s=0; s<=t; s++) {
	      mat(j,i)=(*ik)(t)*(*ik)(s);
	      i++;
	  }
      }
      mat(j,6)=1.;
  }
  char carrier[2]={'h','e'};
  for (int b=0; b<2; b++) {
      Array<Real> vec(E_list.get_size());
      LinkedListIterator<rVector2d> iE(E_list);
      for (int j=0; iE; iE++,j++) {
	  vec(j)=(*iE)(b);
      }
//      cerr << mat << endl;
//      cerr << vec << endl;
      Array<Real> imass;
      calc_ols(&imass,mat,vec);
      Array<Real> pvec,dvec;
      product(&pvec,mat,imass);
      diff(&dvec,vec,pvec);

      Real n=(Real)(vec.get_size());
      Real e2=inner_product(dvec,dvec)/n;
      logfile << "Carrier: " << carrier[b] << " rms= " << sqrt(e2);
      Real R2=1.-e2/(inner_product(vec,vec)/n-sqr(mean(vec)));
      logfile << " R^2= " << R2 << endl;
      logfile << "k_x k_y k_z predicted_E actual_E diff_e" << endl;
      LinkedListIterator<rVector3d> ik(kpoint_list);
      for (int j=0; ik; ik++,j++) {
	  logfile << (rVector3d)(*ik) << " " << pvec(j) << " " << vec(j) << " " << dvec(j) << endl;
      }
      logfile << endl << "end" << endl;

      rMatrix3d mass;
      for (int t=0; t<3; t++) {
	  for (int s=0; s<3; s++) {
	      mass(s,t)=imass(unroll_tri(s,t));
	  }
      }
      mass=!mass;
      {
	  ostrstream filename;
	  filename << carrier[b] << "mass" << '\0';
	  ofstream file(filename.str());
	  file.setf(ios::fixed);
	  file.precision(sigdig);
	  for (int t=0; t<3; t++) {
	      for (int s=0; s<3; s++) {
		  file << mass(t,s) << delim;
	      }
	      file << endl;
	  }
      }
  }
}
