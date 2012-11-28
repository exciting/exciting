#include <fstream.h>
#include "array.h"
#include "getvalue.h"
#include "xtalutil.h"
#include "parse.h"
#include "version.h"

Real angle(const rVector3d &a, const rVector3d &b) {
  return (acos(a*b/(norm(a)*norm(b))));
}

int add_if_nonequiv(LinkedList<Array<rVector3d> > *pclist, const Array<rVector3d> &cells, const Array<rMatrix3d> &point_op1, const Array<rMatrix3d> &point_op2) {
  LinkedListIterator<Array<rVector3d> > i(*pclist);
  for ( ; i; i++) {
    int same1=0;
    for (int op1=0; op1<point_op1.get_size(); op1++) {
      if (near_zero(norm((*i)(0)-point_op1(op1)*cells(0))) && near_zero(norm((*i)(1)-point_op1(op1)*cells(1)))) {
	same1=1;
	break;
      }
    }
    int same2=0;
    for (int op2=0; op2<point_op2.get_size(); op2++) {
      if (near_zero(norm((*i)(2)-point_op2(op2)*cells(2))) && near_zero(norm((*i)(3)-point_op2(op2)*cells(3)))) {
	same2=1;
	break;
      }
    }
    if (same1 && same2) break;
  }
  if (!i) {
    (*pclist) << new Array<rVector3d>(cells);
    return 1;
  }
  else {
    return 0;
  }
}

Real get_2d_strain(const rMatrix2d s) {
  Real d=det(s);
  Real t=trace(s);
  return (fabs(t)+sqrt(t*t-4*d))/2.;
}

int main(int argc, char *argv[]) {
  char *cellfilename1="";
  char *cellfilename2="";
  char *planefilename="";
  Real rmax=0;
  Real eps=zero_tolerance;
  AskStruct options[]={
    {"","Coincidence Site Lattice " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-l1","Input file defining the first lattice",STRINGVAL,&cellfilename1},
    {"-l2","Input file defining the second lattice",STRINGVAL,&cellfilename2},
    {"-p","Input file defining the planes",STRINGVAL,&planefilename},
    {"-ma","maximum lattice parameter",REALVAL,&rmax},
    {"-e","maximum strain",REALVAL,&eps}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    return 1;
  }
  rMatrix3d cell1,cell2;
  rMatrix3d axes1,axes2;
  {
    ifstream cellfile(cellfilename1);
    if (!cellfile) ERRORQUIT("Unable to open cell file 1");
    read_cell(&axes1,cellfile);
    read_cell(&cell1,cellfile);
    cell1=axes1*cell1;
  }
  {
    ifstream cellfile(cellfilename2);
    if (!cellfile) ERRORQUIT("Unable to open cell file 2");
    read_cell(&axes2,cellfile);
    read_cell(&cell2,cellfile);
    cell2=axes2*cell2;
  }
  rVector3d plane1,plane2;
  {
    ifstream planefile(planefilename);
    if (!planefile) ERRORQUIT("Unable to open plane file");
    planefile >> plane1;
    planefile >> plane2;
    plane1=(!~axes1)*plane1;
    plane2=(!~axes2)*plane2;
  }
  Array<rMatrix3d> point_op1,point_op2;
  find_pointgroup(&point_op1, cell1);
  find_pointgroup(&point_op2, cell2);

  LinkedList<Array<rVector3d> > vlist;
  LatticePointIterator a1(cell1,1);
  for (; norm(a1)<=rmax; a1++) {
    if (near_zero(plane1*a1)) {
      Real r1=norm(a1);
      Real rmax1=r1*(1+eps);
      LatticePointIterator a2(cell2,1);
      for (; norm(a2)<=rmax1; a2++) {
	if (near_zero(plane2*a2)) {
	  if (fabs(norm(a2)-r1)/r1 <= 2.*eps) {
	    Array<rVector3d> *pv=new Array<rVector3d>(2);
	    (*pv)(0)=a1;
	    (*pv)(1)=a2;
	    vlist << pv;
	  }
	}
      }
    }
  }
  cout.setf(ios::fixed);
  cout.precision(6);
  rMatrix2d id;
  id.identity();
  rMatrix3d iaxes1=!axes1;
  rMatrix3d iaxes2=!axes2;
  LinkedList<Array<rVector3d> > clist;
  LinkedListIterator<Array<rVector3d> > i(vlist);
  for ( ; i; i++) {
    Real r1=norm((*i)(0));
    LinkedListIterator<Array<rVector3d> > j(vlist);
    for ( ; j && norm((*j)(0))<=r1+zero_tolerance; j++) {
      Real p=((*j)(0)*(*i)(0))/((*i)(0)*(*i)(0));
      if (fabs(p)<=0.5) {
	Real ang1=angle((*j)(0),(*i)(0));
	if (!near_zero(ang1) && !near_zero(ang1-M_PI)) {
	  rMatrix2d supcell1,supcell2;
	  supcell1.set_column(0,norm((*i)(0))*rVector2d(cos(ang1/2), sin(ang1/2)));
	  supcell1.set_column(1,norm((*j)(0))*rVector2d(cos(ang1/2),-sin(ang1/2)));
	  Real ang2=angle((*j)(1),(*i)(1));
	  supcell2.set_column(0,norm((*i)(1))*rVector2d(cos(ang2/2), sin(ang2/2)));
	  supcell2.set_column(1,norm((*j)(1))*rVector2d(cos(ang2/2),-sin(ang2/2)));
	  Real strain=get_2d_strain((!supcell2)*supcell1-id);
	  if (strain<=eps) {
	    Array<rVector3d> lats(4);
	    lats(0)=(*i)(0);
	    lats(1)=(*j)(0);
	    lats(2)=(*i)(1);
	    lats(3)=(*j)(1);
	    if (add_if_nonequiv(&clist,lats,point_op1,point_op2)) {
	      // cout << supcell1 << endl;
	      // cout << supcell2 << endl;
	      cout << strain << endl;
	      cout << norm((iaxes1*(*i)(0))^(iaxes1*(*j)(0))) << " " << norm((iaxes2*(*i)(1))^(iaxes2*(*j)(1))) << endl;
	      cout << iaxes1*(*i)(0) << " " << iaxes2*(*i)(1) << endl; 
	      cout << iaxes1*(*j)(0) << " " << iaxes2*(*j)(1) << endl; 
	      cout << endl;
	    }
	  }
	}
      }
    }
  }
}
