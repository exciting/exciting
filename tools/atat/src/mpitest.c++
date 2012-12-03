#include "mpiinterf.h"
//#include "xtalutil.h"
//#include "parse.h"

void mysum(int *a, const int &b, const int &c) {
  *a=b+c;
}

int main(int argc, char **argv) {
  MyMPIobj.init(argc,argv);

  int v=MyMPIobj.id+3;
  MyMPI_BcastStream(&v,0);

  //  MyMPI_Reduce(&v,mysum);
  cout << MyMPIobj.id << " " << v << endl;
}

/*
int main(int argc, char **argv) {
  MyMPIobj.init(argc,argv);
  int n=15;
  Array<Structure> a(n);
  {
      MPISynchronizer<Structure> sync;
      for (int i=0; i<n; i++) {
	  if (sync.is_my_job()) {
	      a(i).cell.identity();
	      a(i).cell=a(i).cell*(Real)i;
	      a(i).atom_pos.resize(i);
	      a(i).atom_type.resize(i);
	      for (int j=0; j<i; j++) {
		  a(i).atom_pos(j)=rVector3d(1.,2.,3.);
		  a(i).atom_type(j)=j;
	      }
	  }
	  sync.sync(&(a(i)));
      }
  }
  rMatrix3d axes;
  axes.identity();
  Array<AutoString> labels(n);
  for (int i=0; i<n; i++) {
      labels(i).set(1);
      labels(i)[0]='a'+i;
  }  
  if (MyMPIobj.id==0) {
      for (int i=0; i<n; i++) {
	  write_structure(a(i), labels , axes, cout);
      }
  }
}
*/

/*
  Array<Array<Real> > x(n);
  for (int i=0; i<n; i++) {
      x(i).resize(i);
      zero_array(&x(i));
  }
  {
    MPISynchronizer<Array<Real> > sync;
    for (int i=0; i<n; i++) {
      if (sync.is_my_job(i)) {
	for (int j=0; j<i; j++) {
	  x(i)(j)=(Real)i+2.5*j+3.0;
	}
      }
      sync.sync(&(x(i)),i);
      if (MyMPIobj.is_root()) {
	cout << "===(" << endl << x << ")===" << endl;
      }
    }
  }
  if (MyMPIobj.is_root()) {
    for (int i=0; i<n; i++) {
      cout << "id=" << MyMPIobj.id << " " << i << " " << x(i) << endl;
    }
  }
*/

