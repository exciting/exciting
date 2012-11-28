#include <iostream.h>
#include "vectmac.h"
#include "array.h"
#include "findsym.h"
#include "getvalue.h"

void read_matrix(rMatrix3d *m, istream &infile) {
  for (int i=0; i<3; i++) {
    rVector3d tmp;
    infile >> tmp;
    m->set_column(i,tmp);
  }
}

void write_matrix(const rMatrix3d &m, ostream &outfile) {
  for (int i=0; i<3; i++) {
    outfile << m.get_column(i) << endl;
  }
}

int main(int argc, char *argv[]) {
  int alldefault=0,quiet=0,round=0, even=0;
  AskStruct options[]={
    {"-d","Use all default options",BOOLVAL,&alldefault},
    {"-q","Quiet!",BOOLVAL,&quiet},
    {"-r","Round output to next largest integer",BOOLVAL,&round},
    {"-e","Force nb of kpoints to be an even number",BOOLVAL,&even}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    return 1;
  }
  if (even) {round=1;}

  if (!quiet) {cout << "Number of k-points desired:" << endl;}
  Real nb_kpts;
  cin >> nb_kpts;
  rMatrix3d cell;
  if (!quiet) {cout << "Cell: (in cartesian, one lattice vector per line with coordinates separated by spaces)" << endl;}
  read_matrix(&cell,cin);
  rMatrix3d reciprocal=~(!cell);
  rVector3d proj;
  for (int i=0; i<3; i++) {
    rVector3d a=reciprocal.get_column( i     );
    rVector3d b=reciprocal.get_column((i+1)%3);
    rVector3d c=reciprocal.get_column((i+2)%3);
    proj(i)=fabs(a*((b^c)/norm(b^c)));
  }
  {
    iVector3d equals(0,0,0);
    rMatrix3d ireciprocal=!reciprocal;
    Array<rMatrix3d> point_op;
    find_pointgroup(&point_op,reciprocal);
    for (int op=0; op<point_op.get_size(); op++) {
      rMatrix3d m=ireciprocal*point_op(op)*reciprocal;
      for (int i=0; i<3; i++) {
	if (!near_zero(m(i,(i+1)%3))) {
	  equals(i)=1;
	}
      }
    }
    if (l1_norm(equals)>1) {
      for (int i=0; i<3; i++) {
	proj(i)=1.;
      }
    } else {
      for (int i=0; i<3; i++) {
	if (equals(i)==1) {
	  Real avg=sqrt(proj(i)*proj((i+1)%3));
	  proj(i)=avg;
	  proj((i+1)%3)=avg;
	}
      }
    }

  }
  Real normalizer=pow(nb_kpts/(proj(0)*proj(1)*proj(2)),1./3.);
  if (!quiet) {cout << "Suggested division in k point mesh for each axis:" << endl;}
  rVector3d mesh=proj*normalizer;
  if (round) {
    iVector3d imesh=to_int(mesh);
    rVector3d fmesh=mesh-to_real(imesh);
    if (even) {
      for (int j=0; j<3; j++) {
	if (imesh(j)%2) {
	  imesh(j)++;
	  fmesh(j)=0.;
	}
      }
    }
    while (imesh(0)*imesh(1)*imesh(2)<nb_kpts) {
      Real bestf=0.;
      for (int i=0; i<3; i++) {
	if (fmesh(i)>bestf) {
	  bestf=fmesh(i);
	}
      }
      for (int i=0; i<3; i++) {
	if (fabs(fmesh(i)-bestf)<zero_tolerance) {
	  imesh(i)++;
	  if (even) imesh(i)++;
	  fmesh(i)=0.;
	}
      }
    }
    cout << imesh << endl;
  }
  else {
    cout << mesh << endl;
  }
}
