#include "parse.h"
#include "getvalue.h"
#include "version.h"
#include "linalg.h"
#include "linsolve.h"
#include <fstream.h>

void mat_sqrt(rMatrix3d *pS, const rMatrix3d &S2) {
  Array2d<Real> a2;
  convert_matrix(&a2,S2);
  Array<Real> lambda;
  Array2d<Real> vect;
  diagonalize_symmetric_matrix(&lambda,&vect, a2);
  for (int i=0; i<3; i++) {
    lambda(i)=sqrt(lambda(i));
  }
  Array2d<Real> tmp,ivect,a;
  invert_matrix(&ivect,vect);
  product_diag(&tmp,lambda,ivect);
  product(&a,vect,tmp);
  convert_matrix(pS,a);
}

int main(int argc, char *argv[]) {
  char *strfilename="str.out";
  char *strfilename_rel="str_relax.out";
  int dummy=0;
  AskStruct options[]={
    {"","check cell distortion " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-d","default",BOOLVAL,&dummy},
    {"-su","unrelaxed structure file name",STRINGVAL,&strfilename},
    {"-sr","relaxed structure file name",STRINGVAL,&strfilename_rel}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    return 1;
  }
  Array<Arrayint> site_type_list;
  Array<AutoString> label;
  rMatrix3d axes;

  Structure str;
  {
    ifstream strfile(strfilename);
    parse_lattice_file(&str.cell, &str.atom_pos, &str.atom_type, &site_type_list, &label, strfile, &axes);
    fix_atom_type(&str,site_type_list);
  }
  Structure str_rel;
  {
    ifstream strfile(strfilename_rel);
    parse_lattice_file(&str_rel.cell, &str_rel.atom_pos, &str_rel.atom_type, &site_type_list, &label, strfile);
    fix_atom_type(&str_rel,site_type_list);
  }

  rMatrix3d T=str_rel.cell*(!str.cell);
  rMatrix3d S2=(~T)*T;
  rMatrix3d S;
  mat_sqrt(&S,S2);
  rMatrix3d Id;
  Id.identity();
  S=S-Id;
  cout << S;
}
