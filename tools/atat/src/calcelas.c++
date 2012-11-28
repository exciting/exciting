#include <fstream.h>
#include <sys/stat.h>
#include "getvalue.h"
#include "parse.h"
#include "version.h"
#include "findsym.h"
#include "arraylist.h"
#include "tensorsym.h"
#include "lstsqr.h"

//extern char *helpstring;
char *helpstring="";

void convert(rMatrix3d *mat, const rTensor &ten) {
  Array<int> i(2);
  for (i(0)=0; i(0)<3; i(0)++) {
    for (i(1)=0; i(1)<3; i(1)++) {
      (*mat)(i(0),i(1))=ten(i);
    }
  }
}

int main(int argc, char *argv[]) {
  rMatrix3d Id;
  Id.identity();
  char *strfilename="str_relax.out";
  Real max_strain=0.01;
  int nb_strain=1;
  int dohelp=0;
  int sigdig=5;
  int fit=0;
  int dummy=0;
  AskStruct options[]={
    {"","CALCulate ELAStic constants " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-f","Fit elastic constant (instead of the default, generate perturbations)",BOOLVAL,&fit},
    {"-s","Input file defining the relaxed structure (default: str_relax.out)",STRINGVAL,&strfilename},
    {"-ms","Maximum strain magnitude imposed (default: 0.01)",REALVAL,&max_strain},
    {"-ns","Number of strain values (default: 1)",INTVAL,&nb_strain},
    {"-sig","Number of significant digits printed (default: 5)",INTVAL,&sigdig},
    {"-z","Tolerance for finding symmetry operations (default: 1e-3)",REALVAL,&zero_tolerance},
    {"-h","Display more help",BOOLVAL,&dohelp},
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

  Structure str;
  Array<AutoString> atom_label;
  rMatrix3d axes;
  {
    Array<Arrayint> site_type_list;
    ifstream file(strfilename);
    if (!file) ERRORQUIT("Unable to open relaxed structure file.");
    parse_lattice_file(&str.cell, &str.atom_pos, &str.atom_type, &site_type_list, &atom_label, file, &axes);
    if (fabs(det(str.cell))<zero_tolerance) ERRORQUIT("Lattice vectors are coplanar.");
    wrap_inside_cell(&str.atom_pos,str.atom_pos,str.cell);
    fix_atom_type(&str,site_type_list);
  }

  SpaceGroup spacegroup;
  spacegroup.cell=str.cell;
  find_spacegroup(&spacegroup.point_op,&spacegroup.trans,str.cell,str.atom_pos,str.atom_type);
  
  SpaceGroup pointgroup;
  pointgroup_from_spacegroup(&pointgroup.point_op, spacegroup.point_op);

  if (!fit) {
    Array<Array<int> > symflip(1);
    symflip(0).resize(2);  symflip(0)(0)=0;  symflip(0)(1)=1;
    Array<rTensor> pert_strain_t;
    find_symmetry_breaking_basis(&pert_strain_t, 2,pointgroup, symflip);
    
    for (int pert=0; pert<pert_strain_t.get_size(); pert++) {
      for (int m=1; m>=-1; m-=2) {
	rMatrix3d strain;
	convert(&strain,pert_strain_t(pert));
        strain=max_strain*((Real)m)*strain;
	Structure sstr;
	strain_str(&sstr,str,strain);
	// check if negative pert equivalent;
	if (m==-1) {
          rMatrix3d ipluscell=!((Id-strain)*str.cell);
	  int op=0;
	  for ( ; op<pointgroup.point_op.get_size(); op++) {
            if (is_int(ipluscell*pointgroup.point_op(op)*sstr.cell)) break;
	  }
	  if (op<pointgroup.point_op.get_size()) break;
	}
	
	ostrstream pertname;
	pertname << "p" << (m==-1 ? '-' : '+' ) << max_strain << "_" << pert << '\0';
	cerr << " " << pertname.str() << endl;
	mkdir(pertname.str(),S_IRWXU | S_IRWXG | S_IRWXO);
	chdir_robust(pertname.str());
	{
	  ofstream unpertstrfile("str.out");
	  unpertstrfile.setf(ios::fixed);
	  unpertstrfile.precision(sigdig);
	  write_structure(sstr, atom_label, (Id+strain)*axes, unpertstrfile);
	  ofstream waitfile("wait");
	}
	{
	  ofstream strainfile("strain.out");
	  for (int i=0; i<3; i++) {
	    for (int j=0; j<3; j++) {
	      strainfile << strain(i,j) << " ";
	    }
	    strainfile << endl;
	  }
	}
	chdir_robust("..");
      }
    }
  }
  else { // fit instead of generate pert;
    cerr << "Reading..." << endl;
    system("ls p[-+]*/str_relax.out 2> /dev/null | sed 's+/str_relax.out++g' > pertlist.out");
    LinkedList<rMatrix3d> stress_list;
    LinkedList<rMatrix3d> strain_list;
    {
      ifstream volfile("pertlist.out");
      while (skip_delim(volfile)) {
        AutoString volname;
        get_string(&volname,volfile);
	cerr << volname;
        chdir_robust(volname);

        Structure sstr;
        ifstream file("str_relax.out");
        parse_structure_file(&(sstr.cell), &(sstr.atom_pos), &(sstr.atom_type), atom_label, file);
	if (file_exists("stress.out")) {
	  rMatrix3d strain=sstr.cell*(!str.cell)-Id;
	  strain_list << new rMatrix3d(strain);
          ifstream stressfile("stress.out");
	  rMatrix3d stress;
	  for (int i=0; i<3; i++) {
	    for (int j=0; j<3; j++) {
	      stressfile >> stress(i,j);
	    }
	  }
          stress_list << new  rMatrix3d(stress);
	  cerr << endl;
        }
	else {
	  cerr << " cannot read stress.out file" << endl;
	}

	chdir_robust("..");
      }
    }
    unlink("pertlist.out");

    Array<Array<int> > symflip(3);
    symflip(0).resize(2);  symflip(0)(0)=0;  symflip(0)(1)=1;
    symflip(1).resize(2);  symflip(1)(0)=2;  symflip(1)(1)=3;
    symflip(2).resize(4);  symflip(2)(0)=0;  symflip(2)(1)=2;  symflip(2)(2)=1;  symflip(2)(3)=3;
    Array<rTensor> elas_basis;
    find_symmetric_basis(&elas_basis, 4,pointgroup, symflip);
    
    int block=9;
    if (stress_list.get_size()==0) {ERRORQUIT("No stress data found in p[-+]*/stress.out");}
    Array<Real> big_stress(stress_list.get_size()*block);
    Array2d<Real> big_eqnmat(stress_list.get_size()*block,elas_basis.get_size());
    zero_array(&big_stress);
    zero_array(&big_eqnmat);
    
    Array<int> size33(2);
    size33(0)=3; size33(1)=3;
    Array<int> size3333(4);
    size3333(0)=3; size3333(1)=3; size3333(2)=3; size3333(3)=3;
    LinkedListIterator<rMatrix3d> istress(stress_list);
    LinkedListIterator<rMatrix3d> istrain(strain_list);
    for (int p=0; istress; p++,istress++,istrain++) {
      MultiDimIterator<iVector2d> i(iVector2d(3,3));
      iVector2d &ai=(iVector2d &)i;
      for (; i; i++) {
	big_stress(p*block+ai(0)+ai(1)*3)=(*istress)(ai(0),ai(1));
      }
      for (int b=0; b<elas_basis.get_size(); b++) {
	MultiDimIterator<Array<int> > i(size3333);
	Array<int> &ai=(Array<int> &)i;
	for (; i; i++) {
	  big_eqnmat(p*block+ai(0)+ai(1)*3,b)+=elas_basis(b)(ai)*(*istrain)(ai(2),ai(3));
	}
      }
    }
    Array<Real> beta;
    calc_ols(&beta,big_eqnmat,big_stress);
    	  
    Array<Real> stress_hat;
    product(&stress_hat,big_eqnmat,beta);
	
    {
      ofstream logfile("calcelas.log");
      logfile << "system of equations:";
      logfile << big_eqnmat << endl;

      logfile << "actual stress   predicted stress   difference" << endl;
      for (int i=0; i<big_stress.get_size(); i++) {
	logfile << big_stress(i) << " " << stress_hat(i) << " " << (big_stress(i)-stress_hat(i)) << endl;
      }
      logfile << "end" << endl;
    }

    rTensor elas(size3333);
    elas.zero();
    for (int b=0; b<elas_basis.get_size(); b++) {
      MultiDimIterator<Array<int> > i(size3333);
      for (; i; i++) {
	elas(i)+=-elas_basis(b)(i)*beta(b);
      }
    }

    MultiDimIterator<Array<int> > i(size3333);
    for (; i; i++) {
      cout << elas(i) << endl;
    }
  }
}
