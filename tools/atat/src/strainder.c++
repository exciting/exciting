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
  char *tensorfilename="tensor.out";
  Real max_strain=0.01;
  int nb_strain=1;
  int dohelp=0;
  int sigdig=5;
  int fit=0;
  int rank=1;
  int nozero=0;
  char *sym_indices="";
  int proper=0;
  int dummy=0;
  AskStruct options[]={
    {"","STRAIN DERivative calculator " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-f","Fit (instead of the default, generate perturbations)",BOOLVAL,&fit},
    {"-s","Input file defining the relaxed structure (default: str_relax.out)",STRINGVAL,&strfilename},
    {"-t","Name of the tensor file",STRINGVAL,&tensorfilename},
    {"-r","Rank of the tensor",INTVAL,&rank},
    {"-si","Symmetric indexes (comma-separated list; indices start at 0)",STRINGVAL,&sym_indices},
    {"-ms","Maximum strain magnitude imposed (default: 0.01)",REALVAL,&max_strain},
    {"-ns","Number of strain values (default: 1)",INTVAL,&nb_strain},
    {"-nz","Do not assume tensor is zero at zero strain",BOOLVAL,&nozero},
    {"-pr","Compute proper piezoelectric constant",BOOLVAL,&proper},
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
	  strainfile.setf(ios::fixed);
	  strainfile.precision(sigdig);

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
    Array<int> tensorsize(rank);
    if (rank==0) {
      tensorsize.resize(1);
      tensorsize(0)=1;
    }
    else {
      fill_array(&tensorsize,3);
    }

    cerr << "Reading..." << endl;

    rTensor atzero(tensorsize);
    atzero.zero();
    if (nozero) {
      ifstream tensorfile(tensorfilename);
      if (!tensorfile) {
	cerr << " Unable to open file " << tensorfilename << " for zero strain" << endl;
	ERRORQUIT("Aborting");
      }
      MultiDimIterator<Array<int> > i(tensorsize);
      for (; i; i++) {
	tensorfile >> atzero(i);
      }
    }

    system("ls p[-+]*/str_relax.out 2> /dev/null | sed 's+/str_relax.out++g' > pertlist.out");
    LinkedList<rTensor> tensor_list;
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
	if (file_exists(tensorfilename)) {
	  rMatrix3d strain=sstr.cell*(!str.cell)-Id;
	  strain_list << new rMatrix3d(strain);
          ifstream tensorfile(tensorfilename);
	  rTensor tensor(tensorsize);
	  MultiDimIterator<Array<int> > i(tensorsize);
	  for (; i; i++) {
	    tensorfile >> tensor(i);
	    tensor(i)-=atzero(i);
	  }
          tensor_list << new rTensor(tensor);
	  cerr << endl;
        }
	else {
	  cerr << " Unable to open file " << tensorfilename << endl;
	}

	chdir_robust("..");
      }
    }
    unlink("pertlist.out");

    Array<Array<int> > symflip;
    {
      //create symflip array
      istrstream idxstring(sym_indices);
      LinkedList<int> indexlist;
      while (1) {
	int num=-1;
	idxstring >> num;
	if (num==-1) break;
	indexlist << new int(num);
      }
      if (indexlist.get_size()%2!=0) {
	ERRORQUIT("Must specify an even number of indices");
      }
      
      symflip.resize(indexlist.get_size()/2+1);
      LinkedListIterator<int> it(indexlist);
      int j=0;
      for (; it; j++) {
	symflip(j).resize(2);
	symflip(j)(0)=*it;
	it++;
	symflip(j)(1)=*it;
	it++;
      }
      symflip(j).resize(2);
      symflip(0)(0)=rank+0;
      symflip(0)(1)=rank+1;
      // end
    }

    Array<rTensor> tensor_basis;
    find_symmetric_basis(&tensor_basis, rank+2,pointgroup, symflip);
    
    int block=ipow(3,rank);
    Array<Real> big_tensor(tensor_list.get_size()*block);
    Array2d<Real> big_eqnmat(tensor_list.get_size()*block,tensor_basis.get_size());
    zero_array(&big_tensor);
    zero_array(&big_eqnmat);
    
    Array<int> resptensorsize(rank+2);
    fill_array(&resptensorsize,3);
    LinkedListIterator<rTensor> itensor(tensor_list);
    LinkedListIterator<rMatrix3d> istrain(strain_list);
    for (int p=0; itensor; p++,itensor++,istrain++) {
      MultiDimIterator<Array<int> > i(tensorsize);
      for (int j=0; i; i++,j++) {
	big_tensor(p*block+j)=(*itensor)((Array<int>&)i);
      }
      for (int b=0; b<tensor_basis.get_size(); b++) {
	MultiDimIterator<Array<int> > i(resptensorsize);
	Array<int> &ai=(Array<int> &)i;
	for (int j=0; i; i++,j++) {
	  big_eqnmat(p*block+(j%block),b)+=tensor_basis(b)(ai)*(*istrain)(ai(rank+0),ai(rank+1));
	}
      }
    }
    Array<Real> beta;
    calc_ols(&beta,big_eqnmat,big_tensor);
    	  
    Array<Real> tensor_hat;
    product(&tensor_hat,big_eqnmat,beta);

    {
      ofstream logfile("strainder.log");
      logfile.setf(ios::fixed);
      logfile.precision(sigdig);
      logfile << "system of equations:";
      logfile << big_eqnmat << endl;

      logfile << "actual tensor   predicted tensor   difference" << endl;
      for (int i=0; i<big_tensor.get_size(); i++) {
	logfile << big_tensor(i) << " " << tensor_hat(i) << " " << (big_tensor(i)-tensor_hat(i)) << endl;
      }
      logfile << "end" << endl;
    }

    rTensor strainder(resptensorsize);
    strainder.zero();
    for (int b=0; b<tensor_basis.get_size(); b++) {
      MultiDimIterator<Array<int> > i(resptensorsize);
      for (; i; i++) {
	strainder(i)+=tensor_basis(b)(i)*beta(b);
      }
    }

    if (proper) {
      if (rank!=1) {ERRORQUIT("rank must be 1 for this option.");}
      MultiDimIterator<Array<int> > i(resptensorsize);
      Array<int> &ai=(Array<int> &)i;
      for (; i; i++) {
	strainder(i)+=( ai(1)==ai(2) ? 1. : 0.)*atzero.vectorize()(ai(0))-( ai(0)==ai(1) ? 1. : 0.)*atzero.vectorize()(ai(2));
      }
    }


    cout.setf(ios::fixed);
    cout.precision(sigdig);
    MultiDimIterator<Array<int> > i(resptensorsize);
    for (; i; i++) {
      cout << strainder(i) << endl;
    }
  }
}
