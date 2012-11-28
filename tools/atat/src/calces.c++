#include "getvalue.h"
#include "parse.h"
#include "phonlib.h"
#include "version.h"
#include <fstream.h>

Real calc_ewald(const Structure &str, const Array<Real> &charges, Real eta=-1.0, Real prec=1e-25) {
    rMatrix3d rcell=2.*M_PI*!(~str.cell);
    LatticePointIterator lat(str.cell,1);
    Real rmin=norm((rVector3d)lat);
    LatticePointIterator rlat(rcell,1);
    Real kmin=norm((rVector3d)rlat);
    Real omega=fabs(det(str.cell));
    if (eta<0) {
      eta=sqrt(0.5*(kmin/rmin)*sqrt(log(prec)/log((prec*omega)/(4*M_PI))));
      //cerr << eta << endl;
    }
    Real E=0.;
    for (int i=0; i<str.atom_pos.get_size(); i++) {
	for (int j=0; j<=i; j++) {
	    Real qq=charges(i)*charges(j);
	    Real doublecount=(j==i ? 0.5 : 1.);
	    rVector3d dr=str.atom_pos(i)-str.atom_pos(j);
	    lat.init(str.cell,(j==i));
	    Real maxr=sqrt(-log(prec))/eta+norm(dr);
	    int nt=0;
	    while (1) {
		Real r=norm((rVector3d)lat+dr);
		if (r>maxr) break;
		E+=qq*doublecount*erfc(eta*r)/r;
		lat++;
		nt++;
	    }
	    //cerr << "nReal=" << nt << endl;
	    rlat.init(rcell,1);
	    Real maxk2=4.*sqr(eta)*(-log(omega*prec/(4.*M_PI)));
	    nt=0;
	    while (1) {
		rVector3d k=(rVector3d)rlat;
		Real k2=norm2(k);
		if (k2>maxk2) break;
		E+=qq*cos(k*dr)*doublecount*exp(-k2/(4.*sqr(eta)))*4.*M_PI/omega/k2;
		rlat++;
		nt++;
	    }
	    //cerr << "nReciprocal=" << nt << endl;
	    E-=doublecount*qq*M_PI/sqr(eta)/omega;
	    if (i==j) {
		E-=qq*eta/sqrt(M_PI);
	    }
	}
    }
    return E*14.39964886; // Coulomb constant in eV*A/e^2;
}

int main(int argc, char *argv[]) {
  char *unrelfilename="str.out";
  char *chargefilename="charges.in";
  Real eta=-1.0;
  Real prec=1e-25;
  int defval=0;
  AskStruct options[]={
    {"-s","structure file (default: str.out)",STRINGVAL,&unrelfilename},
    {"-c","Charge file (default: charges.in)",STRINGVAL,&chargefilename},
    {"-e","eta",REALVAL,&eta},
    {"-p","prec",REALVAL,&prec},
    {"-d","Use all defaults",BOOLVAL,&defval}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    return 1;
  }

  Array<AutoString> species_name;
  Array<Real> species_charge;
  {
    ifstream chargefile(chargefilename);
    if (!chargefile) ERRORQUIT("Unable to open charges.in file.");
    LinkedList<AutoString> specieslist;
    LinkedList<Real> chargelist;
    while (! chargefile.eof()) {
	AutoString specie;
	int charge;
	skip_delim(chargefile," \t=\n");
	if (chargefile.eof()) break;
	get_string(&specie,chargefile," \t=");
	skip_delim(chargefile," \t=");
	chargefile >> charge;
	specieslist << new AutoString(specie);
	chargelist << new Real(charge);
    }
    LinkedList_to_Array(&species_name,specieslist);
    LinkedList_to_Array(&species_charge,chargelist);
  }
  //cerr << species_name << endl;
  //cerr << species_charge << endl;

  ifstream file(unrelfilename);
  if (!file) ERRORQUIT("Unable to open unrelaxed file.");

  while (1) {
    if (file.eof()) break;
    Structure ustr;
    Array<Arrayint> site_type_list;
    Array<AutoString> atom_label;
    rMatrix3d axes;
    parse_lattice_file(&ustr.cell, &ustr.atom_pos, &ustr.atom_type, &site_type_list, &atom_label, file, &axes);
    skip_to_next_structure(file);
    Array<Real> charges(ustr.atom_pos.get_size());
    for (int i=0; i<ustr.atom_pos.get_size(); i++) {
      int j=index_in_array(species_name,atom_label(site_type_list(ustr.atom_type(i))(0)));
      if (j==-1) {
	cerr << "Unable to find specie" << atom_label(site_type_list(ustr.atom_type(i))(0)) << endl;
	ERRORQUIT("Aborting");
      }
      charges(i)=species_charge(j);
    }
    //cerr << charges << endl;
    Real E=calc_ewald(ustr,charges,eta,prec);
    cout << E << endl;
  }
}
