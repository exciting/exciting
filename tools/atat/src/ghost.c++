#include <fstream.h>
#include "xtalutil.h"
#include "parse.h"
#include "getvalue.h"

int main(int argc, char *argv[]) {

  AskStruct options[]={
    {"","Chemical potential via the GHOST particle method, by Axel van de Walle",TITLEVAL,NULL},
    {"-r","Maximum radius for calculating pair correlations",REALVAL,&rmax},
    {"-w","Width for smoothing pair correlation",REALVAL,&width},
    {"-s","Maximum shell for neighbor",INTVAL,&shell},
    {"-cs","Maximum shell for Common Neighbor Analysis",INTVAL,&cshell}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    return 1;
  }
  LinkedList<Structure> md_str;
  while (!infile.eof()) {
    Structure str;
    parse_lattice_file(&(str.cell), &(str.atom_pos), &(str.atom_type), &site_type_list, &atom_label, infile);
    wrap_inside_cell(&(str.atom_pos), str.atom_pos, str.cell);
    fix_atom_type(&str, site_type_list);
    calc_pair_corr(&corr, str, rmax);
    md_str << new Structure(str);
    skip_to_next_structure(infile);
  }

    Array<Real> smocorr;
    Real dr=rmax/((Real)n_bin-1.);
    smooth_pair_corr(&smocorr, corr, width/dr);
    {
      ofstream file("pc.out");
      for (int i=0; i<smocorr.get_size(); i++) {
	file << ((Real)i)*dr << " " << smocorr(i) << endl;
      }
    }
    Array<Real> cutoff;
    find_valleys(&cutoff, smocorr);
    {
      ofstream file("shells.out");
      for (int i=0; i<cutoff.get_size(); i++) {
	file << cutoff(i)*dr << endl;
      }
      file << rmax << endl;
    }
  }
  else {
    shell--;
    cshell--;
    Array<Real> cutoff;
    {
      ifstream file("shells.in");
      if (!file) ERRORQUIT("Unable to open shells.in");
      LinkedList<Real> cutofflist;
      while (!file.eof()) {
	Real c;
	file >> c;
	cutofflist << new Real(c);
	skip_delim(file);
      }
      LinkedList_to_Array(&cutoff,cutofflist);
    }
    Array<Array<int> > site_type_list;
    Array<AutoString> atom_label;
    LinkedList<Array<int> > all_common_neighbor;
    istream &infile=cin;
    while (!infile.eof()) {
      Structure str;
      parse_lattice_file(&(str.cell), &(str.atom_pos), &(str.atom_type), &site_type_list, &atom_label, infile);
      wrap_inside_cell(&(str.atom_pos), str.atom_pos, str.cell);
      fix_atom_type(&str, site_type_list);
      common_neighbor_analysis(&all_common_neighbor, str, cutoff, shell, cshell);
      skip_to_next_structure(infile);
    }
    Real total=0;
    LinkedListIterator<Array<int> > icn(all_common_neighbor);
    for (; icn; icn++) {
      total+=(Real)(*icn)(0);
    }
    icn.init(all_common_neighbor);
    for (; icn; icn++) {
      int i;
      for (i=1; i<icn->get_size()-2; i++) {
	cout << resetiosflags(ios::fixed) << fabs((*icn)(i)) << " ";
      }
      cout << atom_label((*icn)(i)) << " " << atom_label((*icn)(i+1)) << " ";
      cout << setiosflags(ios::fixed) << setprecision(6) << (Real)((*icn)(0))/total << endl;
    }
  }
}
