#include <fstream.h>
#include "xtalutil.h"
#include "parse.h"
#include "getvalue.h"

int find_shell(const Array<Real> &cutoff, Real r) {
  for (int s=0; s<cutoff.get_size(); s++) {
    if (cutoff(s)>r) return s;
  }
  return cutoff.get_size();
}

class FindSpan {
  rVector3d threshold;
  rMatrix3d icell;
public:
  FindSpan(const rMatrix3d &cell, Real rmax): threshold(), icell() {
    icell=!cell;
    rMatrix3d rcell=~icell;
    for (int i=0; i<3; i++) {
      threshold(i)=rmax*norm(rcell.get_column(i))+zero_tolerance;
    }
  }
  void get_span(iVector3d *pminc, iVector3d *pmaxc, const rVector3d pos) {
    rVector3d f=icell*pos;
    for (int i=0; i<3; i++) {
      (*pminc)(i)=min((int)(floor(f(i)-threshold(i))),0);
      (*pmaxc)(i)=max((int)(floor(f(i)+threshold(i))),0);
    }
  }
};

class AtomAndType {
public:
  rVector3d pos;
  int type;
  int index;
  AtomAndType(void): pos() {}
  AtomAndType(const AtomAndType &a): pos(a.pos), type(a.type), index(a.index) {}
  AtomAndType(const rVector3d &_pos, int _type, int _index): pos(_pos), type(_type), index(_index) {}
};

void calc_pair_corr(Array<Real> *pcorr, const Structure &str, Real rmax) {
  Real nb=(Real)(pcorr->get_size()-1)-zero_tolerance;
  FindSpan find_span(str.cell,rmax);
  for (int i=0; i<str.atom_pos.get_size(); i++) {
    iVector3d minc,maxc;
    find_span.get_span(&minc,&maxc,str.atom_pos(i));
    MultiDimIterator<iVector3d> c(minc,maxc);
    for (; c; c++) {
      for (int j=0; j<str.atom_pos.get_size(); j++) {
	rVector3d x=str.atom_pos(j)+str.cell*to_real((iVector3d &)c);
	Real r=norm(x-str.atom_pos(i));
	if (r>zero_tolerance && r<=rmax) {
	  int b=(int)(r*nb/rmax);
	  (*pcorr)(b)+=1./sqr(r);
	}
      }
    }      
  }
}

void smooth_pair_corr(Array<Real> *psmocorr, const Array<Real> &corr, Real width) {
  Real h=2.*width;
  Real k=1./sqrt(2.*M_PI)/width;
  psmocorr->resize(corr.get_size());
  zero_array(psmocorr);
  for (int i=0; i<corr.get_size(); i++) {
    for (int j=0; j<corr.get_size(); j++) {
      (*psmocorr)(i)+=corr(j)*k*exp(-(Real)sqr(i-j)/h);
    }
  }
}

void find_valleys(Array<Real> *pcutoff, const Array<Real> &corr) {
  LinkedList<Real> cutofflist;
  for (int i=1; i<corr.get_size()-1; i++) {
    if (corr(i-1)>corr(i) && corr(i+1)>corr(i)) {
      cutofflist << new Real((Real)i);
    }
  }
  LinkedList_to_Array(pcutoff,cutofflist);
}

int array_order(const Array<int> &a, const Array<int> &b, int start) {
  for (int i=start; i<a.get_size(); i++) {
    if (a(i)<b(i)) return -1;
    if (a(i)>b(i)) return 1;
  }
  return 0;
}

void common_neighbor_analysis(LinkedList<Array<int> > *pall_common_neighbor, const Structure &str, const Array<Real> &cutoff, int shell, int cshell) {
  Array<LinkedList<AtomAndType> > neighbor_list(str.atom_pos.get_size());
  Real rmax=cutoff(shell);
  FindSpan find_span(str.cell,rmax);
  for (int i=0; i<str.atom_pos.get_size(); i++) {
    //cerr << str.atom_pos(i) << endl;    
    iVector3d minc,maxc;
    find_span.get_span(&minc,&maxc,str.atom_pos(i));
    MultiDimIterator<iVector3d> c(minc,maxc);
    for (; c; c++) {
      for (int j=0; j<str.atom_pos.get_size(); j++) {
	rVector3d x=str.atom_pos(j)+str.cell*to_real((iVector3d &)c);
	Real r=norm(x-str.atom_pos(i));
	if (r>zero_tolerance && r<=rmax) {
	  neighbor_list(i) << new AtomAndType(x,str.atom_type(j),j);
	  //cerr << "    " << x << endl;
	}
      }
    }      
  }
  Array<int> done(str.atom_pos.get_size());
  zero_array(&done);
  for (int i=0; i<str.atom_pos.get_size(); i++) {
    done(i)=1;
    LinkedListIterator<AtomAndType> rj(neighbor_list(i));
    for (; rj; rj++) {
      if (!done(rj->index)) {
	LinkedList<AtomAndType> common_neighbor;
	common_neighbor << new AtomAndType(str.atom_pos(i),str.atom_type(i),i);
	//cerr << str.atom_pos(i) << endl;
	common_neighbor << new AtomAndType(*rj);
	//cerr << rj->pos << endl;
	Real crmax=cutoff(cshell);
	LinkedListIterator<AtomAndType> rk(neighbor_list(i));
	for (; rk; rk++) {
	  if (norm(rk->pos - str.atom_pos(i))<=crmax) {
	    Real d=norm(rk->pos - rj->pos);
	    if ( d>zero_tolerance && d<=crmax) {
	      common_neighbor << new AtomAndType(*rk);
	      //cerr << "  " << rk->pos << endl;	  
	    }
	  }
	}
	Array<AtomAndType> cn;
	LinkedList_to_Array(&cn,common_neighbor);
	Array<int> cn_code(1+2+2+cutoff.get_size()+1);
	zero_array(&cn_code);
	cn_code(0)=1;
	cn_code(1)=find_shell(cutoff,norm(cn(0).pos-cn(1).pos))+1;
	cn_code(2)=cn.get_size()-2;
	for (int a=2; a<cn.get_size(); a++) {
	  for (int b=a+1; b<cn.get_size(); b++) {
	    int s=find_shell(cutoff,norm(cn(a).pos-cn(b).pos));
	    cn_code(3+s)--;
	  }
	}
	int k=3+cutoff.get_size()+1;
	cn_code(k)=str.atom_type(i);
	cn_code(k+1)=rj->type;
	if (cn_code(k)<cn_code(k+1)) { swap(&(cn_code(k)),&(cn_code(k+1)));}
	LinkedListIterator<Array<int> > icn(*pall_common_neighbor);
	for (; icn; icn++) {
	  if (array_order(*icn,cn_code,1)!=-1) break;
	}
	if (!icn) {
	  pall_common_neighbor->add(new Array<int>(cn_code),icn);
	} else {
	  if (array_order(*icn,cn_code,1)==0) {
	    (*icn)(0)++;
	  }
	  else {
	    pall_common_neighbor->add(new Array<int>(cn_code),icn);
	  }
	}
      }
    }
  }
}

int main(int argc, char *argv[]) {
  Real rmax=0.;
  int n_bin=500;
  Real width=0.5;
  int shell=1;
  int cshell=1;
  AskStruct options[]={
    {"","Common Neighbor Analysis code, by Axel van de Walle",TITLEVAL,NULL},
    {"","  First run with -r option, copy shells.out to shells.in and then run with -s",TITLEVAL,NULL},
    {"-r","Maximum radius for calculating pair correlations",REALVAL,&rmax},
    {"-w","Width for smoothing pair correlation",REALVAL,&width},
    {"-s","Maximum shell for neighbor",INTVAL,&shell},
    {"-cs","Maximum shell for Common Neighbor Analysis",INTVAL,&cshell}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    return 1;
  }
  if (shell<cshell) ERRORQUIT("Error: must have s>=cs");
  if (rmax>0.) {
    Array<Array<int> > site_type_list;
    Array<AutoString> atom_label;
    istream &infile=cin;
    Array<Real> corr(n_bin);
    zero_array(&corr);
    while (!infile.eof()) {
      Structure str;
      parse_lattice_file(&(str.cell), &(str.atom_pos), &(str.atom_type), &site_type_list, &atom_label, infile);
      wrap_inside_cell(&(str.atom_pos), str.atom_pos, str.cell);
      fix_atom_type(&str, site_type_list);
      calc_pair_corr(&corr, str, rmax);
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
