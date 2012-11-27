#include <fstream.h>
#include "parse.h"

int is_eol(istream &file) {
  char c;
  do {
    if (file.eof()) return(1);
    file.get(c);
  } while (c==' ' || c=='\t');
  file.putback(c);
  return (c=='\n');
}

void read_cell(rMatrix3d *pcell, istream &file) {
  rVector3d tmp;
  file >> tmp;
  char c;
  do {
    file.get(c);
  } while (c==' ' || c=='\t');
  file.putback(c);
  if (c!='\n') {
    rVector3d length,angle;
    length=tmp;
    file >> angle;
    calc_lattice_vectors(pcell,length,angle);
  }
  else {
    pcell->set_column(0,tmp);
    file >> tmp;
    pcell->set_column(1,tmp);
    file >> tmp;
    pcell->set_column(2,tmp);
  }
}

void parse_lattice_file(rMatrix3d *pcell, Array<rVector3d> *patom_pos, Array<int> *psite_type, Array<Arrayint> *psite_type_list, Array<AutoString> *patom_label, istream &file, rMatrix3d *paxes) {
  char *delim=" \t,;/";
  // read cell;
  rMatrix3d axes;
  read_cell(&axes,file);
  for (int i=0; i<3; i++) {
    rVector3d v;
    file >> v;
    pcell->set_column(i,axes*v);
  }
  if (paxes) {*paxes=axes;}

  // read in atoms;
  LinkedList<rVector3d> atom_pos_list;
  LinkedList<ArrayAutoString> long_site_type_list;
  while (1) {
    rVector3d pos(MAXFLOAT,MAXFLOAT,MAXFLOAT);
    file >> pos;
    if (file.eof() || pos(0)==MAXFLOAT || !file) break;
    pos=axes*pos;
    char buf[MAX_LINE_LEN];
    file.get(buf,MAX_LINE_LEN-1);
    LinkedList<AutoString> site_types;
    char *atom_begin=buf;
    char *string_end=buf+strlen(buf);
    while (atom_begin<string_end) {
      while (strchr(delim,*atom_begin) && atom_begin<string_end) {atom_begin++;}
      if (atom_begin==string_end) break;
      char *atom_end=atom_begin;
      while (!strchr(delim,*atom_end) && atom_begin<string_end) {atom_end++;}
      *atom_end=0;
/*
      LinkedListIterator<AutoString> i(site_types);
      while (i && strcmp(*i,atom_begin)<0) {i++;}
      site_types.add(new AutoString(atom_begin),i);
*/
      site_types << new AutoString(atom_begin);
      atom_begin=atom_end+1;
    }
    Array<AutoString> *parraystring=new Array<AutoString>();
    LinkedList_to_Array(parraystring,site_types);

    LinkedListIterator<ArrayAutoString> i_atom(long_site_type_list);
    LinkedListIterator<rVector3d> i_pos(atom_pos_list);
    while (i_atom && i_atom->get_size() >= parraystring->get_size()) {
      i_atom++;
      i_pos++;
    }
    long_site_type_list.add(parraystring,i_atom);
    atom_pos_list.add(new rVector3d(pos),i_pos);
  }

  // make atom label list;
  LinkedList<AutoString> atom_label_list;
  LinkedListIterator<ArrayAutoString> i_atom(long_site_type_list);
  for (; i_atom; i_atom++) {
    for (int at=0; at<i_atom->get_size(); at++) {
      LinkedListIterator<AutoString> i_label(atom_label_list);
      while (i_label && strcmp(*i_label,(*i_atom)(at))<0) {i_label++;}
      if (!i_label || strcmp(*i_label,(*i_atom)(at))!=0) {
        atom_label_list.add(new AutoString((*i_atom)(at)),i_label);
      }
    }
  }

  // convert string to int in the long site type list;
  LinkedListIterator<ArrayAutoString> il(long_site_type_list);
  LinkedList<Arrayint> long_int_site_type_list;
  for ( ; il; il++) {
    Array<int> *ptype_list=new Array<int>(il->get_size());
    for (int at=0; at<ptype_list->get_size(); at++) {
      LinkedListIterator<AutoString> i_label(atom_label_list);
      int which_type=0;
      while (i_label && strcmp(*i_label,(*il)(at))!=0) {i_label++; which_type++;}
      (*ptype_list)(at)=which_type;
    }
    long_int_site_type_list << ptype_list;
  }

  // create short site type list by removing duplicates;
  LinkedList<Arrayint> short_int_site_type_list;
  LinkedListIterator<Arrayint> ili(long_int_site_type_list);
  for ( ; ili; ili++) {
    LinkedListIterator<Arrayint> is(short_int_site_type_list);
    while (is && (*is)!=(*ili)) {is++;}
    if (!is) {
      short_int_site_type_list << new Array<int>(*ili);
    }
  }

  // convert site type to int;
  LinkedList<int> int_site_type;
  ili.init(long_int_site_type_list);
  for ( ; ili; ili++) {
    LinkedListIterator<Arrayint> is(short_int_site_type_list);
    int which_type=0;
    while ((*is)!=(*ili)) {is++; which_type++;}
    int_site_type << new int(which_type);
  }

  // make arrays from linked lists;
  LinkedList_to_Array(patom_pos,atom_pos_list);
  LinkedList_to_Array(psite_type,int_site_type);
  LinkedList_to_Array(psite_type_list,short_int_site_type_list);
  LinkedList_to_Array(patom_label,atom_label_list);
}

int parse_structure_file(rMatrix3d *pcell, Array<rVector3d> *patom_pos, Array<int> *patom_type, const Array<AutoString> &atom_label, istream &file, rMatrix3d *paxes) {
  // read cell;
  rMatrix3d axes;
  read_cell(&axes,file);
  for (int i=0; i<3; i++) {
    rVector3d v;
    file >> v;
    pcell->set_column(i,axes*v);
  }
  if (paxes) {*paxes=axes;}

  // read in atoms;
  LinkedList<rVector3d> atom_pos_list;
  LinkedList<int> atom_type_list;
  while (1) {
    rVector3d pos(MAXFLOAT,MAXFLOAT,MAXFLOAT);
    file >> pos;
    if (file.eof() || pos(1)==MAXFLOAT || !file) break;
    pos=axes*pos;
    atom_pos_list << new rVector3d(pos);
    char c;
    do {
      file.get(c);
    } while (strchr(" \t",c));
    file.putback(c);
    char buf[MAX_LINE_LEN];
    file.get(buf,MAX_LINE_LEN-1);
    char *end_label=buf;
    while (end_label<buf+strlen(buf) && !strchr(" \t",*end_label)) {end_label++;}
    *end_label=0;
    int i;
    for (i=0; i<atom_label.get_size(); i++) {
      if (strcmp(buf,atom_label(i))==0) {
	atom_type_list << new int(i);
	break;
      }
    }
    if (i==atom_label.get_size()) {
      cerr << "Unknown atom label: " << buf << endl;
      return 0;
    }
  }
  // make arrays from linked lists;
  LinkedList_to_Array(patom_pos,atom_pos_list);
  LinkedList_to_Array(patom_type,atom_type_list);
  return 1;
}

int fix_atom_type(Structure *pstr, const Structure &lat, const Array<Arrayint> &site_type_list, int drop_fixed_site) {
  for (int i=0; i<pstr->atom_type.get_size(); i++) {
    int site_index=which_atom(lat.atom_pos,pstr->atom_pos(i),!lat.cell);
    if (site_index==-1) return 0;
    int site_type=lat.atom_type(site_index);
    int j;
    for (j=0; j<site_type_list(site_type).get_size(); j++) {
      if (pstr->atom_type(i)==site_type_list(site_type)(j)) {
        pstr->atom_type(i)=j;
        break;
      }
    }
    if (j==site_type_list(site_type).get_size()) return 0;
  }
  if (drop_fixed_site) {
    LinkedList<rVector3d> pos_list;
    LinkedList<int> type_list;
    for (int i=0; i<pstr->atom_pos.get_size(); i++) {
      if (site_type_list(lat.atom_type(which_atom(lat.atom_pos,pstr->atom_pos(i),!lat.cell))).get_size()>1) {
	pos_list << new rVector3d(pstr->atom_pos(i));
	type_list << new int(pstr->atom_type(i));
      }
    }
    LinkedList_to_Array(&pstr->atom_pos,pos_list);
    LinkedList_to_Array(&pstr->atom_type,type_list);
  }
  return 1;
}

int fix_atom_type(Structure *pstr, const Array<Arrayint> &site_type_list) {
  int bad=0;
  for (int i=0; i<pstr->atom_type.get_size(); i++) {
    if (site_type_list(pstr->atom_type(i)).get_size()!=1) {bad=1;}
    pstr->atom_type(i)=site_type_list(pstr->atom_type(i))(0);
  }
  return bad;
}

void write_axes(const rMatrix3d &axes, ostream &file, int doabc) {
  if (doabc) {
    rVector3d l,ang;
    calc_lattice_vectors(&l,&ang,axes);
    file << l << " " << ang << endl;
  }
  else {
    for (int i=0; i<3; i++) {
      file << axes.get_column(i) << endl;
    }
  }
}

void write_structure(const Structure &str, const Structure &lat,
		     const Array<Arrayint> &site_type_list, const Array<AutoString> &atom_label, const rMatrix3d &axes, ostream &file, int doabc) {
  rMatrix3d iaxes=!axes;
  write_axes(axes,file,doabc);
  rMatrix3d frac_cell=iaxes*str.cell;
  for (int i=0; i<3; i++) {
    file << frac_cell.get_column(i) << endl;
  }
  for (int i=0; i<str.atom_pos.get_size(); i++) {
    int site_in_cell=which_atom(lat.atom_pos,str.atom_pos(i),!lat.cell);
    file << (iaxes*str.atom_pos(i)) << " " << atom_label(site_type_list(lat.atom_type(site_in_cell))(str.atom_type(i))) << endl;
  }
}

void write_structure(const Structure &str, const Array<AutoString> &atom_label, const rMatrix3d &axes, ostream &file, int doabc) {
  rMatrix3d iaxes=!axes;
  write_axes(axes,file,doabc);
  rMatrix3d frac_cell=iaxes*str.cell;
  for (int i=0; i<3; i++) {
    file << frac_cell.get_column(i) << endl;
  }
  for (int i=0; i<str.atom_pos.get_size(); i++) {
    file << (iaxes*str.atom_pos(i)) << " " << atom_label(str.atom_type(i)) << endl;
  }
}

void skip_to_next_structure(istream &strfile) {
  strfile.clear();
  char c;
  while (!isdigit(c=strfile.get()) && !strfile.eof());
  strfile.putback(c);
}

void write_lattice(const Structure &lat,
		     const Array<Arrayint> &site_type_list, const Array<AutoString> &atom_label, const rMatrix3d &axes, ostream &file, int doabc) {
  rMatrix3d iaxes=!axes;
  write_axes(axes,file,doabc);
  rMatrix3d frac_cell=iaxes*lat.cell;
  for (int i=0; i<3; i++) {
    file << frac_cell.get_column(i) << endl;
  }
  for (int i=0; i<lat.atom_pos.get_size(); i++) {
    file << (iaxes*lat.atom_pos(i)) << " ";
    for (int j=0; j<site_type_list(lat.atom_type(i)).get_size(); j++) {
      if (j!=0) {
	file  << ",";
      }
      file << atom_label(site_type_list(lat.atom_type(i))(j));
    }
    file << endl;
  }
}

Real read_clusters_and_eci(LinkedList<Cluster> *clusterlist, LinkedList<Real> *ecilist,
                           istream &clusterfile, istream &ecifile, const rMatrix3d &axes) {
   Real maxlen=0.;
    while (1) {
      Real mult;
      clusterfile >> mult;
      if (clusterfile.eof()) break;
      Real len;
      clusterfile >> len;
      int nbpt;
      clusterfile >> nbpt;
      clusterfile.get();
      Cluster *pcluster=new Cluster(nbpt);
      for (int j=0; j<nbpt; j++) {
        char buf[MAX_LINE_LEN];
	buf[0]=0;
	while (strlen(buf)==0) {
	  clusterfile.get(buf,MAX_LINE_LEN-1);
	  clusterfile.get();
	}
	  istrstream line(buf);
        line >> (*pcluster)(j);
        (*pcluster)(j)=axes*((*pcluster)(j));
      }
      (*clusterlist) << pcluster;
      maxlen=MAX(maxlen,get_cluster_length(*pcluster));
      if (ecilist) {
	Real cur_eci;
	ecifile >> cur_eci;
	(*ecilist) << new Real(cur_eci);
      }
    }
  return maxlen;
}

Real read_clusters_and_eci(LinkedList<MultiCluster> *clusterlist, LinkedList<Real> *ecilist,
                           istream &clusterfile, istream &ecifile, const rMatrix3d &axes) {
   Real maxlen=0.;
    while (1) {
      Real mult;
      clusterfile >> mult;
      if (clusterfile.eof()) break;
      Real len;
      clusterfile >> len;
      int nbpt;
      clusterfile >> nbpt;
      clusterfile.get();
      MultiCluster *pcluster=new MultiCluster(nbpt);
      for (int j=0; j<nbpt; j++) {
        char buf[MAX_LINE_LEN];
	buf[0]=0;
	while (strlen(buf)==0) {
	  clusterfile.get(buf,MAX_LINE_LEN-1);
	  clusterfile.get();
	}
	istrstream line(buf);
        line >> pcluster->clus(j);
        pcluster->clus(j)=axes*(pcluster->clus(j));
        line >> pcluster->site_type(j) >> pcluster->func(j);
	  if (!line) {
          pcluster->site_type(j)=0;
          pcluster->func(j)=0;
	  }
      }
      (*clusterlist) << pcluster;
      maxlen=MAX(maxlen,get_cluster_length(pcluster->clus));
      if (ecilist) {
	Real cur_eci;
	ecifile >> cur_eci;
	(*ecilist) << new Real(cur_eci);
      }
    }
  return maxlen;
}

int file_exists(const char *filename) {
  ifstream file(filename);
  if (file) {return 1;} else {return 0;}
}
