#include <strstream.h>
#include "getvalue.h"
#include "linklist.h"

char *equal_delim=":= \t\n";
char *blanks_delim=" \t\n";

LinkedList<AutoString> getvalue_string_buffer;

int get_values(istream &s, int nb, AskStruct *label) {
  char c;
  while (1) {
    while (1) {
      c=s.get();
      if (!strchr(blanks_delim,c) || s.eof()) break;
    }
    if (s.eof()) return 1;
    s.putback(c);
    ostrstream cur_label;
    while (1) {
      c=s.get();
      if (strchr(equal_delim,c) || s.eof()) break;
      cur_label << c;
    }
    cur_label << '\0';
    while (1) {
      c=s.get();
      if (!strchr(equal_delim,c) || s.eof()) break;
    }
    s.putback(c);
    int i;
    for (i=0; i<nb; i++) {
      if (strcmp(label[i].shortname,cur_label.str())==0) {
	ostrstream value_str;
	if (label[i].vartype!=BOOLVAL) {
	  while (1) {
	    c=s.get();
	    if (strchr(blanks_delim,c) || s.eof()) break;
	    value_str << c;
	  }
	}
	value_str << '\0';
	istrstream value(value_str.str());
	switch (label[i].vartype) {
        case INTVAL:
          value >> (*(int *)(label[i].outvar));
	  break;
        case REALVAL:
          value >> (*(Real *)(label[i].outvar));
	  break;
        case BOOLVAL:
          (*(int *)(label[i].outvar))=1;
	  break;
        case STRINGVAL:
	  {
	    AutoString *pstr=new AutoString(value_str.str());
	    getvalue_string_buffer << pstr;
	    *((const char **)(label[i].outvar))=(const char *)(*pstr);
	  }
	  break;
        case CHOICEVAL:
          value >> (*(char *)(label[i].outvar));
	  break;
	}
	break;
      }
    }
    if (i==nb) return 0;
  }
}

int get_values(int argc, char* argv[], int nb, AskStruct *label) {
  if (argc==1) return 0;
  ostrstream os;
  for (int i=1; i<argc; i++) {
    os << argv[i] << ' ';
  }
  os << '\0';
  istrstream is(os.str());
  return get_values(is, nb, label);
}

void display_help(int nb, AskStruct *label) {
  int maxw=0;
  for (int i=0; i<nb; i++) {
    maxw=MAX(maxw,strlen(label[i].shortname));
  }
  for (int i=0; i<nb; i++) {
    if (label[i].vartype==TITLEVAL) {
      cerr << label[i].longname << endl;
    }
    else {
      cerr.width(maxw);
      cerr << label[i].shortname;
      switch (label[i].vartype){
        case INTVAL:    cerr << "=[int]    ";break;
        case REALVAL:   cerr << "=[real]   ";break;
        case BOOLVAL:   cerr << "          ";break;
        case STRINGVAL: cerr << "=[string] ";break;
        case CHOICEVAL: cerr << "=[choice] ";break;
      }
      cerr << label[i].longname << endl;
    }
  }
}

void chdir_robust(const char *dir) {
  if (chdir(dir)!=0) {
    cerr << "Cannot cd into " << dir << endl;
    ERRORQUIT("Aborting");
  }
}

int get_string(AutoString *ps, istream &file, char *delim) {
  AutoString accum;
  char c;
  while (1) {
     file.get(c);
     if (strchr(delim,c)!=NULL || file.eof()) break;
     accum+=c;
  }
  *ps=accum;
  if (file.eof()) return 0;
  file.putback(c);
  return 1;
}

int skip_delim(istream &file, char *delim) {
  char c;
  do {
     file.get(c);
  } while (strchr(delim,c)!=NULL && !file.eof());
  if (file.eof()) return 0;
  file.putback(c);
  return 1;
}
