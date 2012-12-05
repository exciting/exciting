#ifndef __PLUGIN_H__
#define __PLUGIN_H__

#include "linklist.h"

template<class B>
class GenericPlugIn {
 protected:
  const char *label;
  GenericPlugIn<B> *next;
  static GenericPlugIn<B> *list;
 public:
  static B *create(const char *_label) {
    GenericPlugIn<B> *i=list;
    while (i && strcmp(i->label,_label)!=0) i=i->next;
    if (!i) return NULL;
    return i->make_new();
  }
  virtual B *make_new(void) {return NULL;}
};

template<class B,class T>
class SpecificPlugIn: public GenericPlugIn<B> {
 public:
  SpecificPlugIn(const char *_label): GenericPlugIn<B>() {
    GenericPlugIn<B>::label=_label;
    GenericPlugIn<B>::next=GenericPlugIn<B>::list;
    GenericPlugIn<B>::list=this;
  }
  virtual B *make_new(void) {
    return new T;
  }
};

template<class B>
void make_plug_in_list(LinkedList<B> *list, const char *labels, char delim) {
  int len=strlen(labels);
  char *lab=new char[len+1];
  strcpy(lab,labels);
  int i=0;
  int j=0;
  while (1) {
    j=i;
    if (j>=len) break;
    while (i<len  && lab[i]!=delim) i++;
    lab[i]=0;
    B *p=GenericPlugIn<B>::create(lab+j);
    if (p) {
      (*list) << p;
    }
    i++;
  }
  delete[] lab;
  return;
}

template<class B>
void make_plug_in_list(LinkedList<B> *list, const char *labels) {
  make_plug_in_list(list,labels,'_');
}

template<class B>
int check_plug_in(const B &dummy, const char *labels, char delim) {
  int len=strlen(labels);
  char *lab=new char[len+1];
  strcpy(lab,labels);
  int isok=1;
  int i=0;
  int j=0;
  while (1) {
    j=i;
    if (j>=len) break;
    while (i<len  && lab[i]!=delim) i++;
    lab[i]=0;
    B *p=GenericPlugIn<B>::create(lab+j);
    if (p) {
      delete p;
    }
    else {
      cerr << "Unable to find plug-in '" << ((const char *)lab)+j << "'." << endl;
      isok=0;
    }
    i++;
  }
  delete[] lab;
  return isok;
}

template<class B>
int check_plug_in(const B &dummy, const char *labels) {
  return check_plug_in(dummy,labels,'_');
}

#endif
