#ifndef __STRINGO_H__
#define __STRINGO_H__

#include <string.h>
#include <iostream.h>
#include <iomanip.h>
#include <strstream.h>
#include "misc.h"

class AutoString {
protected:
	char *buf;
#ifdef DEBUG
  int truelen;
#endif
	static char *empty;
public:
	AutoString(const char *str) {
	  buf=NULL;
	  set(str);
	}
	AutoString(const AutoString &a) {
	  buf=NULL;
	  set(a);
	}
	AutoString(int len) {
	  buf=NULL;
	  set(len);
	}
	void set(const char *str) {
	  if (buf) delete[] buf;
	  if (str) {
	    buf=new char[strlen(str)+1];
	    strcpy(buf,str);
#ifdef DEBUG
	    truelen=strlen(str);
#endif
	  }
	  else
	    buf=NULL;
	}
	void set(int len) {
	  if (buf) delete[] buf;
	  buf=new char[len+1];
	  MEMCLR(buf,len+1);
#ifdef DEBUG
	  truelen=len;
#endif
	}
	AutoString(void) {
	  buf=NULL;
#ifdef DEBUG
          truelen=0;
#endif
	}
	int len(void) const {
	  if (buf) {
	    return strlen(buf);
	  }
	  else {
	    return 0;
	  }
	}
	void operator=(const AutoString &a) {
	  set(a);
	}
#ifndef STRING_FIX
	operator char * () {
	  if (buf)
	    return buf;
	  else
	    return empty;
	}
#endif
	operator const char * () const {
	  if (buf)
	    return buf;
	  else
	    return empty;
	}
	char& operator [](int i) {
#ifdef DEBUG
	  if ( i<0 || i>=truelen ) {
	    cerr << "AutoString out of range: " << i << "/" << truelen << endl;
	    return *empty;
	  }
	  else
#endif
	    return buf[i];
	}
	const char& operator [](int i) const {
#ifdef DEBUG
	  if ( i<0 || i>=truelen ) {
	    cerr << "AutoString out of range: " << i << "/" << truelen << endl;
	    return *empty;
	  }
	  else
#endif
	    return buf[i];
	}
	void operator += (const AutoString &s) {
	  char *newbuf=new char[len()+s.len()+1];
	  strcpy(newbuf,*this);
	  strcpy(newbuf+len(),s);
	  delete[] buf;
	  buf=newbuf;
#ifdef DEBUG
          truelen=strlen(newbuf);
#endif
	}
	void operator += (char c) {
	  char *newbuf=new char[len()+2];
	  strcpy(newbuf,*this);
	  newbuf[len()]=c;
	  newbuf[len()+1]=0;
	  delete[] buf;
	  buf=newbuf;
#ifdef DEBUG
          truelen=strlen(newbuf);
#endif
	}
	int operator == (const AutoString &s) const {
	  return (strcmp(*this,s)==0);
	}
	~AutoString() {
	  if (buf) delete[] buf;
	}
	friend istream& operator >> (istream &file, AutoString &str);
};

Real to_real(const AutoString &s);

#include "binstream.h"

inline ostream & bin_ostream(ostream &file, AutoString &str) {
  bin_ostream(file, str.len());
  for (int i=0; i<str.len(); i++) {
    bin_ostream(file, str[i]);
  }
  return file;
}

inline istream & bin_istream(istream &file, AutoString &str) {
  int len;
  bin_istream(file, len);
  str.set(len);
  for (int i=0; i<len; i++) {
    bin_istream(file, str[i]);
  }
  return file;
}

#endif
