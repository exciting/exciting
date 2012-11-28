#ifndef __MPIINTERF_H__
#define __MPIINTERF_H__

#include <strstream.h>
#include <fstream.h>
#include "arraylist.h"
#include "binstream.h"

#ifdef ATAT_MPI
#include <mpi.h>

#ifdef DEBUG
extern ofstream mpidebugfile;
#endif

class MyMPIobj_class {
public:
  int numproc;
  int id;
public:
  MyMPIobj_class(void) {}
  void init(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
#ifdef DEBUG
    ostrstream debugfilename;
    debugfilename << "debug" << id << ".out" << '\0';
    mpidebugfile.open(debugfilename.str());
#endif
  }
  int is_root(void) { return id==0; }
  void barrier(void) { MPI_Barrier(MPI_COMM_WORLD); }
  int file_exists(const char *filename) {
    int ex;
    if (is_root()) {
      ifstream file(filename);
      ex=(file ? 1 : 0);
    }
    MPI_Bcast(&ex,1,MPI_INT,0,MPI_COMM_WORLD);
    return ex;
  }
  ~MyMPIobj_class(void) {
    MPI_Finalize();
  }
};

extern MyMPIobj_class MyMPIobj;

template<class T>
inline int MyMPI_SendStream(T *pbuf, int dest, int tag, MPI_Comm comm=MPI_COMM_WORLD) {
    ostrstream line;
    bin_ostream(line,*pbuf);
    return MPI_Send((void *)line.str(),line.tellp(),MPI_CHAR,dest,tag,comm);
}

template<class T>
inline int MyMPI_RecvStream(T *pbuf, int src,  int tag, MPI_Comm comm=MPI_COMM_WORLD, MPI_Status *status=NULL) {
  int len;
  MPI_Status tmpstat;
  if (!status) {status=&tmpstat;}
  MPI_Get_count(status, MPI_CHAR, &len);
  char *ptmpbuf=new char[len];
  int retcode=MPI_Recv((void *)ptmpbuf,len,MPI_CHAR,src,tag,comm,status);
  strstream line;
  line.write(ptmpbuf,len);
  bin_istream(line,*pbuf);
  delete ptmpbuf;
  return retcode;
}

template<class T>
inline int MyMPI_BcastStream(T *pbuf, int root=0, MPI_Comm comm=MPI_COMM_WORLD) {
  int retcode,len;
  if (MyMPIobj.id==root) {
    ostrstream line;
    bin_ostream(line,*pbuf);
    len=line.tellp();
    MPI_Bcast(&len,1,MPI_INT,root,comm);
    retcode=MPI_Bcast((void *)line.str(),line.tellp(),MPI_CHAR,root,comm);
  }
  else {
    MPI_Bcast(&len,1,MPI_INT,root,comm);
    char *ptmpbuf=new char[len];
    retcode=MPI_Bcast((void *)ptmpbuf,len,MPI_CHAR,root,comm);
    strstream line;
    line.write(ptmpbuf,len);
    bin_istream(line,*pbuf);
    delete ptmpbuf;
  }
  return retcode;
}

template<class T>
class MPISynchronizer {
  int index;
  Array<T *> to_update_list;
  void finish(void) {
    int windex=index % MyMPIobj.numproc;
    for (int i=0; i<=windex; i++) {
	MyMPI_BcastStream(to_update_list(i),i);
    }
  }
public:
  MPISynchronizer(void): to_update_list(MyMPIobj.numproc) {index=0;}
  int is_my_job(void) {
    return (index % MyMPIobj.numproc == MyMPIobj.id);
  }
  void sync(T *pobject) {
    to_update_list(index % MyMPIobj.numproc)=pobject;
    if ((index % MyMPIobj.numproc) == MyMPIobj.numproc-1) {
	finish();
    }
    index++;
  }
  ~MPISynchronizer() {
    if ( (index % MyMPIobj.numproc) != 0 ) {
      index--;
      finish();
    }
    MyMPIobj.barrier();
  }
};

/*
template<class T>
class MPISynchronizer {
  int index;
  Array<LinkedList<T *> > to_update_list;
  void finish(int wait=1) {
    for (int i=0; i<MyMPIobj.numproc; i++) {
      if (i!=MyMPIobj.id) {
	while (1) {
	  LinkedListIterator<T *> head(to_update_list(i));
	  if (!head) break;
	  int flag;
	  MPI_Status status;
	  if (wait) {
	    MPI_Probe(i,i,MPI_COMM_WORLD,&status);
	  }
	  else {
	    MPI_Iprobe(i,i,MPI_COMM_WORLD,&flag,&status);
	    if (!flag) break;
	  }
	  T **pobj=to_update_list(i).detach(head);
	  MyMPI_RecvStream(*pobj,i,i,MPI_COMM_WORLD,&status);
	  delete pobj;
	}
      }
    }
  }
public:
  MPISynchronizer(void): to_update_list(MyMPIobj.numproc) {index=0;}
  int is_my_job(void) {
    return (index % MyMPIobj.numproc == MyMPIobj.id);
  }
  void sync(T *pobject) {
    finish(0);
    int windex=index % MyMPIobj.numproc;
    if (windex==MyMPIobj.id) {
      for (int i=0; i<MyMPIobj.numproc; i++) {
	if (i!=MyMPIobj.id) {
	  MyMPI_SendStream(pobject,i,MyMPIobj.id,MPI_COMM_WORLD);
	}
      }
    }
    else {
      to_update_list(windex) << new (T *)(pobject);
      //      cout << "ID=" << MyMPIobj.id << " " << nb_pending << endl;
    }
    finish(0);
    index++;
  }
  ~MPISynchronizer() {
    finish(1);
    MyMPIobj.barrier();
  }
};
*/

template<class T>
void MyMPI_Reduce(T *pobj, void combine(T *, const T &, const T &), MPI_Comm comm=MPI_COMM_WORLD ) {
  MPI_Status status;
  if (!MyMPIobj.is_root()) {
    MyMPI_SendStream(pobj,0,MyMPIobj.id,comm);
  }
  else {
    T a;
    for (int i=1; i<MyMPIobj.numproc; i++) {
      MyMPI_RecvStream(&a,i,i,comm,&status);
      combine(pobj,*pobj,a);
    }
  }
  MyMPI_BcastStream(pobj,0,comm);
}

#undef ERRORQUIT
#define ERRORQUIT(s) {if (MyMPIobj.is_root()) {cerr << s << endl;} exit(1);}

#else

class MyMPIobj_class {
public:
  int numproc;
  int id;
public:
    MyMPIobj_class(void) {}
    void init(int argc, char **argv) {numproc=1; id=0;}
    int is_root(void) { return 1; }
    void barrier(void) {}
    int file_exists(const char *filename) {
      ifstream file(filename);
      return (file ? 1 : 0);
    }
};

extern MyMPIobj_class MyMPIobj;

template<class T>
inline int MyMPI_BcastStream(T *pbuf, int root=0) {
  return 1;
}

template<class T>
class MPISynchronizer {
public:
    MPISynchronizer(void) {}
    int is_my_job(void) {return (1);}
    void sync(T *pobject) {}
};

template<class T>
void MyMPI_Reduce(T *pobj, void combine(T *, const T &, const T &)) {}

#endif

#endif
