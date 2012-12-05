#include "mpiinterf.h"

MyMPIobj_class MyMPIobj;

#ifdef DEBUG
ofstream mpidebugfile;
#endif

/*
template<class T>
class MPISynchronizerStream {
  Array<LinkedList<T *> > to_update_list;
public:
  MPISynchronizerStream(void): to_update_list(MyMPIobj.numproc) {}
  int is_my_job(int index) {
    return (index % MyMPIobj.numproc == MyMPIobj.id);
  }
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
	  int len;
	  MPI_Get_count(&status, MPI_CHAR, &len);
	  char *pbuf=new char[len];
	  MPI_Recv(pbuf,len,MPI_CHAR,i,i,MPI_COMM_WORLD,&status);
	  T **pobj=to_update_list(i).detach(head);
	  istrstream buf(pbuf);
	  buf >> (**pobj);
	  T *tmp=*pobj;
	  delete pobj;
	  delete pbuf;
	}
      }
    }
  }
  void sync(T *pobject, int index) {
    int windex=index % MyMPIobj.numproc;
    if (windex==MyMPIobj.id) {
      ostrstream buf;
      buf << *pobject << '\0';
      char *pbuf=buf.str();
      int len=strlen(pbuf)+1;
      for (int i=0; i<MyMPIobj.numproc; i++) {
	if (i!=MyMPIobj.id) {
	  MPI_Send(pbuf,len,MPI_CHAR,i,MyMPIobj.id,MPI_COMM_WORLD);
	}
      }
    }
    else {
      to_update_list(windex) << new (T *)(pobject);
    }
    finish(0);
  }
  ~MPISynchronizerStream() {
    finish(1);
  }
};
*/

