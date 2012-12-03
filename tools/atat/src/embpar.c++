#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <fstream.h>
#include "misc.h"
#include "array.h"

main(int argc, char **argv)
{
    int numproc, my_id;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

    Array<int> done(numproc);
    zero_array(&done);

    for (int i=0; i<my_id; i++) {
      robust_access(&done,i)=1;
    }      

    char buf[MAX_LINE_LEN];
    char c;
    int l=0;
    ifstream file(argv[1]);
    while (1) {
      int flag;
      MPI_Status status;
      while (1) {
	MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);
	if (!flag) {
	  break;
	}
	int source = status.MPI_SOURCE;
	int ll;
	MPI_Recv(&ll,1,MPI_INT,source,1,MPI_COMM_WORLD,&status);
	robust_access(&done,ll)=1;
//	cout << "proc " << my_id << " recv from " << source << " number " << ll << endl;
      }
     
      file.get(buf,MAX_LINE_LEN-1);
      if (strlen(buf)==0) {break;}
      file.get(c);
      if (robust_access(&done,l)==0) {
	for (int j=0; j<numproc; j++) {
	  if (j!=my_id) {MPI_Send(&l,1,MPI_INT,j,1,MPI_COMM_WORLD);
//	  cout << "proc " << my_id << " send to " << j << " number " << l << endl;
	  }
	}
	cout << "splitjob: Process " << my_id+1 << " of " << numproc << " is running " << buf << endl;
	system(buf);
	robust_access(&done,l)=1;
      }
      if (file.eof()) {break;}
      l++;
    }
    cout << "Process " << my_id+1 << " done." << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    cout << "All done." << endl;
    MPI_Finalize();
}
