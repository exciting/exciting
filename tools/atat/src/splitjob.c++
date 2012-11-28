#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <fstream.h>
#include "misc.h"

main(int argc, char **argv)
{
    int numproc, my_id;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

    char buf[MAX_LINE_LEN];
    char c;
    ifstream file(argv[1]);
    for (int i=0; file; i++) {
	file.get(buf,MAX_LINE_LEN-1);
	file.get(c);
	if (strlen(buf)==0) break;
	if ((i%numproc)==my_id) {
	  cout << "splitjob: Process " << my_id+1 << " of " << numproc << " is running " << buf << endl;
	  system(buf);
	}
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}
