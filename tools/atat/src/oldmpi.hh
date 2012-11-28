/*
inline int MyMPI_Send(Real *pbuf, int dest, int tag, MPI_Comm comm) {return MPI_Send(pbuf,1,MPI_DOUBLE,dest,tag,comm);}
inline int MyMPI_Send(int  *pbuf, int dest, int tag, MPI_Comm comm) {return MPI_Send(pbuf,1,MPI_INT,   dest,tag,comm);}
inline int MyMPI_Recv(Real *pbuf, int src,  int tag, MPI_Comm comm, MPI_Status *status) {return MPI_Recv(pbuf,1,MPI_DOUBLE,src,tag,comm,status);}
inline int MyMPI_Recv(int  *pbuf, int src,  int tag, MPI_Comm comm, MPI_Status *status) {return MPI_Recv(pbuf,1,MPI_INT,   src,tag,comm,status);}

inline int MyMPI_Send(Array<Real> *pbuf, int dest, int tag, MPI_Comm comm) {return MPI_Send(pbuf->get_buf(),pbuf->get_size(),MPI_DOUBLE,dest,tag,comm);}
inline int MyMPI_Send(Array<int>  *pbuf, int dest, int tag, MPI_Comm comm) {return MPI_Send(pbuf->get_buf(),pbuf->get_size(),MPI_INT,   dest,tag,comm);}

inline int MyMPI_Recv(Array<Real> *pbuf, int src,  int tag, MPI_Comm comm, MPI_Status *status) {
  MPI_Probe(src,tag,comm,status);
  int len;
  MPI_Get_count(status, MPI_DOUBLE, &len);
  pbuf->resize(len);
  return MPI_Recv(pbuf->get_buf(),len,MPI_DOUBLE,src,tag,comm,status);
}

inline int MyMPI_Recv(Array<int>  *pbuf, int src,  int tag, MPI_Comm comm, MPI_Status *status) {
  MPI_Probe(src,tag,comm,status);
  int len;
  MPI_Get_count(status, MPI_INT, &len);
  pbuf->resize(len);
  return MPI_Recv(pbuf->get_buf(),len,MPI_INT,   src,tag,comm,status);
}
*/
