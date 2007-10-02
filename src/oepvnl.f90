
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine oepvnl(vnlcv,vnlvv)
use modmain
use modmpi
implicit none
! arguments
complex(8), intent(out) :: vnlcv(ncrmax,natmtot,nstsv,nkpt)
complex(8), intent(out) :: vnlvv(nstsv,nstsv,nkpt)
! local variables
integer ik,i
 integer::mpireccnts(procs),mpirecdispls(procs),mpisndcnts(procs),&
       mpisnddispls(procs),kgatherdispls(procs),kgatherrecvcnts(procs)
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
#ifdef MPIEXX
do ik=firstk(rank),lastk(rank)
#endif
#ifndef MPIEXX
do ik=1,nkpt
#endif
!$OMP CRITICAL
  write(*,'("Info(oepvnl): ",I6," of ",I6," k-points on proc:",I6)') ik,nkpt,rank
!$OMP END CRITICAL
  call oepvnlk(ik,vnlcv(1,1,1,ik),vnlvv(1,1,ik))
end do
!$OMP END DO
!$OMP END PARALLEL

#ifdef MPIEXX
  do i=0,procs-1
     kgatherdispls(i+1)=firstk(i)-1
     kgatherrecvcnts(i+1)=nofk(i)
  end do
  mpireccnts=kgatherrecvcnts*ncrmax*natmtot*nstsv
  mpirecdispls=kgatherdispls*ncrmax*natmtot*nstsv
  mpisndcnts(:)=nofk(rank)*ncrmax*natmtot*nstsv
  mpisnddispls(:)=(firstk(rank)-1)*ncrmax*natmtot*nstsv
  call MPI_Alltoallv(vnlcv, mpisndcnts,  mpisnddispls, MPI_DOUBLE_COMPLEX, vnlcv,&
       mpireccnts,  mpirecdispls,  MPI_DOUBLE_COMPLEX ,  MPI_COMM_WORLD,ierr)
  mpireccnts=kgatherrecvcnts*nstsv*nstsv
  mpirecdispls=kgatherdispls*nstsv*nstsv
  mpisndcnts(:)= nofk(rank)*nstsv*nstsv
  mpisnddispls(:)=(firstk(rank)-1)*nstsv*nstsv
  call MPI_Alltoallv(vnlvv,mpisndcnts,  mpisnddispls, MPI_DOUBLE_COMPLEX, vnlvv,&
       mpireccnts,  mpirecdispls, MPI_DOUBLE_COMPLEX , MPI_COMM_WORLD,ierr)
#endif

return
end subroutine

