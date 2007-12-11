
! Copyright (C) 2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine findocclims(iq,iocc0,iocc,iunocc0,iunocc,io0,io,iu0,iu)
  use modmain
  use modxs
  use m_genfilname
  ! arguments
  integer, intent(in) :: iq
  integer, intent(out) :: iocc0,iocc,iunocc0,iunocc
  integer, intent(out) :: io0(nkpt),io(nkpt),iu0(nkpt),iu(nkpt)
  ! local variables
  integer :: ik,ikq,i0,i
  do ik=1,nkpt
     ! k-point set (q=0)
     ! read eigenvectors, eigenvalues and occupancies for G+k (q=0)
     call genfilname(iq=0,setfilext=.true.)
     call getoccsv0(vkl0(1,ik),occsv0(1,ik))
     ! change back file extension
     call genfilname(iq=iq,setfilext=.true.)
     do i0=1,nstsv
        if (occsv0(i0,ik).lt.(1.d0-epsocc)) exit
     end do
     io0(ik)=i0-1
     do i0=nstsv,1,-1
        if (occsv0(i0,ik).gt.epsocc) exit
     end do
     iu0(ik)=i0
     ! k+q-point set
     ikq=ikmapikq(iq,ik)     
     call getoccsv(vkl(1,ikq),occsv(1,ikq))
     do i=1,nstsv
        if (occsv(i,ik).lt.(1.d0-epsocc)) exit
     end do
     io(ik)=i-1
     do i=nstsv,1,-1
        if (occsv(i,ik).gt.epsocc) exit
     end do
     iu(ik)=i
  end do
  ! overall highest (partially) occupied state
  iocc0=maxval(io0)
  iocc=maxval(io)
  ! overall lowest (partially) unoccupied state
  iunocc0=maxval(iu0)
  iunocc=maxval(iu)

  ! debug output
  if (dbglev.gt.0) then
     write(*,*) 'Debug(findocclims):'
     write(*,*) ' iocc0',iocc0
     write(*,*) ' iocc',iocc
     write(*,*) ' iunocc0',iunocc0
     write(*,*) ' iunocc',iunocc
     write(*,*) ' io0',io0
     write(*,*) ' io',io
     write(*,*) ' iu0',iu0
     write(*,*) ' iu',iu
     write(*,*)
  end if

end subroutine findocclims
