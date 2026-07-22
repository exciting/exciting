! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: getevalsv
! !INTERFACE:
!
subroutine getevalsv(vpl, evalsvp)
! !USES:
  use modmain
  use modinput
  use modmpi
  use mod_large_io, only: inquire_large, open_direct_unformatted_large
  use precision, only: i32, dp, long_int, str_256
! !DESCRIPTION:
!   The file where the (second-variational) eigenvalues are stored is
!    {\tt EVALSV.OUT}.
!   It is a direct-access binary file, the record length of which can be
!   determined
!   with the help of the array sizes and data type information.
!   One record of this file has the following structure
!
!   \begin{tabular}{|l|l|l|}
!   \hline
!   $k_{\rm lat}$ & $N_{\rm stsv}$ & $E$ \\
!   \hline
!   \end{tabular}\newline\newline
!   The following table explains the parts of the record in more detail
!
!   \begin{tabular}{|l|l|l|l|}
!   \hline
!   name & type & shape & description\\
!   \hline \hline
!   $k_{\rm lat}$ & real(8) & 3 & k-point in lattice coordinates \\ \hline
!   $N_{\rm stsv}$ & integer & 1 & number of (second-variational) states \\
!    &  &  & (without core states) \\ \hline
!   $E$ & real(8) & $N_{\rm stsv}$ & (second-variational) eigenvalue array \\
!   \hline
!   \end{tabular}\newline\newline
!
! !REVISION HISTORY:
!   Documentation added, Dec 2009 (S. Sagmeister)
!   Consmetic changes, added more comments. 2016 (Aurich)
!EOP
!BOC

  implicit none

  ! Arguments
  real(dp), intent(in) :: vpl(3)
  real(dp), intent(out) :: evalsvp(nstsv)

  ! Local variables
  logical :: lexist
  integer(i32) :: io_unit, isym, ik, koffset, i, nstsv_
  integer(long_int) :: reclen
  real(dp) :: vkl_(3), t1
  character(str_256) :: filetag

  ! External functions
  character(str_256), external :: outfilenamestring


  ! Added feature to access arrays for only a subset of bands
  real(dp), allocatable :: evalsv_(:)

  ! Find the k-point number
  call findkpt(vpl, isym, ik)

  ! Find the record length
  call inquire_large( reclen, vkl_, [nstsv_] )

  ! mod_names:filetag_evalsv is 'EVALSV'
  filetag = trim(filetag_evalsv)
!write(*,*) "getevalsv: Reading from ", outfilenamestring(filetag, ik) 

  ! Try to open 'EVALSV(<krange>)<modmisc:filext>' 
  ! (krange applies only for specific task in mpi mode)
  do i = 1, 100
    inquire(file=outfilenamestring(filetag, ik), exist=lexist)
    if(lexist) then
      call open_direct_unformatted_large( io_unit, outfilenamestring(filetag, ik), "read", reclen, "old" )
      exit
    else
      call system('sync')
      write(*,*) "Waiting for other process to write"&
        & // ":getevalsv:" // trim(outfilenamestring(filetag, ik))
      call sleep(5)
    end if
  end do

  if(splittfile) then
     koffset = ik - firstofset(procofindex(ik, nkpt), nkpt) + 1
  else
     koffset = ik
  end if

  ! Get dimensions of stored file
  read(io_unit, rec=1) vkl_, nstsv_
  close(io_unit)

  if(nstsv .gt. nstsv_) then
     write(*,*)
     write(*, '("Error(getevalsv): invalid nstsv for k-point ",i8)') ik
     write(*, '(" current    : ",i8)') nstsv
     write(*, '(" evalsv.out : ",i8)') nstsv_
     write(*, '(" file       : ",a      )') trim(outfilenamestring(filetag, ik))
     write(*,*)
     stop
  end if

  allocate(evalsv_(nstsv_))

  call inquire_large( reclen, vkl_, [nstsv_], evalsv_ )

  call open_direct_unformatted_large( io_unit, outfilenamestring(filetag, ik), "read", reclen, "old" )

  read(io_unit, rec=koffset) vkl_, nstsv_, evalsv_

  ! Retrieve subset
  evalsvp(:) = evalsv_ (:nstsv)
  deallocate(evalsv_)

  close(io_unit)

  t1 = abs(vkl(1, ik)-vkl_(1)) + abs(vkl(2, ik)-vkl_(2)) + abs(vkl(3, ik)-vkl_(3))

  if(t1 .gt. input%structure%epslat) then
    write(*,*)
    write(*, '("Error(getevalsv): differing vectors for k-point ",i8)') ik
    write(*, '(" current    : ",3g18.10)') vkl(:, ik)
    write(*, '(" evalsv.out : ",3g18.10)') vkl_
    write(*, '(" file       : ",a      )') trim(outfilenamestring(filetag, ik))
    write(*,*)
    stop
  end if

end subroutine getevalsv
!EOC

module m_getevalsvr

      implicit none

  contains

    subroutine getevalsvr(fname, isti, istf, vpl, evalsvp)
       use modmain
       use precision, only: i32, dp, str_256

       implicit none

       ! Arguments
       character(*), intent(in) :: fname
       integer(i32), intent(in) :: isti, istf
       real(dp), intent(in) :: vpl(3)
       real(dp), intent(out) :: evalsvp(:)

       ! Local variables
       integer(i32) :: chkerr
       real(dp), allocatable :: evalsvt(:)
       character(str_256) :: tmpstr

       ! Check correct shapes
       chkerr = 0

       if((isti .lt. 1) .or. (istf .gt. nstsv) .or. (istf .le. isti)) then
         write(*,*)
         write(*, '("Error(getevalsvr): inconsistent limits for bands:")')
         write(*, '(" band limits  : ",2i6)') isti, istf
         write(*, '(" maximum value: ",i6)') nstsv
         write(*,*)
         chkerr = chkerr + 1
       end if

       if(size(evalsvp, 1) .ne. (istf-isti+1)) then
         write(*,*)
         write(*, '("Error(getevalsvr): output array does not match for bands:")')
         write(*, '(" band limits              : ",2i6)') isti, istf
         write(*, '(" requested number of bands: ",i6)') istf - isti + 1
         write(*, '(" array size               : ",i6)') size(evalsvp, 1)
         write(*,*)
         chkerr = chkerr + 1
       end if

       if(chkerr .ne. 0) stop

       allocate(evalsvt(nstsv))

       filetag_evalsv = trim(fname)

       tmpstr = trim(filext)
       filext = ''

       call getevalsv(vpl, evalsvt)

       filetag_evalsv = 'EVALSV'
       filext = trim(tmpstr)

       evalsvp(:) = evalsvt(isti:istf)

       deallocate(evalsvt)
    end subroutine getevalsvr
end module m_getevalsvr
