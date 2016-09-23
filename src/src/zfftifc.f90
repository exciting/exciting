!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: zfftifc
! !INTERFACE:
!
!
Subroutine zfftifc (nd, n, sgn, z)
! !INPUT/OUTPUT PARAMETERS:
!   nd   : number of dimensions (in,integer)
!   n    : grid sizes (in,integer(nd))
!   sgn  : FFT direction, -1: forward; 1: backward (in,integer)
!   z    : array to transform (inout,complex(n(1)*n(2)*...*n(nd)))
! !DESCRIPTION:
!   Interface to the double-precision complex fast Fourier transform routine.
!   This is to allow machine-optimised routines to be used without affecting the
!   rest of the code. See routine {\tt nfftifc}.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: nd
      Integer, Intent (In) :: sgn
      Integer, Intent (In) :: n (nd)
      Complex (8), Intent (Inout) :: z (*)
!
!-------------------------------------!
!     interface to HP MLIB Z3DFFT     !
!-------------------------------------!
!integer ier
!if (nd.eq.3) then
!  call z3dfft(z,n(1),n(2),n(3),n(1),n(2),sgn,ier)
!end if
!
!-------------------------------------!
!     interface to FFTW version 3     !
!-------------------------------------!

#ifdef FFTW
integer, parameter :: FFTW_ESTIMATE=64
integer i,p
integer(8) plan
real(8) t1
real(8) ta,tb
!!$OMP CRITICAL
!call timesec(ta)
call dfftw_plan_dft(plan,nd,n,z,z,sgn,FFTW_ESTIMATE)
!call timesec(tb)
!write(*,*) 'fftw',tb-ta
!!$OMP END CRITICAL
!call timesec(ta)
call dfftw_execute(plan)
!call timesec(tb)
!write(*,*) 'fftw',tb-ta
!!$OMP CRITICAL
!call timesec(ta)
call dfftw_destroy_plan(plan)
!call timesec(tb)
!write(*,*) 'fftw',tb-ta
!read(*,*)
!!$OMP END CRITICAL
if (sgn.eq.-1) then
  p=1
  do i=1,nd
    p=p*n(i)
  end do
  t1=1.d0/dble(p)
  call zdscal(p,t1,z,1)
end if
#else
!
!----------------------------------!
!     interface to MKL 8.1/9.1     !
!----------------------------------!
! (with thanks to Torbjorn Bjorkman)
!use MKL_DFTI ! this module required by MKL
!integer dftistatus,i,p
!real(8) dftiscale
!type(DFTI_DESCRIPTOR), POINTER :: handle
!p=1
!do i=1,nd
!  p=p*n(i)
!end do
!dftiscale=1.d0/dble(p)
!dftistatus=DftiCreateDescriptor(handle,DFTI_DOUBLE,DFTI_COMPLEX,nd,n)
!dftistatus=DftiSetValue(handle, DFTI_FORWARD_SCALE,dftiscale)
!dftistatus=DftiCommitDescriptor(handle)
!if (sgn.eq.-1) then
!  dftistatus=DftiComputeForward(handle,z)
!else
!  dftistatus=DftiComputeBackward(handle,z)
!end if
!dftistatus=DftiFreeDescriptor(handle)
!
!----------------------------------------!
!     interface to modified FFTPACK5     !
!----------------------------------------!
      Call cfftnd (nd, n, sgn, z)
#endif
!
      Return
End Subroutine
!EOC
