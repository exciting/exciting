
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

Module modmixermsec

      use precision, only: i32, long_int, dp, str_256

      real(dp), allocatable :: residual(:), last_outputp(:), work2(:), work3(:)
      real(dp), allocatable :: PWHIST(:), FHIST(:), CLMHIST(:), yhist(:)
      integer(i32) :: noldstepsin_file, noldsteps, MUSE, IDSCALE, noldstepsmax
      integer(long_int) :: record_of_last_iter
      integer(i32), parameter :: icond = 1
      real(dp), parameter :: dbase = 0.005_dp
      real(dp) :: scl_plane, qmx, RedOld, RedPred, qmx_input, PM1, DIAG, dmix_last, dmixout(4)
      real(dp) :: MSECINFO(20), rtrap, SCHARGE, TCharge, splane, tplane, qtot
      real(dp) :: dmix

Contains

      Subroutine initmixermsec (n,nmax)
         use modmain, only: CHGIR, CHGMTTOT, input
         implicit none
         integer(long_int), Intent (In) :: n
         integer(i32), Intent (In) :: nmax
         integer(i32) :: niter
         noldstepsmax=nmax
!         noldstepsmax=input%groundstate%msecStoredSteps
         if (allocated(residual)) deallocate(residual)
         allocate(residual(n))
         if (allocated(last_outputp)) deallocate(last_outputp)
         allocate(last_outputp(n))
         if (allocated(work2)) deallocate(work2)
         allocate(work2(n))
         if (allocated(work3)) deallocate(work3)
         allocate(work3(n))
         if (allocated(PWHIST)) deallocate(PWHIST)
         allocate(PWHIST(noldstepsmax))
         if (allocated(FHIST)) deallocate(FHIST)
         allocate(FHIST(noldstepsmax))
         if (allocated(CLMHIST)) deallocate(CLMHIST)
         allocate(CLMHIST(noldstepsmax))
         if (allocated(yhist)) deallocate(yhist)
         allocate(yhist(noldstepsmax))
         record_of_last_iter = 0
         residual = 0
         last_outputp = 0
         work2 = 0
         work3 = 0
         PWHIST = 0
         FHIST = 0
         CLMHIST = 0
         yhist = 0
         scl_plane = 4
         RedOld = 1
         RedPred = 1
         qmx_input = 0.2_dp
         qmx = qmx_input
         PM1 = 1
         IDSCALE = 1
         DIAG = 5e-4_dp
         noldstepsin_file = 0
         noldsteps = 0
         rtrap = 0.1_dp
         SCHARGE = CHGIR
         TCharge = CHGMTTOT
         splane = 0.000001_dp
         tplane = 0.000001_dp
         MSECINFO = 1
         dmix = 0.5_dp
         dmix_last=0.5_dp
      End Subroutine


      Subroutine freearraysmixermsec ()
         use mod_misc, only: scrpath
         Character (str_256), External :: outfilenamestring
         Character (str_256) :: filetag
         filetag = "BROYDEN"
         if (allocated(residual)) deallocate(residual)
         if (allocated(last_outputp)) deallocate(last_outputp)
         If (allocated(work2)) deallocate (work2)
         If (allocated(work3)) deallocate (work3)
         if (allocated(PWHIST)) deallocate(PWHIST)
         if (allocated(FHIST)) deallocate(FHIST)
         if (allocated(CLMHIST)) deallocate(CLMHIST)
         if (allocated(yhist)) deallocate(yhist)
         Open (23, File=trim(scrpath)//trim(filetag)//'.OUT')
         Close (23, Status='DELETE')
      End Subroutine

End Module
