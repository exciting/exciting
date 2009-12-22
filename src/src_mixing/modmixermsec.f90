!
!
Module modmixermsec
!
      Real (8), Allocatable :: residual (:), last_outputp (:), work2 &
     & (:), work3 (:)
      Real (8), Allocatable :: PWHIST (:), FHIST (:), CLMHIST (:), &
     & yhist (:)
      Integer :: record_of_last_iter, noldstepsin_file, noldsteps, &
     & MUSE, IDSCALE
      Integer, Parameter :: icond = 1, noldstepsmax = 8, dbase = &
     & 0.005D0
      Real (8) :: scl_plane, qmx, RedOld, RedPred, qmx_input, PM1, &
     & DIAG, dmix_last, dmixout (4)
      Real (8) :: MSECINFO (20), rtrap, SCHARGE, TCharge, splane, &
     & tplane, qtot
      Real (8) :: dmix
Contains
!
!
      Subroutine initmixermsec (n)
         Use modmain, Only: CHGIR, CHGMTTOT
         Integer, Intent (In) :: n
!
         Integer :: niter
         Allocate (residual(n), last_outputp(n), work2(n), work3(n))
!
         Allocate (PWHIST(noldstepsmax), FHIST(noldstepsmax), &
        & CLMHIST(noldstepsmax), yhist(noldstepsmax))
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
         qmx_input = .2
         qmx = qmx_input
         PM1 = 1
         IDSCALE = 1
         DIAG = 5D-4
         noldstepsin_file = 0
         noldsteps = 0
         rtrap = 0.1
         SCHARGE = CHGIR
         TCharge = CHGMTTOT
         splane = .000001
         tplane = .000001
         MSECINFO = 1
!
         dmix = .5
      End Subroutine
!
!
      Subroutine freearraysmixermsec ()
         Character (256), External :: outfilenamestring
         Character (256) :: filetag
         filetag = "BROYDEN"
         Deallocate (residual, last_outputp)
         If (allocated(work2)) deallocate (work2)
         If (allocated(work3)) deallocate (work3)
         Deallocate (PWHIST, FHIST, CLMHIST, yhist)
         Open (23, File=outfilenamestring(filetag, 1))
         Close (23, Status='DELETE')
      End Subroutine
!
End Module
