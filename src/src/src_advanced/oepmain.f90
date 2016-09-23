!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
!
!BOP
! !ROUTINE: oepmain
! !INTERFACE:
Subroutine oepmain
! !USES:
      Use modmain
      Use modinput
      use modmpi
! !DESCRIPTION:
!   Main routine for the calculation of the optimized effective potential.
!
! !REVISION HISTORY:
!  Created ()
!  Modified Feb 2014 (UW)
!  Modified March 2014 (UW)
!EOP
!BOC

      Implicit None
! local variables
      Integer :: is, ia, ias, ik
      Integer :: ir, irc, it, idm
      integer :: msize
      Real (8) :: tau, resp, t1 
! allocatable arrays
      Real (8), Allocatable :: rflm (:)
      Real (8), Allocatable :: rfmt (:, :, :)
      Real (8), Allocatable :: rfir (:)
      Real (8), Allocatable :: rvfmt (:, :, :, :)
      Real (8), Allocatable :: rvfir (:, :)
      Real (8), Allocatable :: dvxmt(:,:,:)
      Real (8), Allocatable :: dvxir(:)
      Real (8), Allocatable :: dbxmt(:,:,:,:)
      Real (8), Allocatable :: dbxir(:,:)
      Complex (8), Allocatable :: vnlcv (:, :, :, :)
      Complex (8), Allocatable :: vnlvv (:, :, :)
      Complex (8), Allocatable :: zflm (:)
! external functions
      Real (8) :: rfinp
      Complex (8) zfint
      External rfinp, zfint
      Real (8) :: t_0,t_1
      If (iscl .Lt. 1) Return
! calculate nonlocal matrix elements
      Allocate (vnlcv(ncrmax, natmtot, nstsv, nkpt))
      Allocate (vnlvv(nstsv, nstsv, nkpt))
      Call timesec(t_0) 
      Call oepvnl (vnlcv, vnlvv)
      Call timesec(t_1)
      time_oepvnl=time_oepvnl+(t_1-t_0)
! allocate local arrays
      Allocate (rflm(lmmaxvr))
      Allocate (rfmt(lmmaxvr, nrmtmax, natmtot))
      Allocate (rfir(ngrtot))
      Allocate (dvxmt(lmmaxvr, nrcmtmax, natmtot))
      Allocate (dvxir(ngrtot))
      Allocate (zflm(lmmaxvr))
      If (associated(input%groundstate%spin)) Then
         Allocate (rvfmt(lmmaxvr, nrmtmax, natmtot, ndmag))
         Allocate (rvfir(ngrtot, ndmag))
         Allocate (dbxmt(lmmaxvr, nrcmtmax, natmtot, ndmag))
         Allocate (dbxir(ngrtot, ndmag))
      End If
! zero the complex potential
      zvxmt (:, :, :) = 0.d0
      zvxir (:) = 0.d0
      If (associated(input%groundstate%spin)) Then
         zbxmt (:, :, :, :) = 0.d0
         zbxir (:, :) = 0.d0
      End If
      resp = 0.d0
! initial step size
      tau = input%groundstate%OEP%tauoep(1)
! start iteration loop
      Call timesec(t_0)
      Do it = 1, input%groundstate%OEP%maxitoep
         If ((Mod(it, 50) .Eq. 0).and.(rank==0)) Then
            Write (*, '("Info(oepmain): done ", I4, " iterations of ", I4)') &
           &  it, input%groundstate%OEP%maxitoep
         End If
! zero the residual
         dvxmt(:,:,:) = 0.d0
         dvxir(:) = 0.d0
         If (associated(input%groundstate%spin)) Then
            dbxmt(:,:,:,:) = 0.d0
            dbxir(:,:) = 0.d0
         End If
! calculate the k-dependent residuals
#ifdef MPI
         Do ik = firstk (rank), lastk (rank)
#else
         Do ik = 1, nkpt
#endif            
             Call oepresk(ik, vnlcv, vnlvv, dvxmt, dvxir, dbxmt, dbxir)
         End Do
#ifdef MPI
         msize = lmmaxvr*nrcmtmax*natmtot
         call MPI_ALLREDUCE(MPI_IN_PLACE, dvxmt, msize, MPI_DOUBLE_PRECISION, &
        &  MPI_SUM, MPI_COMM_WORLD, ierr)
         msize = ngrtot
         call MPI_ALLREDUCE(MPI_IN_PLACE, dvxir, msize, MPI_DOUBLE_PRECISION, &
        &  MPI_SUM, MPI_COMM_WORLD, ierr)
         if (associated(input%groundstate%spin)) Then
            msize = lmmaxvr*nrcmtmax*natmtot*ndmag
            call MPI_ALLREDUCE(MPI_IN_PLACE, dbxmt, msize, MPI_DOUBLE_PRECISION, &
           &  MPI_SUM, MPI_COMM_WORLD, ierr)
            msize = ngrtot*ndmag
            call MPI_ALLREDUCE(MPI_IN_PLACE, dbxir, msize, MPI_DOUBLE_PRECISION, &
           &  MPI_SUM, MPI_COMM_WORLD, ierr)
         end if
#endif
! convert muffin-tin residuals to spherical harmonics
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
               irc = 0
               Do ir = 1, nrmt (is), input%groundstate%lradstep
                  irc = irc + 1
                  Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rfshtvr, &
                 & lmmaxvr, dvxmt(:, irc, ias), 1, 0.d0, rfmt(:, ir, &
                 & ias), 1)
                  Do idm = 1, ndmag
                     Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rfshtvr, &
                    & lmmaxvr, dbxmt(:, irc, ias, idm), 1, 0.d0, &
                    & rvfmt(:, ir, ias, idm), 1)
                  End Do
               End Do
            End Do
         End Do
! symmetrise the residuals
         Call symrf (input%groundstate%lradstep, rfmt, dvxir)
         If (associated(input%groundstate%spin)) Call symrvf &
        & (input%groundstate%lradstep, rvfmt, dbxir)
! magnitude of residuals
         resoep = Sqrt (Abs(rfinp(input%groundstate%lradstep, rfmt, &
        & rfmt, dvxir, dvxir)))
         Do idm = 1, ndmag
            t1 = rfinp (input%groundstate%lradstep, rvfmt(:, :, :, &
           & idm), rvfmt(:, :, :, idm), dbxir(:, idm), dbxir(:, idm))
            resoep = resoep + Sqrt (Abs(t1))
         End Do
         resoep = resoep / omega
! in case convergencece is reached go to end of iteration loop
         If (abs(resoep).le.input%groundstate%OEP%convoep) exit
! adjust step size
         If (it .Gt. 1) Then
            If (resoep .Gt. resp) Then
               tau = tau * input%groundstate%OEP%tauoep(2)
            Else
               tau = tau * input%groundstate%OEP%tauoep(3)
            End If
         End If
         resp = resoep
!--------------------------------------------!
!     update complex potential and field     !
!--------------------------------------------!
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
               irc = 0
               Do ir = 1, nrmt (is), input%groundstate%lradstep
                  irc = irc + 1
! convert residual to spherical coordinates and subtract from complex potential
                  Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rbshtvr, &
                 &  lmmaxvr, rfmt(:, ir, ias), 1, 0.d0, rflm, 1)
                  zvxmt (:, irc, ias) = zvxmt (:, irc, ias) - tau * rflm (:)
                  Do idm = 1, ndmag
                     Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rbshtvr, &
                    & lmmaxvr, rvfmt(:, ir, ias, idm), 1, 0.d0, rflm, 1)
                     zbxmt (:, irc, ias, idm) = zbxmt (:, irc, ias, idm) - tau * rflm (:)
                  End Do
               End Do
            End Do
         End Do
         zvxir (:) = zvxir (:) - tau * dvxir (:)
         Do idm = 1, ndmag
            zbxir (:, idm) = zbxir (:, idm) - tau * dbxir (:, idm)
         End Do
! end iteration loop
      End Do
      Call timesec(t_1)
      time_oep_iter=time_oep_iter+(t_1-t_0)
! generate the real potential and field
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            irc = 0
            Do ir = 1, nrmt (is), input%groundstate%lradstep
               irc = irc + 1
! convert to real spherical harmonics
               rflm (:) = dble (zvxmt(:, irc, ias))
               Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rfshtvr, &
              & lmmaxvr, rflm, 1, 0.d0, rfmt(:, ir, ias), 1)
               Do idm = 1, ndmag
                  rflm (:) = dble (zbxmt(:, irc, ias, idm))
                  Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rfshtvr, &
                 & lmmaxvr, rflm, 1, 0.d0, rvfmt(:, ir, ias, idm), 1)
               End Do
            End Do
         End Do
      End Do
! convert potential and field from a coarse to a fine radial mesh
      Call rfmtctof (rfmt)
      Do idm = 1, ndmag
         Call rfmtctof (rvfmt(:, :, :, idm))
      End Do
! add to existing correlation potential and field
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ir = 1, nrmt (is)
               vxcmt (:, ir, ias) = vxcmt (:, ir, ias) + ex_coef*rfmt (:, ir, ias)
               ! check if ex_coef needed also here for hybrids
               Do idm = 1, ndmag
                  bxcmt (:, ir, ias, idm) = bxcmt (:, ir, ias, idm) + &
                 & rvfmt (:, ir, ias, idm)
               End Do
            End Do
         End Do
      End Do
      vxcir (:) = vxcir (:) + ex_coef* (dble (zvxir(:)))
      Do idm = 1, ndmag
! check if ex_coef also needed here for hybrids
         bxcir (:, idm) = bxcir (:, idm) + dble (zbxir(:, idm))
      End Do
! symmetrise the exchange potential and field
      Call symrf (1, vxcmt, vxcir)
      If (associated(input%groundstate%spin)) Then
         Call symrvf (1, bxcmt, bxcir)
      End If
      Deallocate (rflm, rfmt, rfir, vnlcv, vnlvv)
      Deallocate (dvxmt, dvxir, zflm)
      If (associated(input%groundstate%spin)) Then
         Deallocate (rvfmt, rvfir)
         Deallocate (dbxmt, dbxir)
      End If

!********************************************
! Calculation of the potential discontinuity
!         after last iteration
!********************************************
      if (tlast) call deltax     

      Return
End Subroutine
