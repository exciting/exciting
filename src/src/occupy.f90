!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: occupy
! !INTERFACE:
!
!
Subroutine occupy
! !USES:
      Use modinput
      Use modmain
      use mod_opt_tetra

! !DESCRIPTION:
!   Finds the Fermi energy and sets the occupation numbers for the
!   second-variational states using the routine {\tt fermi}.
!
! !REVISION HISTORY:
!   Created February 2004 (JKD)
!   Modifiactions for tetrahedron method, November 2007 (RGA alias
!     Ricardo Gomez-Abal)
!   Modifications for tetrahedron method, 2007-2010 (Sagmeister)
!   Modifications for tetrahedron method, 2011 (DIN)
!   Simplicistic method for systems with gap added, 2013 (STK)
!EOP
!BOC
      Implicit None
! local variables
      Integer, Parameter :: maxit = 1000
      real(8), parameter :: de0=1.d0
      Integer :: ik, ist, it, nvm
      Real (8) :: e0, e1, chg, x, t1
! external functions
      Real (8) :: sdelta, stheta
      real(8) :: egap
      real(8) :: dfde(nstsv,nkpt)
      
      character(1024) :: message
      
      external sdelta, stheta
      real(8), external :: dostet_exciting
      
      if ( input%groundstate%stypenumber .ge. 0 ) then
         t1 = 1.d0 / input%groundstate%swidth

!     next lines taken in part from libbzint (STK)
!
!!    nvm is the number of bands for an insulating system 
!!    since for a system with gap, the procedure to determine the
!!    band gap can be unstable, just try first whether it is an
!!    insulating system, but such a simplicistic way to determine the Fermi energy
!!    is valid only for no spin polarized cases 
!
         if (.not.associated(input%groundstate%spin)) then
           nvm  = nint(chgval/occmax)
           e0 = maxval(evalsv(nvm,:))
           e1 = minval(evalsv(nvm+1,:))
           efermi = 0.5*(e0 + e1)

           fermidos = 0.d0
           chg = 0.d0
           Do ik = 1, nkpt
              Do ist = 1, nstsv
                 x = (evalsv(ist, ik)-efermi) * t1
                 fermidos = fermidos + wkpt (ik) * sdelta (input%groundstate%stypenumber, x) * t1
                 occsv (ist, ik) = occmax * stheta (input%groundstate%stypenumber, -x)
                 chg = chg + wkpt (ik) * occsv (ist, ik)
              End Do
           End Do
           fermidos = fermidos * occmax
           if ((e1 .ge. e0) .and. (abs(chg - chgval) .lt. input%groundstate%epsocc)) then
!            Write (*, '("Info(occupy): System has gap, simplicistic method used in determining efermi and occupation")')
             goto 10
           endif
         end if

! find minimum and maximum eigenvalues
         e0 = evalsv (1, 1)
         e1 = e0
         Do ik = 1, nkpt
            Do ist = 1, nstsv
               e0 = Min (e0, evalsv(ist, ik))
               e1 = Max (e1, evalsv(ist, ik))
            End Do
         End Do
!	 	 
! determine the Fermi energy using the bisection method
!
         Do it = 1, maxit
            efermi = 0.5d0 * (e0+e1)
            chg = 0.d0
            Do ik = 1, nkpt
               Do ist = 1, nstsv
                  x = (efermi-evalsv(ist, ik)) * t1
                  occsv (ist, ik) = occmax * stheta &
                   & (input%groundstate%stypenumber, x)
                  chg = chg + wkpt (ik) * occsv (ist, ik)
               End Do
            End Do
            If (chg .Lt. chgval) Then
               e0 = efermi
            Else
               e1 = efermi
            End If
            If ((e1-e0) .Lt. input%groundstate%epsocc) Go To 10
         End Do
         Write (*,*)
         Write (*, '("Error(occupy): could not find Fermi energy")')
         Write (*,*)
         Stop
10       Continue
! find the density of states at the Fermi surface in units of
! states/Hartree/spin/unit cell
         fermidos = 0.d0
         Do ik = 1, nkpt
            Do ist = 1, nstsv
               x = (evalsv(ist, ik)-efermi) * t1
               fermidos = fermidos + wkpt (ik) * sdelta &
                 & (input%groundstate%stypenumber, x) * t1
            End Do
            If (occsv(nstsv, ik) .Gt. input%groundstate%epsocc) Then
               call warning('Warning(occupy):')
               Write (message, '(" Not enough empty states for k-point ", I6)') ik
               call warning(message)
            End If
         End Do
         fermidos = fermidos * occmax

      else  if (input%groundstate%stypenumber==-1) then
         !------------------------------------------------------------
         ! Use the tetrahedron integration method (LIBBZINT library)
         !------------------------------------------------------------

         ! Calculate the Fermi energy
         call fermi_exciting(input%groundstate%tevecsv, &
         &                   chgval, &
         &                   nstsv,nkpt,evalsv, &
         &                   ntet,tnodes,wtet,tvol, &
         &                   efermi,egap,fermidos)
         
         ! Calculate state occupation numbers
         call tetiw(nkpt,ntet,nstsv,evalsv,tnodes,wtet,tvol,efermi,occsv)
         do ik = 1, nkpt
           do ist = 1, nstsv
             occsv(ist,ik) = occmax/wkpt(ik)*occsv(ist,ik)
           end do
         end do
         
      else  if (input%groundstate%stypenumber==-2) then
         !--------------------------------------
         ! Use the improved tetrahedron method
         !--------------------------------------
         call opt_tetra_efermi(chgval/dble(occmax),nkpt,nstsv,evalsv,efermi,occsv)
         do ik = 1, nkpt
           occsv(:,ik) = dble(occmax)/wkpt(ik)*occsv(:,ik)
         end do
         !write(*,*) 'occsv=', occsv(:,1)
      
         call opt_tetra_dos(nkpt,nstsv,evalsv,efermi,dfde)
         fermidos = sum(dfde)
         !write(*,*) 'dos at Ef=', fermidos
         
      End If ! modified tetrahedron integration method
      
      Return
End Subroutine
!EOC
